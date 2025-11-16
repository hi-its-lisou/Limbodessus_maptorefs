#!/usr/bin/env bash
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=24:00:00
#SBATCH --mem=10G
#SBATCH --cpus-per-task=10
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=elisabethwilliamsonn@gmail.com

set -euo pipefail # exit on first error


########################################
# CONFIG: inputs and outpots
########################################
REF="/hpcfs/users/a1707473/allgenes_noUTR_concat.fasta" #reference to map against
IN_DIR="/hpcfs/users/a1707473/raw_reads/" # FASTQs
THREADS="${SLURM_CPUS_PER_TASK:-4}" 
REMOVE_DUPLICATES="true" #Picard behaviour, see https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates

OUT_ROOT="/hpcfs/users/a1707473/map_to_refs" 
PEAR_DIR="${OUT_ROOT}/pear"
ALN_DIR="${OUT_ROOT}/aligned"
MERGE_DIR="${OUT_ROOT}/merged"
TMPDIR="${OUT_ROOT}/tmp"

mkdir -p "$PEAR_DIR" "$ALN_DIR" "$MERGE_DIR" "$TMPDIR"

########################################
# MODULES
########################################
# Needed to include private module path because I installed PEAR locally.
module use /hpcfs/users/a1707473/.local/easybuild/modules/all
module use /home/a1707473/.local/easybuild/modules/all

module load BWA/0.7.17-GCCcore-11.2.0
module load SAMtools/1.17-GCC-11.2.0
module load PEAR/0.9.11-GCC-11.2.0
module load Java/17.0.6
module load picard/3.3.0-Java-17.0.6

### The two steps below are only neccessary running on pheonix
# Resolve Picard command
if command -v picard >/dev/null 2>&1; then
  PICARD_CMD="picard"
elif [[ -n "${EBROOTPICARD:-}" && -f "$EBROOTPICARD/picard.jar" ]]; then
  PICARD_CMD="java -jar $EBROOTPICARD/picard.jar"
else
  echo "ERROR: Picard not found (no 'picard' and no \$EBROOTPICARD/picard.jar)" >&2
  exit 2
fi

# Help Java not fill shared /tmp by forcing picard to use my own tmp dir 
export _JAVA_OPTIONS="-Djava.io.tmpdir=$TMPDIR"

module list
echo "Time is $(date)"
echo "Using PICARD_CMD: $PICARD_CMD"
echo

########################################
# Debugging checks
########################################
# Checks for reference 
if [[ ! -f "$REF" ]]; then
  echo "Reference FASTA not found: $REF" >&2
  exit 1
fi

# Index the reference if needed
if [[ ! -f "${REF}.bwt" ]]; then
  echo "[INFO] Indexing reference..."
  bwa index "$REF"
fi

########################################
# HELPERS: wrap BWA and SAM
########################################

#Align single-end (assmebled/merged) reads, sort them, and index

align_single() {
  local fq="$1" sample="$2"
  local bam="${ALN_DIR}/${sample}.assembled.SE.sorted.bam"
  local RG="@RG\\tID:${sample}\\tSM:${sample}\\tPL:ILLUMINA\\tLB:${sample}\\tPU:${sample}"

#Map with BWA-MEM, convert back to BAM, sort, index.
  bwa mem -t "$THREADS" -B 3 -L 10,10 -R "$RG" "$REF" "$fq" \
    | samtools view -@ "$THREADS" -b - \
    | samtools sort -@ "$THREADS" -o "$bam" -
  samtools index -@ "$THREADS" "$bam"
  echo "$bam"
}

# Align paired-end (unmerged) reads, sort, index.
align_paired() {
  local fq1="$1" fq2="$2" sample="$3"
  local bam="${ALN_DIR}/${sample}.unassembled.PE.sorted.bam"
  local RG="@RG\\tID:${sample}\\tSM:${sample}\\tPL:ILLUMINA\\tLB:${sample}\\tPU:${sample}"

#Map pairs with BWA-MEM, convert to BAM, sort, index.
  bwa mem -t "$THREADS" -B 3 -L 10,10 -R "$RG" "$REF" "$fq1" "$fq2" \
    | samtools view -@ "$THREADS" -b - \
    | samtools sort -@ "$THREADS" -o "$bam" -
  samtools index -@ "$THREADS" "$bam"
  echo "$bam"
}

########################################
# FIND INPUT PAIRS (_1/_2, fastq/f/gz or not)
########################################
shopt -s nullglob

#Collect and list all candidate R1 files
readarray -t R1_LIST < <(printf '%s\n' \
  "${IN_DIR}"/*_1*.fastq \
  "${IN_DIR}"/*_1*.fq \
  "${IN_DIR}"/*_1*.fastq.gz \
  "${IN_DIR}"/*_1*.fq.gz 2>/dev/null | sort -u)

echo "Scanning inputs in: $IN_DIR"
printf 'R1 candidates (%d):\n' "${#R1_LIST[@]}"
printf '  %s\n' "${R1_LIST[@]:-<none>}"
echo

# Bail out early if nothing matched the expected patterns
if (( ${#R1_LIST[@]} == 0 )); then
  echo "ERROR: No FASTQs matched *_1*.fastq/.fq[.gz] in $IN_DIR" >&2
  exit 3
fi

########################################
# MAIN: process discovered samples
########################################
for R1 in "${R1_LIST[@]}"; do
  # Derive R2 path by swapping the first _1 → _2 in the filename
  R2="${R1/_1/_2}"
  if [[ ! -f "$R2" ]]; then
    echo ">> Skipping: mate not found for $R1 (expected $R2)" >&2
    continue
  fi

# Extract sample name as the prefix before the first _1 occurrence
  base="$(basename "$R1")"
  sample="${base%%_1*}"

# Progress banner for the current sample
  echo "========== SAMPLE: ${sample} =========="
  echo "R1: $R1"
  echo "R2: $R2"

  # 1) PEAR: collapse overlaps (accepts gz or plain fastq)
  PFX="${PEAR_DIR}/${sample}"
  pear -f "$R1" -r "$R2" -o "$PFX" -j "$THREADS"

# Define expected PEAR outputs for merged and unmerged reads
  FQ_MERGED="${PFX}.assembled.fastq"
  FQ_UNM_F="${PFX}.unassembled.forward.fastq"
  FQ_UNM_R="${PFX}.unassembled.reverse.fastq"

  # 2) Align SE (assembled) and PE (unassembled) separately
  mapped_bams=()

# Align merged single-end reads
  if [[ -s "$FQ_MERGED" ]]; then
    bam_se="$(align_single "$FQ_MERGED" "$sample")"
    mapped_bams+=("$bam_se")
  else
    echo ">> No assembled (SE) reads for ${sample}"
  fi

# Align unmerged paired-end reads
  if [[ -s "$FQ_UNM_F" && -s "$FQ_UNM_R" ]]; then
    bam_pe="$(align_paired "$FQ_UNM_F" "$FQ_UNM_R" "$sample")"
    mapped_bams+=("$bam_pe")
  else
    echo ">> No unassembled (PE) reads for ${sample}"
  fi

# If nothing mapped for this sample, skip further processing
  if (( "${#mapped_bams[@]}" == 0 )); then
    echo "!!! No BAMs produced for ${sample}; skipping"
    continue
  fi

# 3) Merge the sorted BAMs, re-sort, then MarkDuplicates (index pre-dedup)
  merged_pre="${MERGE_DIR}/${sample}.preDedup.merged.bam"
  if (( "${#mapped_bams[@]}" == 1 )); then
    # Copy to maintain consistent filenames
    samtools view -@ "$THREADS" -b "${mapped_bams[0]}" -o "$merged_pre"
  else
    samtools merge -@ "$THREADS" -o "$merged_pre" "${mapped_bams[@]}"
  fi

# Sort merged BAM for Picard and index for convenience
  merged_sorted="${MERGE_DIR}/${sample}.preDedup.merged.sorted.bam"
  samtools sort  -@ "$THREADS" -o "$merged_sorted" "$merged_pre"
  samtools index -@ "$THREADS" "$merged_sorted"

# Define final outputs (deduplicated BAM + duplication metrics)
  final="${MERGE_DIR}/${sample}.final.dedup.bam"
  metrics="${MERGE_DIR}/${sample}.final.dedup.metrics.txt"

# 4) Mark and  remove duplicates with Picard
  if [[ "$REMOVE_DUPLICATES" == "true" ]]; then
    $PICARD_CMD MarkDuplicates \
      I="$merged_sorted" O="$final" M="$metrics" \
      REMOVE_DUPLICATES=true READ_NAME_REGEX='null'
  else
    $PICARD_CMD MarkDuplicates \
      I="$merged_sorted" O="$final" M="$metrics" \
      REMOVE_DUPLICATES=false READ_NAME_REGEX='null'
  fi

# Index final deduplicated BAM for downstream tools
samtools index -@ "$THREADS" "$final"
echo "[DONE] Final BAM for ${sample}: $final"
echo "       Metrics: $metrics"

# Keep only mapped reads for THIS sample
final_mapped="${MERGE_DIR}/${sample}.final.dedup.mapped.bam"
samtools view -@ "$THREADS" -b -F 0x4 "$final" -o "$final_mapped"
samtools index -@ "$THREADS" "$final_mapped"
echo "       Mapped-only BAM: $final_mapped"
echo
done   # ← this closes: for R1 in "${R1_LIST[@]}"; do

echo "Finish time is $(date)"


