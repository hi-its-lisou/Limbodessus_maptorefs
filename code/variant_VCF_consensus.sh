#!/usr/bin/env bash
#SBATCH -J vcf_and_exon_consensus
#SBATCH -p batch
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=04:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=elisabethwilliamsonn@gmail.com

set -euo pipefail

########################################
# CONFIG  (USE SAME BUILD AS MAPPING)
########################################
IN_DIR="/hpcfs/users/a1707473/map_to_refs/merged"
OUT_DIR="/hpcfs/users/a1707473/map_to_refs/filtered"
PATTERN="*.final.dedup.mapped.bam"
MAPQ_THRESH=30
THREADS="${SLURM_CPUS_PER_TASK:-4}"

# EXACT reference + BED used for BWA mapping
REF="/hpcfs/users/a1707473/allgenes_noUTR_concat.fasta"
TARGETS="/hpcfs/users/a1707473/allgenes_noUTR_concat.bed"

# Per-exon consensus options
MIN_COV=2          # mask positions with depth < MIN_COV INSIDE exons
IUPAC="true"       # heterozygotes as IUPAC (bcftools -I). Set "false" for major allele.

# Output layout
POST_DIR="${OUT_DIR%/}/postproc"
VAR_DIR="${POST_DIR}/variants"
COV_DIR="${POST_DIR}/coverage"
MASK_DIR="${POST_DIR}/masks"
EXON_DIR="${POST_DIR}/consensus_per_exon_bcftools"

mkdir -p "$OUT_DIR" "$POST_DIR" "$VAR_DIR" "$COV_DIR" "$MASK_DIR" "$EXON_DIR"

########################################
# MODULES / REQUIREMENTS
########################################
module load SAMtools/1.17-GCC-11.2.0
module load freebayes/1.3.6-foss-2021b-R-4.2.3
module load BCFtools/1.17-GCC-11.2.0 || module load bcftools || true

require_cmd() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: missing '$1'" >&2; exit 1; }; }
require_file() { [[ -s "$1" ]] || { echo "ERROR: missing file: $1" >&2; exit 1; }; }

require_cmd samtools
require_cmd freebayes
require_cmd bcftools
require_cmd bgzip
require_cmd tabix
require_cmd rev   # used for reverse-complement

require_file "$REF"
require_file "$TARGETS"

echo
echo "Start: $(date)"
echo "REF:     $REF"
echo "TARGETS: $TARGETS"
echo "MAPQ ≥ $MAPQ_THRESH | MIN_COV $MIN_COV | IUPAC=$IUPAC"
echo

########################################
# PREP: index FASTA & make 1-based regions for bcftools -R/-r
########################################
[[ -s "${REF}.fai" ]] || samtools faidx "$REF"

TARGETS_1B="${POST_DIR}/$(basename "$TARGETS").1based.tsv"
awk -v OFS='\t' '
  /^#/ || /^track/ || /^browser/ {next}
  NF>=3 {print $1, $2+1, $3}
' "$TARGETS" > "$TARGETS_1B"

# Does BED have strand?
bed_has_strand=0
awk 'BEGIN{s=1} /^#/||/^track/||/^browser/{next} NF>=6{exit 0} END{exit 1}' "$TARGETS" && bed_has_strand=1

# Tag coverage cache by BED hash to avoid stale reuse when TARGETS changes
COV_TAG=$(md5sum "$TARGETS" | awk '{print substr($1,1,8)}')

# small helper for reverse-complement (single line)
rc() { rev | tr 'ACGTNacgtn' 'TGCANtgcan'; }

########################################
# MAIN
########################################
readarray -t BAM_LIST < <(find "$IN_DIR" -maxdepth 1 -type f -name "$PATTERN" | sort)
(( ${#BAM_LIST[@]} > 0 )) || { echo "ERROR: No BAMs in $IN_DIR matching $PATTERN" >&2; exit 1; }

echo "Found ${#BAM_LIST[@]} BAM(s):"
printf '  %s\n' "${BAM_LIST[@]}"
echo

for bam in "${BAM_LIST[@]}"; do
  base="$(basename "$bam")"
  sample="${base%.bam}"                              # e.g., X.final.dedup.mapped
  filt="${OUT_DIR}/${sample}.MQ${MAPQ_THRESH}.bam"   # filtered BAM
  stem="$(basename "${filt%.bam}")"                  # used to name outputs

  echo "========== Processing: $base =========="

  # 1) MAPQ filter BAM (kept for provenance; fast)
  samtools view -@ "$THREADS" -b -q "$MAPQ_THRESH" -o "$filt" "$bam"
  samtools index -@ "$THREADS" "$filt"

  # 2) Call variants (FREEBAYES)  -> keep _freebayes.vcf
  vcf_free="${VAR_DIR}/${stem}_freebayes.vcf"
  echo "[freebayes] -> $vcf_free"
  freebayes \
    --min-base-quality 20 \
    --min-alternate-fraction 0.2 \
    --min-alternate-count 2 \
    -f "$REF" "$filt" > "$vcf_free"

  # bgzip + index TEMPORARILY so bcftools can do region operations
  vcf_free_gz="${vcf_free}.gz"
  bgzip -c "$vcf_free" > "$vcf_free_gz"
  tabix -p vcf "$vcf_free_gz"

  # 3) On-target VCF via bcftools view -R (keep as _vcf_intersectBed.vcf)
  vcf_bed="${VAR_DIR}/${stem}_vcf_intersectBed.vcf"
  echo "[bcftools view -R] -> $vcf_bed"
  bcftools view -R "$TARGETS_1B" -Ov -o "$vcf_bed" "$vcf_free_gz"

  # 4) Build low-coverage mask inside targets (no outside-targets masking)
  cov="${COV_DIR}/${stem}_bases.${COV_TAG}.txt"
  if [[ ! -s "$cov" ]]; then
    echo "Generating coverage table (within TARGETS): $cov"
    samtools depth -@ "$THREADS" -aa -J -b "$TARGETS" "$filt" > "$cov"
  fi

  lowcov="${MASK_DIR}/${stem}.lowcov_lt${MIN_COV}.${COV_TAG}.bed"
  awk -v M="$MIN_COV" '
    $3 < M {
      if (cur=="" || $1!=chr || $2!=end+1) {
        if (cur!="") print chr, start-1, end;
        chr=$1; start=$2; end=$2; cur=1
      } else { end=$2 }
    }
    $3 >= M { if (cur!="") { print chr, start-1, end; cur="" } }
    END { if (cur!="") print chr, start-1, end }
  ' OFS='\t' "$cov" > "$lowcov"

  # 5) Per-exon consensus directly from REF+VCF (no bedtools)
  exonfa="${EXON_DIR}/${stem}.per_exon.fa"
  : > "$exonfa"

  # Header prefix "Genus_species" from sample (strip ".final.dedup...")
  sample_core="${sample%.final.dedup*}"
  prefix="$(printf '%s\n' "$sample_core" | awk -F'_' '{if(NF>=2) print $1"_"$2; else print $1}')"

  awk -v HAS_STR="$bed_has_strand" '
    /^#/ || /^track/ || /^browser/ { next }
    NF>=3 {
      chr=$1; start=$2; end=$3;
      name=(NF>=4 && $4!="") ? $4 : chr":"(start+1)"-"end;
      strand=(HAS_STR && NF>=6) ? $6 : "+";
      printf "%s\t%d\t%d\t%s\t%s\n", chr, start, end, name, strand
    }
  ' "$TARGETS" \
  | while IFS=$'\t' read -r chr start end name strand; do
      beg=$((start+1))
      reg="${chr}:${beg}-${end}"
      iup_opt=(); [[ "$IUPAC" == "true" ]] && iup_opt=(-I)

# Generate consensus for JUST this exon; mask low coverage; FASTA region comes from samtools faidx
seq=$(
  samtools faidx "$REF" "$reg" \
  | bcftools consensus -m "$lowcov" "${iup_opt[@]}" "$vcf_free_gz" \
  | sed '1d' | tr -d '\n'
)

      # Fallback to Ns if nothing emitted (rare)
      want_len=$((end - start))
      if [[ -z "$seq" ]]; then
        seq=$(printf '%.0sN' $(seq 1 "$want_len"))
      fi

      # Respect negative strand if present
      if [[ "$strand" == "-" ]]; then
        seq=$(printf "%s" "$seq" | rc)
      fi

      clean_name=$(printf "%s" "$name" | tr ' \t' '_')
      printf ">%s|%s\n%s\n" "$prefix" "$clean_name" "$seq" >> "$exonfa"
    done

  echo "[DONE] Wrote:"
  echo "  - $vcf_free"
  echo "  - $vcf_bed"
  echo "  - $exonfa"
  echo

  # Remove temporary compressed index to leave only the two VCFs requested
  rm -f "$vcf_free_gz" "${vcf_free_gz}.tbi"

done

echo "Finish: $(date)"

 
