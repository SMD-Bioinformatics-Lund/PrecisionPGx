#!/usr/bin/env bash
#
# Download PGx-region reads of 4 GeT-RM samples from the 1000 Genomes
# Project 30X high-coverage release (NYGC) and convert them to paired
# FASTQ — the pipeline's natural input. Used to validate the pipeline and
# to reproduce issues #26 and #27 (which only manifest with >=2 parallel
# samples).
#
# Strategy: slicing 45 PGx-gene regions in a single samtools call against
# the remote CRAM is dominated by EBI rate-limiting and CRAM-container
# reseek overhead; per-region slices in parallel are orders of magnitude
# faster. Each region completes in ~10-25 s on its own.
#
# Output (gitignored): test_data/
#   fastqs/<SAMPLE>_R1.fq.gz, <SAMPLE>_R2.fq.gz
#   pgx_genes.bed                # derived from the project's panel BED
#   samplesheet.csv              # ready to feed to --input
#   1000G_2504_high_coverage.sequence.index
#
# Requirements: bash, curl, awk, samtools >= 1.18. Run via
#   pixi exec --spec 'samtools>=1.18' -- bash tests/download_validation_data.sh

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA_DIR="${REPO_ROOT}/test_data"
FASTQ_DIR="${DATA_DIR}/fastqs"
TMP_BAM_DIR="${DATA_DIR}/.tmp_bams"     # per-region BAMs (deleted after merge)
PANEL_BED="${REPO_ROOT}/assets/panpgx.TE-97619858.grch38.refseq_mane.whole_genes_snps.annot.v2.bed"
SLICE_BED="${DATA_DIR}/pgx_genes.bed"
SAMPLESHEET="${DATA_DIR}/samplesheet.csv"

INDEX_URL='http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/1000G_2504_high_coverage.sequence.index'
INDEX_FILE="${DATA_DIR}/1000G_2504_high_coverage.sequence.index"

# GeT-RM PGx reference samples that also live in the 1KG 30X 2504 release.
# Spread: CEU, CEU/GIAB, CHB, YRI — enough parallelism to expose #26 and #27.
SAMPLES=(NA07000 NA12878 NA18526 NA19238)

# How many regions to slice concurrently per sample. Higher values get
# throttled by EBI after a short burst; 4 stays inside their per-IP cap.
REGION_CONCURRENCY=4

mkdir -p "${FASTQ_DIR}" "${TMP_BAM_DIR}"

# CRAM reference lazily fetched from EBI's MD5 service; avoids needing a local
# hg38 fasta for sub-region decoding.
export REF_PATH='https://www.ebi.ac.uk/ena/cram/md5/%s'

if [[ ! -s "${INDEX_FILE}" ]]; then
  echo "Fetching 1KG 30X sequence index..."
  curl -fsSL "${INDEX_URL}" -o "${INDEX_FILE}"
fi

# Derive a gene-level slicing BED from the panel BED — group exon entries by
# gene name prefix, take min start / max end per gene, pad 10 kb on each side.
# ~45 regions, sized to cover star-allele defining variants plus regulatory
# flanks and (for chr22) the CYP2D6/CYP2D7 hybrid spacer.
if [[ ! -s "${SLICE_BED}" || "${PANEL_BED}" -nt "${SLICE_BED}" ]]; then
  awk 'BEGIN { OFS = "\t" }
    {
      split($4, a, "_"); g = a[1]
      if (!(g in chrom)) { chrom[g]=$1; start[g]=$2; end[g]=$3 }
      if ($2 < start[g]) start[g] = $2
      if ($3 > end[g])   end[g]   = $3
    }
    END {
      pad = 10000
      for (g in chrom) {
        s = start[g] - pad; if (s < 0) s = 0
        print chrom[g], s, end[g] + pad, g
      }
    }' "${PANEL_BED}" | sort -k1,1 -k2,2n > "${SLICE_BED}"
  echo "Built gene-level slice BED: $(wc -l < "${SLICE_BED}") regions"
fi

# Build header for the pipeline's samplesheet (schema_input.json contract).
printf 'sample,case_id,type,lane,fastq_1,fastq_2,seq_type,genes\n' > "${SAMPLESHEET}"

for sample in "${SAMPLES[@]}"; do
  r1="${FASTQ_DIR}/${sample}_R1.fq.gz"
  r2="${FASTQ_DIR}/${sample}_R2.fq.gz"

  if [[ -s "${r1}" && -s "${r2}" ]]; then
    echo "[${sample}] FASTQ already present, skipping."
  else
    # Resolve the CRAM URL for this sample.
    cram_url=$(awk -F'\t' -v s="${sample}" '
      /^##/        { next }
      /^#/         { for (i=1;i<=NF;i++) col[$i]=i; next }
      $col["SAMPLE_NAME"]==s && $col["#ENA_FILE_PATH"] ~ /\.final\.cram$/ {
        print $col["#ENA_FILE_PATH"]; exit
      }' "${INDEX_FILE}")

    if [[ -z "${cram_url}" ]]; then
      echo "[${sample}] no .final.cram entry found in index — skipping." >&2
      continue
    fi
    cram_url="${cram_url/ftp:\/\//https://}"

    echo "[${sample}] slicing ${cram_url}"
    echo "[${sample}]   $(wc -l < "${SLICE_BED}") regions, ${REGION_CONCURRENCY} concurrent"

    # Per-region parallel slicing. xargs runs samtools view once per BED line,
    # `--fetch-pairs` brings in mates whose primary read falls outside the
    # region so paired-end reconstruction stays correct.
    sample_tmp="${TMP_BAM_DIR}/${sample}"
    mkdir -p "${sample_tmp}"

    # Resume + retry: slice missing regions, validate with `samtools quickcheck`,
    # delete corrupted BAMs, and re-slice. Loop a few times so transient EBI
    # blips don't kill the whole sample.
    for pass in 1 2 3; do
      # Drop any BAM that's truncated (e.g. samtools view killed mid-write
      # in an earlier session) BEFORE deciding what's missing — the
      # `[[ -s "$out" ]]` test below would accept a nonzero-size but
      # truncated file otherwise. quickcheck exits non-zero on failure;
      # swallow it explicitly so `set -e` does not abort here.
      bad=$(find "${sample_tmp}" -name 'region_*.bam' -print0 2>/dev/null \
              | { xargs -0 -r samtools quickcheck 2>&1 || true; } \
              | awk '/was missing EOF block|is not a sequence file/ {print $1}')
      [[ -n "${bad}" ]] && { echo "${bad}" | xargs -r rm -f; }

      missing=$(
        awk '{ printf "%s\t%s:%s-%s\n", NR, $1, $2+1, $3 }' "${SLICE_BED}" \
          | while IFS=$'\t' read -r i region; do
              out="${sample_tmp}/region_$(printf "%03d" "${i}").bam"
              [[ -s "${out}" ]] || printf "%s\t%s\n" "${i}" "${region}"
            done
      )
      [[ -z "${missing}" ]] && break

      n=$(printf '%s\n' "${missing}" | wc -l)
      echo "[${sample}]   pass ${pass}: slicing ${n} region(s), ${REGION_CONCURRENCY} concurrent"

      printf '%s\n' "${missing}" \
        | xargs -n2 -P "${REGION_CONCURRENCY}" \
            bash -c '
              i="$1"; region="$2"
              out="'"${sample_tmp}"'/region_$(printf "%03d" "$i").bam"
              samtools view --threads 1 -h -b --fetch-pairs \
                "'"${cram_url}"'" "$region" -o "$out" 2>/dev/null
            ' _ || true

      # Drop any BAM that failed quickcheck so the next pass re-slices it.
      # quickcheck exits non-zero on any failure; swallow it explicitly so
      # `set -e` does not abort the script here — that is exactly what we
      # want to inspect.
      bad=$(find "${sample_tmp}" -name 'region_*.bam' -print0 \
              | { xargs -0 -r samtools quickcheck 2>&1 || true; } \
              | awk '/was missing EOF block|is not a sequence file/ {print $1}')
      [[ -n "${bad}" ]] && { echo "${bad}" | xargs -r rm -f; }
    done

    # Merge per-region BAMs (drops duplicates from overlapping regions /
    # mate-fetching automatically when names collide).
    echo "[${sample}]   merging $(ls "${sample_tmp}"/*.bam | wc -l) region BAMs"
    samtools merge --threads 2 -f -c -p -o "${sample_tmp}/merged.bam" "${sample_tmp}"/region_*.bam

    # Coordinate-sorted CRAM -> name-collated -> paired FASTQ.
    echo "[${sample}]   collate + fastq"
    samtools collate -O -u --threads 2 "${sample_tmp}/merged.bam" \
      | samtools fastq --threads 2 -1 "${r1}" -2 "${r2}" -s /dev/null -0 /dev/null -

    rm -rf "${sample_tmp}"
    echo "[${sample}]   R1=$(du -h "${r1}" | cut -f1) R2=$(du -h "${r2}" | cut -f1)"
  fi

  # Each sample is its own case; lane = 1 (single-lane after merge).
  printf '%s,%s,N,1,test_data/fastqs/%s,test_data/fastqs/%s,dna,\n' \
    "${sample}" "${sample}" "$(basename "${r1}")" "$(basename "${r2}")" \
    >> "${SAMPLESHEET}"
done

rmdir "${TMP_BAM_DIR}" 2>/dev/null || true

echo
echo "Done. Sample sheet: ${SAMPLESHEET}"
cat "${SAMPLESHEET}"
