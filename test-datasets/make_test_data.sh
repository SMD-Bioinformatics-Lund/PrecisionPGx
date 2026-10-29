#!/usr/bin/env bash
#
# Build a small GRCh38 test dataset for the PrecisionPGx test profiles.
#
# Reads from the GIAB HG002 60x Illumina GRCh38 BAM are extracted over the
# PGx target regions in assets/ (~66 kb, padded by 500 bp), giving a BAM and a
# FASTQ pair of a few MB each. The full GRCh38 reference and the PharmCAT
# reference are also needed.
#
# Usage:
#   bash test-datasets/make_test_data.sh [-o OUTDIR] [-f GRCh38_FASTA]
#
#   -o  Output directory (default: test-datasets/data in this repo)
#   -f  An existing GRCh38 FASTA (chr-prefixed, no alts) to use instead of
#       downloading one. A .fai next to it is reused if present.
#
# Requires: samtools >= 1.10, curl, tar, gzip
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(dirname "${SCRIPT_DIR}")"
OUTDIR="${SCRIPT_DIR}/data"
FASTA=""

while getopts "o:f:h" opt; do
    case ${opt} in
        o) OUTDIR="${OPTARG}" ;;
        f) FASTA="${OPTARG}" ;;
        h) sed -n '2,18p' "$0"; exit 0 ;;
        *) sed -n '2,18p' "$0"; exit 1 ;;
    esac
done

TARGET_BED="${REPO_DIR}/assets/panpgx.TE-97619858.grch38.refseq_mane.whole_genes_snps.annot.v2.bed"
GIAB_BAM="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.60x.1.bam"
GRCH38_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
PHARMCAT_REF_URL="https://zenodo.org/records/7288118/files/GRCh38_reference_fasta.tar"
SAMPLE="HG002"
GENES="CYP2C19,CYP2C9,CYP2D6,DPYD,TPMT,SLCO1B1"

command -v samtools >/dev/null || { echo "ERROR: samtools not found on PATH (try: module load samtools)" >&2; exit 1; }

mkdir -p "${OUTDIR}"
OUTDIR="$(cd "${OUTDIR}" && pwd)"
cd "${OUTDIR}"

#
# 1. GRCh38 reference genome
#
mkdir -p reference
if [[ -n "${FASTA}" ]]; then
    echo "[1/4] Linking existing reference ${FASTA}"
    ln -sf "$(cd "$(dirname "${FASTA}")" && pwd)/$(basename "${FASTA}")" reference/genome.fa
    if [[ -s "${FASTA}.fai" ]]; then
        ln -sf "$(cd "$(dirname "${FASTA}")" && pwd)/$(basename "${FASTA}").fai" reference/genome.fa.fai
    fi
elif [[ ! -s reference/genome.fa ]]; then
    echo "[1/4] Downloading GRCh38 no-alt analysis set (~0.9 GB)"
    curl -fL --retry 3 -o reference/genome.fa.gz "${GRCH38_URL}"
    gunzip reference/genome.fa.gz
else
    echo "[1/4] reference/genome.fa exists, skipping"
fi
[[ -s reference/genome.fa.fai ]] || samtools faidx reference/genome.fa

#
# 2. Reads over the PGx targets
#
# Picard and GATK require the
# BAM header to match the reference exactly, so the header is rebuilt from the
# reference .fai, dropping reads whose mate maps to a contig not in the reference.
#
if [[ ! -s "${SAMPLE}.pgx.bam.bai" ]]; then
    echo "[2/4] Extracting ${SAMPLE} reads over PGx targets"
    awk 'BEGIN{OFS="\t"} {s=$2-500; if (s<0) s=0; print $1, s, $3+500}' "${TARGET_BED}" > regions.bed
    # htslib downloads the remote .bai into the working directory
    samtools view -b -M -L regions.bed -o "${SAMPLE}.giab.bam" "${GIAB_BAM}"
    {
        printf '@HD\tVN:1.6\tSO:unsorted\n'
        awk 'BEGIN{OFS="\t"} {print "@SQ", "SN:"$1, "LN:"$2}' reference/genome.fa.fai
        samtools view -H "${SAMPLE}.giab.bam" | grep -v -e '^@HD' -e '^@SQ'
        samtools view "${SAMPLE}.giab.bam" \
            | awk -F'\t' 'NR==FNR {keep[$1]=1; next} ($3 in keep) && ($7=="=" || $7=="*" || ($7 in keep))' \
                reference/genome.fa.fai -
    } | samtools sort -u - \
      | samtools addreplacerg -r "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tPL:ILLUMINA" -o "${SAMPLE}.pgx.bam" -
    samtools index -o "${SAMPLE}.pgx.bam.bai" "${SAMPLE}.pgx.bam"
    rm -f "${SAMPLE}.giab.bam" HG002.GRCh38.60x.1.bam.bai
else
    echo "[2/4] ${SAMPLE}.pgx.bam exists, skipping"
fi

#
# 3. FASTQ pairs from the same reads (pairs with a mate outside the targets are dropped)
#
# The read pairs are split alternately into two lanes, so the test profile also
# covers merging of multi-lane samples.
#
if [[ ! -s "${SAMPLE}_L002_R2.fastq.gz" ]]; then
    echo "[3/4] Converting to paired FASTQ, split into two lanes"
    samtools collate -u -O "${SAMPLE}.pgx.bam" "${SAMPLE}.collate" \
        | samtools fastq -n -1 "${SAMPLE}_R1.fastq" -2 "${SAMPLE}_R2.fastq" -0 /dev/null -s /dev/null -
    for READ in R1 R2; do
        awk -v out="${SAMPLE}" -v read="${READ}" \
            '{ lane = int((NR - 1) / 4) % 2 + 1; print > (out "_L00" lane "_" read ".fastq") }' \
            "${SAMPLE}_${READ}.fastq"
        rm "${SAMPLE}_${READ}.fastq"
    done
    gzip -f "${SAMPLE}"_L00[12]_R[12].fastq
else
    echo "[3/4] FASTQ files exist, skipping"
fi

#
# 4. PharmCAT reference genome
#
if [[ -z "$(find pharmcat_reference -name '*.bgz' 2>/dev/null)" ]]; then
    echo "[4/4] Downloading PharmCAT reference (~0.9 GB)"
    mkdir -p pharmcat_reference
    curl -fL --retry 3 -o pharmcat_reference.tar "${PHARMCAT_REF_URL}"
    tar -xf pharmcat_reference.tar -C pharmcat_reference
    rm pharmcat_reference.tar
else
    echo "[4/4] PharmCAT reference exists, skipping"
fi
PC_FASTA="$(find "${OUTDIR}/pharmcat_reference" -name '*.bgz' | head -1)"
for f in "${PC_FASTA}" "${PC_FASTA}.fai" "${PC_FASTA}.gzi"; do
    [[ -s "${f}" ]] || { echo "ERROR: missing ${f} in pharmcat_reference/" >&2; exit 1; }
done
ln -sf "${PC_FASTA}"     pharmcat_reference.fna.bgz
ln -sf "${PC_FASTA}.fai" pharmcat_reference.fna.bgz.fai
ln -sf "${PC_FASTA}.gzi" pharmcat_reference.fna.bgz.gzi

#
# Samplesheets (absolute paths)
#
cat > samplesheet_fastq.csv <<EOF
sample,case_id,type,lane,fastq_1,fastq_2,seq_type,genes
${SAMPLE},${SAMPLE},N,L001,${OUTDIR}/${SAMPLE}_L001_R1.fastq.gz,${OUTDIR}/${SAMPLE}_L001_R2.fastq.gz,dna,"${GENES}"
${SAMPLE},${SAMPLE},N,L002,${OUTDIR}/${SAMPLE}_L002_R1.fastq.gz,${OUTDIR}/${SAMPLE}_L002_R2.fastq.gz,dna,"${GENES}"
EOF

cat > samplesheet_bam.csv <<EOF
sample,case_id,type,lane,bam,bai,seq_type,genes
${SAMPLE},${SAMPLE},N,L001,${OUTDIR}/${SAMPLE}.pgx.bam,${OUTDIR}/${SAMPLE}.pgx.bam.bai,dna,"${GENES}"
EOF

echo
echo "Test data ready in ${OUTDIR}"
if [[ "${OUTDIR}" != "${SCRIPT_DIR}/data" ]]; then
    echo "Pass --pipelines_testdata_base_path ${OUTDIR} when running the test profiles."
fi
