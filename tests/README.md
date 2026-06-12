# Tests

## Validation data

`download_validation_data.sh` downloads PGx-region BAM slices of 4 GeT-RM
PharmCAT reference samples from the 1000 Genomes 30X high-coverage release.

| Sample   | Population | Notes                                  |
| -------- | ---------- | -------------------------------------- |
| NA07000  | CEU        | GeT-RM panel                           |
| NA12878  | CEU        | GIAB + GeT-RM (most-validated sample)  |
| NA18526  | CHB        | GeT-RM panel                           |
| NA19238  | YRI        | GeT-RM panel                           |

Four parallel samples is the minimum needed to reliably trigger the
`.fai`-rewrite race (issue #26) and the BED-vs-VCF channel mispairing under
non-deterministic completion order (issue #27).

### Usage

```bash
pixi exec --spec 'samtools>=1.18' -- bash tests/download_validation_data.sh
```

The script slices each sample's CRAM against
`assets/panpgx.TE-97619858.grch38.refseq_mane.whole_genes_snps.annot.v2.bed`
using a remote CRAM reference (EBI's MD5 service) so no local hg38 fasta is
required for the slice step. Output lands in `test_data/` (gitignored):

```
test_data/
  1000G_2504_high_coverage.sequence.index
  bams/<SAMPLE>.pgx.bam(.bai)
  samplesheet.csv
```

`samplesheet.csv` is pipeline-ready (`schema_input.json` contract).

### Validating against GeT-RM consensus

The published GeT-RM consensus diplotypes (Pratt et al., PharmGKB) are
*not* downloaded by this script. After running the pipeline on these
samples, compare PharmCAT's `.report.tsv` calls per gene against the
published consensus for that sample.
