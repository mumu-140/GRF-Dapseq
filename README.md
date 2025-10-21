# DAP-seq Pipeline (DapSeq-84K)

A reproducible, single-node DAP-seq processing pipeline wrapping common best-practice steps:
quality control and trimming, alignment, duplicate handling, peak calling, and annotation.

## Quick Start

```bash
# Run the pipeline
bash DapSeq-84K.sh

## Input / Output
### Inputs
- Reads: Paired-end FASTQ(.gz); replicate layout expected by the script (e.g., SampleX_R1_rep1.fq.gz, SampleX_R2_rep1.fq.gz, etc.).
- Reference: Genome FASTA (indexed on-the-fly if not present).
- GTF/GFF for downstream R-based annotation.
- Paths are read from variables inside DapSeq-84K.sh. Edit those variables to point to your data.

### Outputs
Cleaned FASTQ, alignment (.sam/.bam), sorted/deduplicated BAMs, peak calls (MACS2), peaks and annotation tables.
