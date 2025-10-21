# DAP-seq Pipeline (DapSeq-84K)

A reproducible, single-node DAP-seq processing pipeline wrapping common best-practice steps:
quality control and trimming, alignment, duplicate handling, peak calling, and annotation.

## Quick Start

```bash
# Run the pipeline
bash DapSeq-84K.sh
```
## Input / Output
### Inputs
- Reads: Paired-end FASTQ(.gz); replicate layout expected by the script (e.g., SampleX_R1_rep1.fq.gz, SampleX_R2_rep1.fq.gz, etc.).
- Reference: Genome FASTA (indexed on-the-fly if not present).
- GTF/GFF for downstream R-based annotation.
- Paths are read from variables inside DapSeq-84K.sh. Edit those variables to point to your data.

### Outputs
Cleaned FASTQ, alignment (.sam/.bam), sorted/deduplicated BAMs, peak calls (MACS2), peaks and annotation tables.
- Narrow Peaks Format
The NarrowPeak file format is an extended BED6+4 format commonly produced by MACS2 for narrow peak detection assays such as DAP-seq, ChIP-seq, and transcription factor binding analysis. Each line represents a called peak and contains 10 tab-delimited columns:
| Column | Name          | Description                                                                          |
| ------ | ------------- | ------------------------------------------------------------------------------------ |
| 1      | `chrom`       | Chromosome or scaffold identifier                                                    |
| 2      | `chromStart`  | Peak start coordinate (0-based)                                                      |
| 3      | `chromEnd`    | Peak end coordinate (exclusive)                                                      |
| 4      | `name`        | Peak name or ID                                                                      |
| 5      | `score`       | Integer score (0–1000) for browser display                                           |
| 6      | `strand`      | Strand information (`+`, `-`, or `.` if not applicable)                              |
| 7      | `signalValue` | Enrichment signal or average pileup across the peak                                  |
| 8      | `pValue`      | –log10 transformed p-value for the peak significance                                 |
| 9      | `qValue`      | –log10 transformed q-value (FDR corrected)                                           |
| 10     | `peak`        | Offset (in bp) of the summit position relative to `chromStart` (–1 if not available) |
> chr1    345000    345250    Peak_001    500    .    25.3    7.2    6.8    120

- *IDR Outputs:* reference https://github.com/nboley/idr
