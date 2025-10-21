# DAP-seq Pipeline (DapSeq-84K)

A reproducible, single-node DAP-seq processing pipeline wrapping common best-practice steps:
quality control and trimming, alignment, duplicate handling, peak calling, and annotation.

## Quick Start

```bash
# Run the pipeline
bash DapSeq-84K.sh
```
---

## Inputs
- Reads: Paired-end FASTQ(.gz); replicate layout expected by the script (e.g., SampleX_R1_rep1.fq.gz, SampleX_R2_rep1.fq.gz, etc.).
- Reference: Genome FASTA (indexed on-the-fly if not present).
- GTF/GFF for downstream R-based annotation.
- Paths are read from variables inside DapSeq-84K.sh. Edit those variables to point to your data.


## Outputs
Cleaned FASTQ, alignment (.sam/.bam), sorted/deduplicated BAMs, peak calls (MACS2), peaks and annotation tables.

### Narrow Peaks Format

The **NarrowPeak** file format is an extended **BED6 + 4** convention, widely produced by **MACS2** for narrow-peak assays such as DAP-seq, ChIP-seq, and transcription-factor binding analyses.
Each record (one line per called peak) contains **10 tab-delimited columns** as described below:

| **Column** | **Name**      | **Description**                                                                    |
| ---------- | ------------- | ---------------------------------------------------------------------------------- |
| 1          | `chrom`       | Chromosome or scaffold identifier                                                  |
| 2          | `chromStart`  | Peak start coordinate (0-based)                                                    |
| 3          | `chromEnd`    | Peak end coordinate (exclusive)                                                    |
| 4          | `name`        | Peak name or unique ID                                                             |
| 5          | `score`       | Integer score (0–1000) for genome browser visualization                            |
| 6          | `strand`      | Strand information (`+`, `−`, or `.` if not applicable)                            |
| 7          | `signalValue` | Enrichment signal or average pileup across the peak                                |
| 8          | `pValue`      | −log₁₀-transformed *p*-value representing peak significance                        |
| 9          | `qValue`      | −log₁₀-transformed *q*-value (FDR-corrected)                                       |
| 10         | `peak`        | Offset (in bp) of the summit position relative to `chromStart` (−1 if unavailable) |

**Example:**

```text
chr1    345000    345250    Peak_001    500    .    25.3    7.2    6.8    120
```

This format ensures compatibility with standard genomic tools and browsers (e.g., UCSC Genome Browser, IGV) and facilitates downstream annotation using utilities such as **ChIPseeker**, **HOMER**, or **BEDTools**.


### IDR Outputs

For reproducible peak-ranking and replicate concordance assessment, the pipeline supports **Irreproducible Discovery Rate (IDR)** analysis.
Detailed specifications and implementation are available in the [NBoley IDR repository](https://github.com/nboley/idr).
 Here’s a clear, publication-quality Markdown section describing the **ChIPseeker** output format — suitable for inclusion in your `README.md` or a methods supplement.

### Annotation Output Format

The **ChIPseeker** R package provides structured peak annotation results in both tabular and graphical forms.
The typical columns include:

| **Name**                | **Description**                                                                          |
| ----------------------- | ---------------------------------------------------------------------------------------- |
| `seqnames`              | Chromosome name (e.g., `chr1`, `chrX`)                                                   |
| `start`                 | Start coordinate of the peak (1-based)                                                   |
| `end`                   | End coordinate of the peak                                                               |
| `width`                 | Peak length (`end - start + 1`)                                                          |
| `strand`                | DNA strand information (`+`, `-`, or `*`)                                                |
| `annotation`            | Functional annotation category (e.g., *Promoter (≤1kb)*, *Exon*, *Intron*, *Intergenic*) |
| `geneChr`               | Chromosome of the associated gene                                                        |
| `geneStart` / `geneEnd` | Gene coordinates                                                                         |
| `geneLength`            | Length of the associated gene                                                            |
| `geneId`                | Ensembl or Entrez gene identifier                                                        |
| `transcriptId`          | Transcript identifier linked to the gene                                                 |
| `distanceToTSS`         | Signed distance from peak center to the transcription start site (TSS)                   |
