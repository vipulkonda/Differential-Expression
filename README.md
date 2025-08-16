Differential-Expression

RNA-seq alternative splicing / differential exon usage pipeline using
DEXSeq, ASpli, and edgeR.
This repo contains four R scripts—one per method plus a driver script—so
you can run tools independently or side-by-side and compare results.

------------------------------------------------------------------------

Repository layout

-   AS_DEXseq.R — DEXSeq workflow for differential exon usage
-   AS_ASpli.R — ASpli workflow for splicing events and bin/junction
    usage
-   AS_edgeR.R — edgeR workflow for feature/junction/bin counts
-   Alternative_splicing.R — optional driver to run multiple methods and
    gather outputs

------------------------------------------------------------------------

Requirements

-   R ≥ 4.2
-   R/Bioconductor packages:
    -   Core: DEXSeq, ASpli, edgeR
    -   Helpers: BiocParallel, GenomicFeatures, rtracklayer,
        SummarizedExperiment, tidyverse, data.table

Install:

    install.packages(c("tidyverse", "data.table"))

    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")

    BiocManager::install(c(
      "DEXSeq", "ASpli", "edgeR",
      "GenomicFeatures", "rtracklayer",
      "SummarizedExperiment", "BiocParallel"
    ))

------------------------------------------------------------------------

Inputs

-   Sample metadata (samples.csv):
    -   Required: sample, condition
    -   Optional: batch, replicate, bam_path (if running from BAMs)
-   Annotation: GTF or TxDb compatible with your genome build
-   Counts:
    -   DEXSeq: exon-level counts or flattened GTF counts
    -   ASpli: counts from bins/junctions or BAMs
    -   edgeR: feature count matrix (features × samples)

Example samples.csv:

    sample,condition,batch,bam_path
    S1,Control,B1,/data/bams/S1.bam
    S2,Control,B1,/data/bams/S2.bam
    S3,Treated,B1,/data/bams/S3.bam
    S4,Treated,B1,/data/bams/S4.bam

------------------------------------------------------------------------

Usage

  Edit the parameter block at the top of each script (paths, conditions,
  etc.), then run:

    # DEXSeq
    Rscript AS_DEXseq.R

    # ASpli
    Rscript AS_ASpli.R

    # edgeR
    Rscript AS_edgeR.R

    # Run multiple methods together
    Rscript Alternative_splicing.R

------------------------------------------------------------------------

Outputs

Each method writes to its own subfolder under results/:

-   Tables: statistics per exon/bin/junction (logFC, p-value, FDR)
-   Plots: library sizes, dispersion trends, MA/volcano, per-gene exon
    usage (DEXSeq)
-   Summaries: significant features and per-gene/event roll-ups

------------------------------------------------------------------------

Workflow

1.  Prepare samples.csv
2.  Use consistent genome build + GTF across all tools
3.  Generate counts appropriate to each method
4.  Run scripts and save to separate output dirs
5.  QC: check low-count filtering, dispersion fit, p-value histograms
6.  Integrate results across tools for confidence
7.  Annotate features → gene symbols or event types

------------------------------------------------------------------------

Reproducibility

-   Save sessionInfo() with results
-   Fix random seeds where relevant
-   Archive exact GTF/TxDb used

------------------------------------------------------------------------

