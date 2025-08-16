# Load required libraries
library(ASpli)          # Alternative splicing analysis toolkit
library(GenomicFeatures) # For working with genomic annotations (TxDb objects)
library(writexl)        # For writing Excel files if needed

# Function: Run ASpli alternative splicing analysis
AS_aspli <- function(GTF_File,        # Path to the GTF annotation file
                     Bam_Files,       # Vector of BAM files (aligned RNA-seq reads)
                     folderpath,      # Output directory path
                     condition_names, # Vector of sample conditions (e.g., WT, KO)
                     project_name,    # Project identifier (used in output naming)
                     Define_end)               # Library type: "SE" (single-end) or "PE" (paired-end)
{
  # Track runtime
  ASpliStarttime <- Sys.time()
  
  # Create output directory for ASpli results
  Dexseq_path <- dir.create(paste0(folderpath, project_name, "/ASpli"))
  
  
  # Build transcript database object from GTF
  txdb <- makeTxDbFromGFF(GTF_File)
  
  # Generate genomic bins (exons, introns, junctions, etc.)
  features <- binGenome(txdb)
  
  # Prepare metadata (targets file) linking samples with BAM files & conditions
  targets <- data.frame(row.names = paste0('Sample', c(1:length(Bam_Files))),
                        bam = Bam_Files, 
                        f1 = condition_names,
                        stringsAsFactors = FALSE)
  
  # Count reads in genomic bins
  gbcounts <- gbCounts(features = features,
                       targets = targets,
                       libType = Define_end,         # SE or PE
                       minReadLength = 1,
                       maxISize = 300)
  
  # Count junction reads
  asd <- jCounts(counts = gbcounts,
                 features = features,
                 libType = x,
                 threshold = 10,
                 minReadLength = 1)
  
  # Differential usage analysis for bins
  gb <- gbDUreport(gbcounts, contrast = c(-1, 1))
  
  # Differential usage analysis for junctions
  jdur <- jDUreport(asd, contrast = c(-1, 1))
  
  # Create splicing report combining bins + junctions
  sr <- splicingReport(gb, jdur, counts = gbcounts)
  
  # Integrate bin and junction signals
  is <- integrateSignals(sr, asd)
  
  # Extract p-values / FDR for splicing signals
  pvalues <- data.frame(is@signals$bin, is@signals$b.fdr)
  
  # Track runtime again
  ASpliendtime <- Sys.time()
  totaltimeASpli <- ASpliendtime - ASpliStarttime
  totaltimeASpli
  
  # Save results (RData objects and CSV)
  rdata_file_path <- file.path(folderpath, project_name, "ASpli", "ASpli.RData")
  saveRDS(is, file = rdata_file_path)
  
  rdata_file_path <- file.path(folderpath, project_name, "ASpli", "ASplitime.RData")
  saveRDS(totaltimeASpli, file = rdata_file_path)
  
  rdata_file_path <- file.path(folderpath, project_name, "ASpli", "ASpli.csv")
  write.csv(pvalues, rdata_file_path)
}
