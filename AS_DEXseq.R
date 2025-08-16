# Load required libraries
library(GenomicAlignments)  # For summarizeOverlaps
library(GenomicFeatures)    # For TxDb and exonicParts
library(BiocParallel)       # For parallelization (not directly used here)
library(DEXSeq)             # Main package for differential exon usage
library(writexl)            # For exporting results if needed

# Function: Run DEXSeq alternative splicing analysis
AS_DEXseq <- function(GTF_File,          # Path to GTF annotation file
                      Bam_Files,         # List of BAM files (aligned reads)
                      singleEnd_T_or_F,  # TRUE if single-end, FALSE if paired-end
                      folderpath,        # Output directory path
                      condition_names,   # Vector of sample conditions
                      project_name)      # Project identifier
{
  # Create project and subdirectories for storing results
  dir.create(file.path(folderpath, project_name))
  dir.create(file.path(folderpath, project_name, "DEXSeq"))
  dir.create(file.path(folderpath, project_name, "ASpli"))   
  dir.create(file.path(folderpath, project_name, "edgeR"))   
  
  # Track start time
  DEXSeqStarttime <- Sys.time()
  
  # Build transcript database from GTF file
  txdb <- makeTxDbFromGFF(GTF_File)
  
  # Extract exonic parts (flattened exons grouped per gene)
  ebg <- exonicParts(txdb, linked.to.single.gene.only = TRUE)
  head(ebg) # Preview
  
  # Count reads overlapping exonic parts
  se <- summarizeOverlaps(features = ebg,
                          reads = Bam_Files,
                          mode = "Union",
                          singleEnd = as.logical(singleEnd_T_or_F), 
                          ignore.strand = TRUE,
                          fragments = FALSE)
  head(se)
  
  # Add experimental conditions
  colData(se)$condition <- condition_names
  
  # Build DEXSeq dataset
  dxd <- DEXSeqDataSetFromSE(se, design = ~ sample + exon + condition:exon)
  head(dxd)
  
  # Split columns by exon (internal check)
  split(seq_len(ncol(dxd)), colData(dxd)$exon)
  
  # Inspect feature counts / annotations (optional)
  # head(featureCounts(dxd), 5)
  # head(rowRanges(dxd), 3)
  sampleAnnotation(dxd)
  
  # Normalize data
  dxd <- estimateSizeFactors(dxd)
  head(dxd)
  
  # Estimate dispersion
  dxd <- estimateDispersions(dxd)
  head(dxd)
  
  # Test for differential exon usage (DEU)
  dxd <- testForDEU(dxd)
  head(dxd)
  
  # Estimate exon fold changes relative to condition
  dxd <- estimateExonFoldChanges(dxd, fitExpToVar = "condition")
  head(dxd)
  
  # Extract results
  dxr1 <- DEXSeqResults(dxd)
  dxr1
  
  # Track end time
  DEXSeqendtime <- Sys.time()
  totaltimeDEXseq <- DEXSeqendtime - DEXSeqStarttime
  totaltimeDEXseq
  
  # Extract p-values (adjusted)
  pvalues <- data.frame(dxr1$groupID, dxr1$featureID, dxr1$padj)
  
  # Save outputs
  rdata_file_path <- file.path(folderpath, project_name, "DEXSeq", "Dexseq.RData")
  saveRDS(dxr1, file = rdata_file_path)
  
  rdata_file_path <- file.path(folderpath, project_name, "DEXSeq", "Dexseqtime.RData")
  saveRDS(totaltimeDEXseq, file = rdata_file_path)
  
  rdata_file_path <- file.path(folderpath, project_name, "DEXSeq", "Dexseq_Pvalues.csv")
  write.csv(pvalues, rdata_file_path)
}
