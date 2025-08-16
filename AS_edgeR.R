# Load required libraries
library(Rsubread)   # For featureCounts (counting reads per feature)
library(edgeR)      # For differential expression analysis
library(writexl)    # For saving results (not directly used here)

# Function: Run edgeR differential analysis on featureCounts results
AS_edgeR <- function(Bam_Files,         # Vector of BAM file paths
                     GTF_File,          # Path to GTF annotation file
                     condition_names,   # Sample condition labels (vector)
                     pairedEnd,         # Logical: TRUE if paired-end sequencing
                     GTF_attribute_type,# Attribute to use from GTF (e.g., "gene_id")
                     project_name,      # Project identifier
                     folderpath)        # Output directory path
{
  # Track runtime
  edgeR_Starttime <- Sys.time()
  
  # Step 1: Count reads per feature using Rsubread::featureCounts
  ft_counts <- Rsubread::featureCounts(Bam_Files,
                                       annot.ext = GTF_File,
                                       isGTFAnnotationFile = TRUE,
                                       GTF.attrType = GTF_attribute_type,
                                       isPairedEnd = pairedEnd)
  
  # Step 2: Wrap counts in edgeR's DGEList object
  edger_DGEList <- DGEList(ft_counts$counts, group = condition_names)
  
  # Step 3: Build design matrix for group comparisons
  design <- model.matrix(~0 + group, data = edger_DGEList$samples)
  
  # Step 4: Estimate dispersion
  est_dispersion <- estimateDisp(edger_DGEList, design = design)
  
  # (Optional) Plot biological coefficient of variation (dispersion plot)
  plotBCV(est_dispersion)
  
  # Step 5: Run exact test for differential expression
  exact_test <- exactTest(est_dispersion)
  
  # Step 6: Extract top differential results
  Results <- topTags(exact_test, n = 300)  # Returns top 300 genes
  
  # Step 7: Collect output table
  genes <- rownames(Results$table)
  adj_pvalues <- Results$table$FDR
  final_results <- data.frame("genes" = genes, "pvalues" = adj_pvalues)
  
  # Track end time
  edgeR_endtime <- Sys.time()
  total_time_taken <- edgeR_endtime - edgeR_Starttime
  
  # Save outputs:
  # - Runtime
  rdata_file_path <- file.path(folderpath, project_name, "edgeR", "edgeRtime.RData")
  saveRDS(total_time_taken, file = rdata_file_path)
  
  # - Full edgeR results object
  rdata_file_path <- file.path(folderpath, project_name, "edgeR", "edgeR.RData")
  saveRDS(Results, file = rdata_file_path)
  
  # - CSV with adjusted p-values
  rdata_file_path <- file.path(folderpath, project_name, "edgeR", "edgeR.csv")
  write.csv(final_results, rdata_file_path)
}
