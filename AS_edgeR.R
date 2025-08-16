library(Rsubread)
library(edgeR)
library(writexl)

AS_edgeR <- function(Bam_Files,
                     GTF_File,
                     condition_names,
                     pairedEnd,
                     GTF_attribute_type,
                     project_name,
                     folderpath)
  {
  
  edgeR_Starttime <- Sys.time()
  
      ft_counts <- Rsubread::featureCounts(Bam_Files,
                                     annot.ext = GTF_File,
                                     isGTFAnnotationFile = TRUE,
                                     GTF.attrType = GTF_attribute_type,
                                     isPairedEnd = pairedEnd )
      
      edger_DGEList <- DGEList(ft_counts$counts, group = condition_names)
      
      design <- model.matrix(~0+group, data = edger_DGEList$samples)
      
      est_dispersion <- estimateDisp(edger_DGEList,design = design)
      
      plotBCV(est_dispersion)
      
      exact_test <- exactTest(est_dispersion)
      
      Results <- topTags(exact_test,n = 300)
      
      genes = rownames(Results$table)
      
      adj_pvalues = Results$table$FDR
      
      final_results <- data.frame("genes"= genes,"pvalues"=adj_pvalues)
      
      edgeR_endtime <- Sys.time()
      
      total_time_taken <- edgeR_endtime - edgeR_Starttime
      
      rdata_file_path <- file.path(folderpath,project_name,"edgeR", "edgeRtime.RData")
      saveRDS(total_time_taken, file = rdata_file_path)
      
      rdata_file_path <- file.path(folderpath,project_name,"edgeR", "edgeR.RData")
      saveRDS(Results, file = rdata_file_path)
      
      rdata_file_path <- file.path(folderpath,project_name,"edgeR", "edgeR.csv")
      write.csv(final_results, rdata_file_path)
      }
