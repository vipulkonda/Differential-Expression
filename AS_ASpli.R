library(ASpli)
library(GenomicFeatures)
library(writexl)

#ASpli
AS_aspli <- function(GTF_File,
                     Bam_Files,
                     folderpath,
                     condition_names,
                     project_name,
                     x)
{
  ASpliStarttime <- Sys.time()
  
  Dexseq_path <- dir.create(paste0(folderpath, project_name, "/ASpli"))
  
  print(singleEnd)

  txdb <- makeTxDbFromGFF(GTF_File)
  
  features<-binGenome(txdb)
  
  targets<-data.frame(row.names=paste0('Sample',c(1:length(Bam_Files))),
                      bam=Bam_Files, 
                      f1=condition_names,
                      stringsAsFactors=FALSE)
  
  gbcounts<-gbCounts(features=features,
                     targets=targets,
                     libType = x,
                     minReadLength=1,
                     maxISize=300)
  
  asd<-jCounts(counts=gbcounts,
               features=features,
               libType = x,
               threshold = 10,
               minReadLength=1)
  
  gb<-gbDUreport(gbcounts,contrast=c(-1,1))
  
  jdur<-jDUreport(asd,contrast=c(-1,1))
  
  sr<-splicingReport(gb,jdur,counts=gbcounts)
  
  is<-integrateSignals(sr,asd)
  
  pvalues <- data.frame(is@signals$bin,is@signals$b.fdr)
  
  ASpliendtime <- Sys.time()
  totaltimeASpli <- ASpliendtime - ASpliStarttime
  totaltimeASpli
  
  rdata_file_path <- file.path(folderpath,project_name,"ASpli", "ASpli.RData")
  saveRDS(is, file = rdata_file_path)
  
  rdata_file_path <- file.path(folderpath,project_name,"ASpli", "ASplitime.RData")
  saveRDS(totaltimeASpli, file = rdata_file_path)
  
  rdata_file_path <- file.path(folderpath,project_name,"ASpli", "ASpli.csv")
  write.csv(pvalues, rdata_file_path)
}
