library(GenomicAlignments)
library(GenomicFeatures)
library(BiocParallel)
library(DEXSeq)
library(writexl)
#DEXseq
AS_DEXseq <- function(GTF_File, 
                      Bam_Files, 
                      singleEnd_T_or_F,
                      folderpath,
                      condition_names,
                      project_name)
{
  dir.create(file.path(folderpath, project_name))
  Dexseq_path <- dir.create(file.path(folderpath,project_name,"DEXSeq"))
  Dexseq_path <- dir.create(file.path(folderpath,project_name,"ASpli"))
  Dexseq_path <- dir.create(file.path(folderpath,project_name,"edgeR"))
  DEXSeqStarttime <- Sys.time()
  txdb <- makeTxDbFromGFF(GTF_File)
  ebg <-exonicParts( txdb, linked.to.single.gene.only = TRUE )
  head(ebg)
  se <- summarizeOverlaps(features=ebg, reads=Bam_Files, mode="Union", 
                          singleEnd=as.logical(singleEnd), 
                          ignore.strand=TRUE, fragments=FALSE )
  head(se)
  colData(se)$condition =condition_names
  dxd <- DEXSeqDataSetFromSE( se, design= ~ sample + exon + condition:exon )
  head(dxd)
  split( seq_len(ncol(dxd)), colData(dxd)$exon )
  #head( featureCounts(dxd), 5 )
  #head( rowRanges(dxd), 3 )
  sampleAnnotation( dxd )
  dxd = estimateSizeFactors( dxd )
  head(dxd)
  dxd = estimateDispersions( dxd )
  head(dxd)
  dxd = testForDEU( dxd )
  head(dxd)
  dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition")
  head(dxd)
  dxr1 = DEXSeqResults( dxd )
  dxr1
  DEXSeqendtime <- Sys.time()
  totaltimeDEXseq <- DEXSeqendtime - DEXSeqStarttime
  totaltimeDEXseq
  pvalues <- data.frame(dxr1$groupID,dxr1$featureID,dxr1$padj)

  
  #plotCounts(top_genes, gene =$padj , intgroup = "condition")

###
  rdata_file_path <- file.path(folderpath,project_name,"DEXSeq", "Dexseq.RData")
  saveRDS(dxr1, file = rdata_file_path)
  
  rdata_file_path <- file.path(folderpath,project_name,"DEXSeq", "Dexseqtime.RData")
  saveRDS(totaltimeDEXseq, file = rdata_file_path)
  
  rdata_file_path <- file.path(folderpath,project_name,"DEXSeq", "Dexseq_Pvalues.csv")
  write.csv(pvalues, rdata_file_path)
  
}
