source("AS_DEXSeq_Aspli_rmats/AS_DEXseq.R")

source("AS_DEXSeq_Aspli_rmats/AS_ASpli.R")

source("AS_DEXSeq_Aspli_rmats/AS_edgeR.R")

Alternative_splicing <- function(GTF_File, 
                                 Bam_Files, 
                                 singleEnd,
                                 folderpath,
                                 condition_names,
                                 project_name,
                                 GTF_attribute_type){


 AS_DEXseq(GTF_File, 
            Bam_Files, 
            singleEnd_T_or_F,
            folderpath,
            condition_names,
            project_name)
  pairedEnd <-  ifelse(singleEnd,FALSE,TRUE)
  AS_edgeR(Bam_Files,
           GTF_File,
           condition_names,
           pairedEnd,
           GTF_attribute_type,
           project_name,
           folderpath
           )
           
  x <- ifelse(singleEnd,"SE","PE")
  AS_aspli(GTF_File,
           Bam_Files,
           folderpath,
           condition_names,
           project_name,
           x)

}



