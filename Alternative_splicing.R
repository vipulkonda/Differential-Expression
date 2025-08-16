# Load supporting scripts for different alternative splicing methods
source("AS_DEXseq.R")
source("AS_ASpli.R")
source("AS_edgeR.R")

# Main wrapper function to run multiple alternative splicing analysis tools
Alternative_splicing <- function(GTF_File,        # Path to GTF annotation file
                                 Bam_Files,       # List of BAM files (aligned reads)
                                 singleEnd,       # Logical: TRUE if single-end reads, FALSE if paired-end
                                 folderpath,      # Output directory path
                                 condition_names, # Vector of experimental condition labels
                                 project_name,    # Project identifier for output naming
                                 GTF_attribute_type){ # GTF attribute type (e.g., "gene_id", "transcript_id")

  # Run DEXSeq-based alternative splicing analysis
  AS_DEXseq(GTF_File, 
            Bam_Files, 
            singleEnd_T_or_F,   # <-- Check: this should be 'singleEnd', unless 'singleEnd_T_or_F' is defined elsewhere
            folderpath,
            condition_names,
            project_name)

  # Define pairedEnd flag based on singleEnd
  pairedEnd <- ifelse(singleEnd, FALSE, TRUE)

  # Run edgeR-based differential splicing analysis
  AS_edgeR(Bam_Files,
           GTF_File,
           condition_names,
           pairedEnd,
           GTF_attribute_type,
           project_name,
           folderpath)

  # Set label for single-end or paired-end for ASpli
  Define_end <- ifelse(singleEnd, "SE", "PE")

  # Run ASpli-based splicing analysis
  AS_aspli(GTF_File,
           Bam_Files,
           folderpath,
           condition_names,
           project_name,
           Define_end)
}
