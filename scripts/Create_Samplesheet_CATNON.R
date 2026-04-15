#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Create_Samplesheet_CATNON.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Create samplesheet for panel analysis of CATNON samples
#
# Author: Jurriaan Janssen (j.janssen.1@erasmusmc.nl)
#
# condaenv: 
# Usage: 
#
# TODO:
# 1) 
#
# History:
#  14-04-2026: File creation
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0.1  Load packages
#-------------------------------------------------------------------------------
suppressMessages(library(dplyr))

#-------------------------------------------------------------------------------
# 0.2 Parse command line arguments
#-------------------------------------------------------------------------------
if(exists("snakemake")){
    data_dir <- snakemake@params[["data_dir"]]
    output <- snakemake@output[[1]]
}else{
    data_dir <-  '~/mnt/BIGR_home/CATNON/bam_IDHmt/'
    output <- 'samplesheet.csv'
}

#-------------------------------------------------------------------------------
# 1.1 Read data and fetch file paths
#-------------------------------------------------------------------------------
# Fetch WES .bam files
bam_paths <- list.files(data_dir, pattern = '.bam$', recursive = T, full.names = T)

#-------------------------------------------------------------------------------
# 1.2 Prepare files/sheet1
#-------------------------------------------------------------------------------
bam_files <- data.frame(bam = bam_paths) %>%
    mutate(patient = purrr::map_chr(bam, ~strsplit(basename(.x),'\\.')[[1]][1]),
           sample = paste0(patient,'_tumor1'),
           lane =  1,
           status= 1)
#-------------------------------------------------------------------------------
# 1.2 Merge and write to file
#-------------------------------------------------------------------------------
bam_files %>%  select(patient,sample,status,lane,bam) %>% write.table(file=output, sep = ',', quote = F, row.names = F)
