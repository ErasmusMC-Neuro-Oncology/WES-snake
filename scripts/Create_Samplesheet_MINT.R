#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Create_Samplesheet_MINT.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Create samplesheet for WES analysis of MINT samples
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
#  30-03-2026: File creation
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0.1  Load packages
#-------------------------------------------------------------------------------
suppressMessages(library(dplyr))

#-------------------------------------------------------------------------------
# 0.2 Parse command line arguments
#-------------------------------------------------------------------------------
if(exists("snakemake")){
    data_dir <- snakemake@params[["data_dir"]]
    data_dir2 <- snakemake@params[["data_dir2"]]
    input_sample_overview <- snakemake@params[["sample_overview"]]
    input_sample_overview2 <- snakemake@params[["sample_overview2"]]
    output <- snakemake@output[[1]]
}else{
    data_dir <-  '~/mnt/BIGR_home/MINT/fastq/'
    data_dir2 <-  '~/mnt/BIGR_home/MINT/fastq_batch2/'
    input_sample_overview <- '../data/MINT_db.xlsx'
    input_sample_overview2 <- '../data/1kuvre_fastq_list.csv'

    output <- 'samplesheet.csv'
}

#-------------------------------------------------------------------------------
# 1.1 Read data and fetch file paths
#-------------------------------------------------------------------------------
# Read sample overview
sample_overview <- readxl::read_xlsx(input_sample_overview , skip = 2, sheet = 'DB_Iris')
sample_overview2 <- read.delim(input_sample_overview2, sep = ',')

# Fetch WES .fastq files
fastq_paths <- list.files(data_dir, pattern = '.fastq.gz', recursive = T, full.names = T)
fastq_paths_batch2 <- list.files(data_dir2, pattern = '.fastq.gz', recursive = T, full.names = T)

#-------------------------------------------------------------------------------
# 1.2 Prepare files/sheet1
#-------------------------------------------------------------------------------
# Fetch fastq file IDs
fastq_files <- data.frame(fastq = fastq_paths) %>%
    mutate(GS_ID = purrr::map_chr(fastq, ~strsplit(basename(.x),'_')[[1]][2]),
           read = ifelse(grepl('R1.fastq.gz', fastq),'fastq_1','fastq_2'),
           lane =  purrr::map_chr(fastq, ~strsplit(basename(.x),'_')[[1]][4])) %>%
    tidyr::pivot_wider(names_from = read , values_from = fastq)

# Fetch samplenames
sample_overview <-
    sample_overview %>%
    mutate(patient = gsub('\\(res2\\)','_R2',paste0('MINT',gsub('^M','',`MINT subject nr`)))) %>%
    select(patient, GS_ID, Sentrix_ID) %>%
    group_by(patient) %>%
    mutate(sample = paste0(patient,'_',ifelse(is.na(patient), NA, paste0("tumor", row_number())))) %>%
    ungroup() %>%
    mutate(status = 1)

#-------------------------------------------------------------------------------
# 1.2 Prepare files/sheet\
#-------------------------------------------------------------------------------
fastq_files2 <- data.frame(fastq = fastq_paths_batch2) %>%
    mutate(GS_ID = purrr::map_chr(fastq, ~strsplit(basename(.x),'_')[[1]][2]),
           read = ifelse(grepl('R1.fastq.gz', fastq),'fastq_1','fastq_2'),
           lane =  purrr::map_chr(fastq, ~strsplit(basename(.x),'_')[[1]][4])) %>%
    tidyr::pivot_wider(names_from = read , values_from = fastq) %>%
    mutate(basenameR1 = basename(fastq_1))


sample_overview2 <- sample_overview2 %>%
    mutate(patient = gsub('_M','',purrr::map_chr(RGSM,~strsplit(.x,'-')[[1]][4])),
           sample = paste0(patient,'_',ifelse(is.na(patient), NA, paste0("tumor", row_number()))),
           basenameR1 = basename(Read1File)) %>%
    group_by(patient) %>%
    mutate(sample = paste0(patient,'_',ifelse(is.na(patient), NA, paste0("tumor", row_number())))) %>%
    ungroup() %>%
    filter(!grepl('Normal',patient)) %>%
    mutate(status = 1)


#-------------------------------------------------------------------------------
# 1.2 Merge and write to file
#-------------------------------------------------------------------------------
rbind(
    fastq_files %>% left_join(sample_overview) %>%  select(patient,sample,status,lane,fastq_1,fastq_2),
    fastq_files2 %>% left_join(sample_overview2)  %>%  select(patient,sample,status,lane,fastq_1,fastq_2)
) %>% filter(!is.na(patient)) %>%  write.table(file=output, sep = ',', quote = F, row.names = F)
