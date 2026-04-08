#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Create_Samplesheet_GLASS.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Create samplesheet for WES analysis of GLASS samples
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
#  08-04-2026: File creation
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0.1  Load packages
#-------------------------------------------------------------------------------
suppressMessages(library(dplyr))

#-------------------------------------------------------------------------------
# 0.2 Parse command line arguments
#-------------------------------------------------------------------------------
if(exists("snakemake")){
    output <- snakemake@output[[1]]
}else{
    output <- 'samplesheet.csv'
}

#-------------------------------------------------------------------------------
# 1.1 Read data and fetch file paths
#-------------------------------------------------------------------------------
Tumor_files <- list.files('/trinity/home/r115502/GLASS/bam/Tumor/',pattern = '.bam$', full.names = T)
Normal_files <- list.files('/trinity/home/r115502/GLASS/bam/Normal/',pattern = '.bam$', full.names = T)   

#-------------------------------------------------------------------------------
# 1.2 Prepare files/sheet1
#-------------------------------------------------------------------------------
# Fetch patients
Patients <- purrr::map_chr(Tumor_files, ~strsplit(basename(.x),'_')[[1]][1])
Tumors <-  purrr::map_chr(Tumor_files, ~paste(strsplit(basename(.x),'_')[[1]][c(1,2)], collapse='_'))
Normals <-  purrr::map_chr(Normal_files, ~strsplit(basename(.x),'-')[[1]][1])


# Tumor only 
Tumor_only <- data.frame(
    patient = paste0(Tumors,'_TumorOnly'),
    sample = paste0(Tumors,'_TumorOnly_tumor1'),
    lane = 1,
    status = 1,
    bam = Tumor_files)

# Paired
Paired <- data.frame(
    patientID = Patients,
      patient = paste0(Tumors,'_Paired'),
      sample = paste0(Tumors,'_Paired_tumor1'),
      lane = 1,
      status = 1,
    bam = Tumor_files)

# Normals
Normals <- data.frame(
    patientID = Normals,
    sample = 'normal1',
    lane = 1 ,
    status = 0 ,
    bam = Normal_files) 

Normals <- Paired %>% left_join(Normals, by = 'patientID') %>% select(patient,sample.y,lane.y,status.y,bam.y)
colnames(Normals) <- gsub('\\.y','',colnames(Normals))

Paired <- rbind(Paired %>% select(-patientID), Normals)

          

# Create Samplesheet and remove paired without a normal
Samplesheet <- rbind(Tumor_only, Paired) %>%
  group_by(patient)%>%
  filter(all(!is.na(sample))) %>%
  ungroup()




#-------------------------------------------------------------------------------
# 1.2 Merge and write to file
#-------------------------------------------------------------------------------
Samplesheet %>% write.table(file=output, sep = ',', quote = F, row.names = F)
