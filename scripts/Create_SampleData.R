#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Create_SampleData.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Combine WES pipeline results into sample data
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
#  05-03-2026: File creation
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0.1  Load packages
#-------------------------------------------------------------------------------
suppressMessages(library(dplyr))
suppressMessages(library(stringr))

#-------------------------------------------------------------------------------
# 0.2 Parse command line arguments
#-------------------------------------------------------------------------------
if(exists("snakemake")){
    input_coverage <- snakemake@input[["on_target_coverage"]]
    input_CNA_stats <- snakemake@input[["CNA_stats"]]
    input_ACE_results <- snakemake@input[["ACE_results"]]
    input_CNH_results <- snakemake@input[["CNH_results"]]
    input_TMB <- snakemake@input[["TMB"]]
    input_signatures <- snakemake@input[["signatures"]]
    output_SampleData <- snakemake@output[["SampleData"]]
}else{
    input_coverage <- Sys.glob('../output/sarek/*/reports/mosdepth/*_tumor1/*_tumor1.md.mosdepth.summary.txt')
    input_CNA_stats <- Sys.glob('../output/QDNAseq/*/*/data/CNA_stats.txt')
    input_ACE_results <- Sys.glob('../output/ACE/*/*/ACE_fits.txt')
    input_CNH_results <- Sys.glob('../output/CNH/*/*/CNH_results.txt')
    input_TMB <- Sys.glob('../output/PureCN/*/*/*_tumor1_mutation_burden.csv')
    input_signatures <- Sys.glob('../output/PureCN/*/*/*_tumor1_signatures.csv')

    output_SampleData <- "sampledata/SampleData_WES.txt"
}
#-------------------------------------------------------------------------------
# 1.1 Read data 
#-------------------------------------------------------------------------------
# Read datasets
coverage <- tibble::tibble(sample = purrr::map(input_coverage, ~strsplit(.x,'/')[[1]][4]), data = purrr::map(input_coverage,read.delim)) %>% tidyr::unnest()

CNA_stats <- tibble::tibble(sample = purrr::map(input_ACE_results, ~strsplit(.x,'/')[[1]][5]),binsize = purrr::map(input_ACE_results, ~strsplit(.x,'/')[[1]][4]), data = purrr::map(input_CNA_stats,read.delim)) %>% tidyr::unnest()

ACE_results <- tibble::tibble(sample = purrr::map(input_ACE_results, ~strsplit(.x,'/')[[1]][5]),binsize = purrr::map(input_ACE_results, ~strsplit(.x,'/')[[1]][4]), data = purrr::map(input_ACE_results,read.delim)) %>% tidyr::unnest()

CNH_results <- tibble::tibble(sample = purrr::map(input_CNH_results, ~strsplit(.x,'/')[[1]][5]),binsize = purrr::map(input_CNH_results, ~strsplit(.x,'/')[[1]][4]), data = purrr::map(input_CNH_results,~read.delim(.x, sep =' '))) %>% tidyr::unnest()

TMB <- tibble::tibble(sample = purrr::map(input_TMB, ~strsplit(.x,'/')[[1]][5]),binsize = purrr::map(input_TMB, ~strsplit(.x,'/')[[1]][4]), data = purrr::map(input_TMB,~read.delim(.x, sep =','))) %>% tidyr::unnest()

signatures <- tibble::tibble(sample = purrr::map(input_signatures, ~strsplit(.x,'/')[[1]][5]),binsize = purrr::map(input_signatures, ~strsplit(.x,'/')[[1]][4]), data = purrr::map(input_signatures,~read.delim(.x, sep =','))) %>% tidyr::unnest()

#-------------------------------------------------------------------------------
# 1.2 Reformat and join data
#-------------------------------------------------------------------------------
# Calculate coverage
target_coverage <- coverage %>%
    filter(grepl("_region", chrom)) %>%  
    group_by(sample) %>%
    summarise(mean_target_depth = sum(bases) / sum(length) )

# Fetch best ACE fit
ACE_results <- ACE_results %>% group_by(sample,binsize) %>% filter(error == min(error)) %>% ungroup() %>% select(-error, - minimum)

# Join data
SampleData <- CNA_stats %>%
    left_join(target_coverage) %>%
    left_join(ACE_results) %>%
    left_join(CNH_results) %>%
    left_join(TMB) %>%
    left_join(signatures)

#-------------------------------------------------------------------------------
# 2.1 Write to file
#-------------------------------------------------------------------------------
write.table(SampleData, file =  output_SampleData, sep = '\t', quote = F, row.names = F)


library(ggplot2)
pdf('../plots/MINT_BarPlot_Mean_Target_depth.pdf', height = 4 , width =5)
SampleData %>%
  select(sample, mean_target_depth) %>%
  distinct() %>%
  ggplot(aes(forcats::fct_reorder(sample, mean_target_depth), mean_target_depth)) +
  geom_col(col = 'black',fill = "#4E79A7") +
  scale_y_log10(expand = expansion(mult = c(0, 0.1))) +
  geom_text(
    aes(label = round(mean_target_depth, 0)),
    vjust = -0.5,
    size = 4,
  ) +
    theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x='Sample',y = 'Mean target depth')
dev.off()


head(SampleData)
SampleData %>% filter(sample == 'MINT01_R2') %>% View()
read.delim()

