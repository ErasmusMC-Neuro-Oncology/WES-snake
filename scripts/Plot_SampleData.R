#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Plot_SampleData.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Create SampleData plots
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
#  03-04-2026: File creation
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0.1  Load packages
#-------------------------------------------------------------------------------
suppressMessages(library(dplyr))
suppressMessages(library(stringr))

#-------------------------------------------------------------------------------
# 0.2 Parse command line arguments
#-------------------------------------------------------------------------------
if(exists("snakemake")){
    input_SampleData <- snakemake@input[["SampleData"]]
    output_depth <- snakemake@output[["Barplot_depth"]]
}else{
    input_SampleData <- "../output/sampledata/SampleData_WES.txt"
    output_depth <- "../output/plots/Barplot_target_depth.pdf"
}
#-------------------------------------------------------------------------------
# 1.1 Read data 
#-------------------------------------------------------------------------------
# Read datasets
SampleData <- read.delim(input_SampleData)

#-------------------------------------------------------------------------------
# 1.2 Plot data
#-------------------------------------------------------------------------------
Nsample <- length(unique(SampleData$sample))

# Create bar plot of coverage
pdf(output_depth, height = 4 , width = Nsample * 0.35)
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
