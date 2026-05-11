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
suppressMessages(library(ggplot2))


#-------------------------------------------------------------------------------
# 0.2 Parse command line arguments
#-------------------------------------------------------------------------------
if(exists("snakemake")){
    input_SampleData <- snakemake@input[["SampleData"]]
    output_depth <- snakemake@output[["Barplot_depth"]]
}else{
    input_SampleData <- "../output/WES/sampledata/SampleData_WES.txt"
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



function(x){
clinical <- read.delim('~/mnt/BIGR_home/SSLOWGRADE/data/clinical/Merged_clinical_data.csv', sep = ',') 

SampleData <- SampleData %>%
    left_join(clinical, by = c('sample' = 'patient'))

library(ComplexHeatmap)

# Example data
samples <- forcats::fct_reorder(SampleData$sample, SampleData$TMB)
values <- SampleData$TMB

library(circlize)

SampleData %>% filter(sample == 'SG_022') %>% View()


ha_top <- HeatmapAnnotation(
  `TMZ signature (SBS11)` = SampleData$SBS11,
  col = list(
    `TMZ signature (SBS11)` = circlize::colorRamp2(
      range(SampleData$SBS11, na.rm = TRUE),
      c("white", "#E15759")
    )
  ),

  rect_gp = grid::gpar(col = "white", lwd = 0.8),

  annotation_legend_param = list(
    `TMZ signature (SBS11)` = list(direction = "horizontal")
  )
)

ha_points <- HeatmapAnnotation(
    `Mutations per Mb` = anno_points(
        values,
        pch = 21,
        size = unit(2, "mm"),
        gp = grid::gpar(col = "white", fill = "#4E79A7", lwd = 0.5),
        ylim = c(0, max(values, na.rm = TRUE)),

        panel.fun = function(index, nm) {
            grid::grid.yline(
                      y = 6,
                      gp = grid::gpar(col = "black", lwd = 1, lty = 2)
                  )
        }
    ),
    height = unit(3, "cm")
)



ht <- Heatmap(
  matrix(values, nrow = 1),
  col = c("white", "white"),  # <- trick: constant color
  rect_gp = grid::gpar(col = NA),
  height = unit(1, "mm"),
  show_heatmap_legend = FALSE,
  top_annotation = ha_top,
  bottom_annotation = ha_points
)

pdf('TMB_versus_TMZsignature.pdf', height = 3, width = 6)
draw(ht,annotation_legend_side = "top")
dev.off()





# Main data
df <- data.frame(
  sample = SampleData$sample,
  value = SampleData$TMB
)

# Annotation data (same x, different "track")
annot <- data.frame(
  sample = df$sample,
  TMZ = SampleData$Treatment,
  TMZ_signature = SampleData$SBS11,
  SBS6 = SampleData$SBS6,
  SBS14 = SampleData$SBS14,
  SBS15 = SampleData$SBS15,
  SBS20 = SampleData$SBS20,
  SBS21 = SampleData$SBS21,
  SBS26 = SampleData$SBS26,
  SBS44 = SampleData$SBS44)

# Convert to long format
library(tidyr)
annot_long <- pivot_longer(annot, -sample, names_to = "type", values_to = "value")

# Assign y positions for annotation tracks
annot_long$y <- as.numeric(factor(annot_long$type)) + 1.5 * 100
annot_long

ggplot(df, aes(x = sample, y = value)) +
  geom_point() +
  geom_tile(
    data = annot_long,
    aes(x = sample, y = y, fill = value),
    height = 100
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.3)))

SampleData %>%
  select(sample, TMB) %>%
  distinct() %>%
  ggplot(aes(forcats::fct_reorder(sample, TMB), TMB)) +
  geom_col(col = 'black',fill = "#4E79A7") +
  scale_y_log10(expand = expansion(mult = c(0, 0.1))) +
  geom_text(
    aes(label = round(TMB, 0)),
    vjust = -0.5,
    size = 4,
  ) +
    theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x='Sample',y = 'Mean target depth')
colnames(SampleData)


SampleData %>%
  select(sample, TMB) %>%
  distinct() %>%
  ggplot(aes(forcats::fct_reorder(sample, TMB), TMB)) +
  geom_col(col = 'black',fill = "#4E79A7") +
  scale_y_log10(expand = expansion(mult = c(0, 0.1))) +
  geom_text(
    aes(label = round(TMB, 0)),
    vjust = -0.5,
    size = 4,
  ) +
    theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x='Sample',y = 'Mean target depth')
colnames(SampleData)
}
