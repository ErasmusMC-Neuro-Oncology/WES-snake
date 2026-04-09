#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# CNA_analysis.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Run CopywriteR, Normalize with QDNAseq and export results
#
# Author: Jurriaan Janssen (j.janssen.1@erasmusmc.nl)
#
# condaenv: 
# Usage: 
#
# TODO:
# 1) Add PON corrections
#
# History:
#  13-03-2026: File creation
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0.1  Load packages
#-------------------------------------------------------------------------------
if(!'QDNAseq.hg38' %in% installed.packages()){devtools::install_github("asntech/QDNAseq.hg38@main")}
suppressMessages(library(dplyr))
suppressMessages(library(QDNAseq))
suppressMessages(library(QDNAseq.hg38))
suppressMessages(library(Biobase))
suppressMessages(library(GenomicRanges))

#-------------------------------------------------------------------------------
# 0.2 Parse command line arguments
#-------------------------------------------------------------------------------
if(exists("snakemake")){
    input_bam <- snakemake@input[["bam"]]
    genome <- snakemake@params[["genome"]]
    binsize <- snakemake@wildcards[["binsize"]]
    QDNAseq_output <- snakemake@output[["QDNAseq"]]
    Profile_output <- snakemake@output[["Profile"]]  
    Segments_output <- snakemake@output[["Segments"]]
        
}else{
    input_bam <- '../output/sarek/MINT12/preprocessing/recalibrated/MINT12_tumor1/MINT12_tumor1.recal.bam'
    QDNAseq_output <- '../output/QDNAseq/1000kbp/MINT12/data/QDNAseq_Segments.Rds'
    Segments_output <- '../output/QDNAseq/1000kbp/MINT12/data/QDNAseq_Segments.txt'
    Profile_output <- '../output/QDNAseq/1000kbp/MINT12/data/QDNAseq_Segments.Rds'
    genome <- 'hg38'
    binsize <- '1000kbp'
}

#-------------------------------------------------------------------------------
# 1.1 Read bam file
#-------------------------------------------------------------------------------
bins <- getBinAnnotations(as.integer(gsub('kbp','',binsize)), genome="hg38")


QDNAseqCopyNumbers <- binReadCounts(bins, bamfiles=input_bam, cache=TRUE)

#-------------------------------------------------------------------------------
# 2.1 Perform QDNAseq normalizations
#-------------------------------------------------------------------------------
corrected <- applyFilters(QDNAseqCopyNumbers, residual=TRUE, blacklist=TRUE, mappability=FALSE, bases=FALSE , chromosomes=c('chrY','chrX','X','Y')) %>%
    estimateCorrection() %>%
    correctBins() %>%
    normalizeBins() %>%
    smoothOutlierBins() %>%
    segmentBins() %>%
    normalizeSegmentedBins()

#-------------------------------------------------------------------------------
# 4.2 Plot QDNAseq profile and callBins
#-------------------------------------------------------------------------------
pdf(Profile_output, width = 6 , height = 5)
plot(corrected)
dev.off()

#-------------------------------------------------------------------------------
# 4.3 Retrieve segments
#-------------------------------------------------------------------------------
# Fetch segments and calculate their values
Segments <- fData(corrected) %>%
    mutate(seg.mean =assayData(corrected)$segmented[,1]) %>%
    filter(!is.na(seg.mean)) %>%
    arrange(chromosome, start) %>%
  group_by(chromosome) %>%
  mutate(group = cumsum(seg.mean != lag(seg.mean, default = dplyr::first(seg.mean)))) %>%
  group_by(chromosome, group) %>%
  summarise(
    loc.start = min(start),
    loc.end   = max(end),
    num.mark  = n(),
    seg.mean  = dplyr::first(seg.mean),
    .groups = "drop"
  ) %>%
    mutate(
        chrom_order = case_when(
            chromosome %in% c("X","x") ~ 23,
            chromosome %in% c("Y","y") ~ 24,
            TRUE ~ as.numeric(chromosome)),
        ID = sampleNames(corrected),
        chrom = paste0('chr',chromosome)) %>%
  arrange(chrom_order, loc.start) %>%
    select(ID, chrom, loc.start, loc.end, num.mark, seg.mean) 


#-------------------------------------------------------------------------------
# 5.1 Write to file
#-------------------------------------------------------------------------------
# Save QDNAseq object
saveRDS(corrected,QDNAseq_output)

# Save segments
write.table(Segments, Segments_output, sep = '\t',quote = F, row.names = F)


