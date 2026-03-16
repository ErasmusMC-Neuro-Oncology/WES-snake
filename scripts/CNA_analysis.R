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
# 1) 
#
# History:
#  13-03-2026: File creation
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0.1  Load packages
#-------------------------------------------------------------------------------
suppressMessages(library(dplyr))
suppressMessages(library(CopywriteR))
suppressMessages(library(QDNAseq))
suppressMessages(library(Biobase))
suppressMessages(library(GenomicRanges))

#-------------------------------------------------------------------------------
# 0.2 Parse command line arguments
#-------------------------------------------------------------------------------
if(exists("snakemake")){
    input_bam <- snakemake@input[["bam"]]
    sample <- snakemake@wildcards[["sample"]]
    genome <- snakemake@params[["genome"]]
    binsize <- snakemake@wildcards[["binsize"]]
    cores <- snakemake@params[["cores"]]
    outdir <- snakemake@params[["outdir"]]
    sample_dir <- snakemake@output[["sample_dir"]]
    QDNAseq_output <- snakemake@output[["QDNAseq"]]
    Segments_output <- snakemake@output[["Segments"]]  
}else{
    input_bam <- 'output/sarek/MINT12/preprocessing/mapped/MINT12_tumor1/MINT12_tumor1.sorted.bam'
    sample <- 'MINT12_tumor1'
    sample_dir <- 'output/copywriter/MINT12_tumor1/'
    QDNAseq_output <- 'output/QDNAseq/MINT12/data/QDNAseq_Segments_100000bp.Rds'
    Segments_output <- 'output/QDNAseq/MINT12/data/QDNAseq_Segments_100000bp.Rds'
    genome <- 'hg38'
    binsize <- '100000'
    cores <- 10
    outdir <- 'output/copywriter/'
}

#-------------------------------------------------------------------------------
# 1.1 Define CopyWritR parameters
#-------------------------------------------------------------------------------
# Create annotation files
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
preCopywriteR(output.folder = outdir,
              bin.size = as.integer(binsize),
              ref.genome = genome,
              prefix = "chr")

# get number of kb bins
kbbin <- substring(binsize,1,nchar(binsize)-3)
# Load annotation files
load(file = file.path(outdir, paste0(genome,"_",kbbin,"kb_chr"), "blacklist.rda"))
# set number of workers
bp.param <- SnowParam(workers = cores, type = "SOCK")
# Create sample df
sample.control <- data.frame(samples = input_bam,controls=input_bam)

#-------------------------------------------------------------------------------
# 2.1 Run CopyWritR
#-------------------------------------------------------------------------------
# Run CopyWriteR
if(!"input.Rdata" %in% list.files(paste0(sample_dir,"/CNAprofiles/"))){
    unlink(paste0(sample_dir,"/CNAprofiles/"), recursive = TRUE)
    CopywriteR(sample.control = sample.control,
               destination.folder = sample_dir,
               reference.folder = file.path(outdir, paste0(genome,"_",kbbin,"kb_chr")),
               bp.param=bp.param)
}

#-------------------------------------------------------------------------------
# 3.1 Parse Copywriter output
#-------------------------------------------------------------------------------
read_counts <- read.delim(paste0(sample_dir,'/CNAprofiles/read_counts.txt'))

#-------------------------------------------------------------------------------
# 2.1 Prepare QDNAseq objects
#-------------------------------------------------------------------------------
#---------- create bins file  ----------
kbbin <- substring(binsize,1,nchar(binsize)-3)
load(paste0(outdir,genome,"_",kbbin,"kb_chr/GC_mappability.rda"))

# create dataframe containing fdata fields
fData_all <-
    cbind(as.data.frame(seqnames(GC.mappa.grange)),
          as.data.frame(ranges(GC.mappa.grange)),
          as.data.frame((mcols(GC.mappa.grange))))
# create features as row names
rownames(fData_all) <- paste0(fData_all$value,":",fData_all$start,"-",fData_all$end)

bins <- getBinAnnotations(as.integer(kbbin))
pData_bins <- pData(bins)
features_bins <- paste0('chr',rownames(pData_bins))
# remove weird features from bin data
pData_bins <- pData_bins[!rownames(pData_bins) %in% features_bins[!features_bins %in% rownames(fData_all)], ]
bins  <-pData_bins
#--------------------------------------------
#---------- create copynumber file ----------
colnames(read_counts)[6] <- 'read_counts'
read_counts$samples <- sample

counts  <-
    read_counts %>% 
    dplyr::select(Feature,read_counts,samples) %>%
    tidyr::spread(samples,read_counts,fill=0)

counts <- counts[match(paste0('chr',rownames(bins)), counts$Feature),] %>%
    dplyr::group_by(Feature) %>% dplyr::filter(dplyr::row_number() == 1 ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!is.na(Feature)) %>% 
    tibble::column_to_rownames(.,"Feature") %>%
    #tibble::add_column(bins = rownames(.), .before = colnames(.)[1]) %>%
    as.matrix()

#-------------------------------------------
#---------- create phenodata file ----------
phenodata <- data.frame(name = sample,
           total.reads = colSums(counts),
           used.reads = colSums(counts),
           row.names = sample,
           expected.variance = 0)

#--------------------------------------------------
#--------- create QDNAseqReadCounts object --------
rownames(bins) <- paste0('chr',rownames(bins))
features <- intersect(rownames(bins),rownames(counts))

QDNAseqCopyNumbers <- new("QDNAseqReadCounts",bins=bins[features,],counts=as.matrix(counts[features,]),phenodata=phenodata)

#-------------------------------------------------------------------------------
# 4.1 Perform QDNAseq normalizations
#-------------------------------------------------------------------------------
corrected <- applyFilters(QDNAseqCopyNumbers, residual=TRUE, blacklist=TRUE, mappability=FALSE, bases=FALSE , chromosomes=c('chrY','chrX','X','Y')) %>%
    estimateCorrection() %>%
    correctBins() %>%
    normalizeBins() %>%
    smoothOutlierBins() %>%
    segmentBins() %>%
    normalizeSegmentedBins()

#-------------------------------------------------------------------------------
# 5.1 Write to file
#-------------------------------------------------------------------------------
saveRDS(corrected,QDNAseq_output)

exportBins(corrected,'test.bed', format = 'bed', type = 'segments')
