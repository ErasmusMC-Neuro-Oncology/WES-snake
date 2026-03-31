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
suppressMessages(library(CopywriteR))
suppressMessages(library(QDNAseq))
suppressMessages(library(QDNAseq.hg38))
suppressMessages(library(Biobase))
suppressMessages(library(GenomicRanges))

#-------------------------------------------------------------------------------
# 0.2 Parse command line arguments
#-------------------------------------------------------------------------------
if(exists("snakemake")){
    input_bam <- snakemake@input[["bam"]]
    sample <- snakemake@wildcards[["patient"]]
    genome <- snakemake@params[["genome"]]
    binsize <- snakemake@wildcards[["binsize"]]
    cores <- snakemake@params[["cores"]]
    ACE_purity_penalty <- snakemake@params[["ACE_purity_penalty"]]
    ACE_ploidy_penalty <- snakemake@params[["ACE_ploidy_penalty"]]
    CNH_path <- snakemake@params[["CNH_path"]]
    cytobands <- snakemake@params[["cytobands"]]
    sample_dir <- snakemake@output[["sample_dir"]]
    QDNAseq_output <- snakemake@output[["QDNAseq"]]
    Profile_output <- snakemake@output[["Profile"]]  
    Segments_output <- snakemake@output[["Segments"]]
    CNA_stats_output <- snakemake@output[["CNA_stats"]]
    Segments_igv_output <- snakemake@output[["Segments_igv"]]
    Called_output <- snakemake@output[["Called"]]
    ACE_results_output <- snakemake@output[["ACE_results"]]
    ACE_matrix_output <- snakemake@output[["ACE_matrix"]]
    CNH_results_output <- snakemake@output[["CNH_results"]]
    CNH_plot_output <- snakemake@output[["CNH_plot"]]
    CNH_error_plot_output <- snakemake@output[["CNH_error_plot"]]
    
}else{
    input_bam <- '../output/sarek/MINT12/preprocessing/recalibrated/MINT12_tumor1/MINT12_tumor1.recal.bam'
    sample <- 'MINT12_tumor1'
    sample_dir <- '../output/copywriter/1000kbp/MINT12_tumor1/'
    QDNAseq_output <- '../output/QDNAseq/1000kbp/MINT12/data/QDNAseq_Segments.Rds'
    Segments_output <- '../output/QDNAseq/1000kbp/MINT12/data/QDNAseq_Segments.txt'
    ACE_purity_penalty <- 0
    ACE_ploidy_penalty <- 0.5
    CNH_path <- 'scripts/CNH/'

    Segments_igv_output <-paste0(getwd(), '/../output/QDNAseq/1000kbp/MINT12/data/QDNAseq_Segments.igv')
    Profile_output <- '../output/QDNAseq/1000kbp/MINT12/data/QDNAseq_Segments.Rds'
    CNH_results_output <- paste0(getwd(),'/../output/CNH/1000kbp/MINT12/CNH_results.txt')
    CNH_plot_output <- paste0(getwd(),'/../output/CNH/1000kbp/MINT12/CNH_plot.pdf')
    CNH_error_plot_output <- paste0(getwd(),'/../output/CNH/1000kbp/MINT12/CNH_error_plot.pdf')
    genome <- 'hg38'
    binsize <- '1000kbp'
    cores <- 1
    cytobands <- '/data/Resources/cytobands/hg38/cytoBand.txt'
    
}

#-------------------------------------------------------------------------------
# 1.1 Define CopyWritR parameters
#-------------------------------------------------------------------------------
binsize <- format(as.integer(gsub('kbp','',binsize))*1000, scientific = F)

# Create annotation files
if(!dir.exists(sample_dir)) dir.create(sample_dir, recursive = TRUE)
preCopywriteR(output.folder = sample_dir,
              bin.size = as.integer(binsize),
              ref.genome = genome,
              prefix = "chr")

# get number of kb bins
kbbin <- substring(binsize,1,nchar(binsize)-3)
# Load annotation files
load(file = file.path(sample_dir, paste0(genome,"_",kbbin,"kb_chr"), "blacklist.rda"))
# set number of workers
bp.param <- SnowParam(workers = cores, type = "SOCK")
# Create sample df
sample.control <- data.frame(samples = input_bam,controls=input_bam)

#-------------------------------------------------------------------------------
# 2.1 Run CopyWritR
#-------------------------------------------------------------------------------
# Run CopyWriteR
CopywriteR(sample.control = sample.control,
           destination.folder = sample_dir,
           reference.folder = file.path(sample_dir, paste0(genome,"_",kbbin,"kb_chr")),
           bp.param=bp.param)

#-------------------------------------------------------------------------------
# 3.1 Parse Copywriter output
#-------------------------------------------------------------------------------
read_counts <- read.delim(paste0(sample_dir,'/CNAprofiles/read_counts.txt'))
#-------------------------------------------------------------------------------
# 3.2 Prepare QDNAseq objects
#-------------------------------------------------------------------------------
#---------- create bins file  ----------
kbbin <- substring(binsize,1,nchar(binsize)-3)
load(paste0(sample_dir,'/',genome,"_",kbbin,"kb_chr/GC_mappability.rda"))

# create dataframe containing fdata fields
fData_all <-
    cbind(as.data.frame(seqnames(GC.mappa.grange)),
          as.data.frame(ranges(GC.mappa.grange)),
          as.data.frame((mcols(GC.mappa.grange))))
# create features as row names
rownames(fData_all) <- paste0(fData_all$value,":",fData_all$start,"-",fData_all$end)

bins <- getBinAnnotations(as.integer(kbbin), genome="hg38")

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
# 4.4 Run ACE
#-------------------------------------------------------------------------------
ACE_results <- ACE::squaremodel(corrected, QDNAseqobjectsample = T,
                                penalty = as.numeric(ACE_purity_penalty),
                                penploidy = as.numeric(ACE_ploidy_penalty)) 

# Save matrix plot
pdf(ACE_matrix_output)
ACE_results$matrixplot
dev.off()

#-------------------------------------------------------------------------------
# 4.5 Run CGH call
#-------------------------------------------------------------------------------
# Run CGHcall and extract calls
called <- callBins(corrected, nclass = 5) %>% CGHbase::calls()

# Save calls as df
called <- called %>% as.data.frame() %>% tibble::rownames_to_column()
colnames(called) <- c('bin','call')

#-------------------------------------------------------------------------------
# 4.6 Prepare annotations
#-------------------------------------------------------------------------------
cytobands <- read.delim(cytobands, header = F, col.names = c('chromosome','start','end','cytoband','staining'))

cytobands <-  cytobands %>%
    mutate(
        start = start + 1,
        arm = substr(cytoband, 1, 1)) %>%
    group_by(chromosome, cytoband,arm) %>%
    mutate(
    start = min(start),
    end   = max(end)) %>%
    filter(arm %in% c("p", "q"))

# annotate bins
called_annotated <- cbind(fData(corrected)[,c('chromosome','start','end')],called) %>%
    mutate(chromosome = paste0('chr',chromosome))  %>%
  inner_join(cytobands, by = c("chromosome" = "chromosome")) %>%
  filter(
    end.x >= start.y,
    start.x <= end.y
  )  %>%
    dplyr::rename(start = start.x, end = end.x) %>%
    select(chromosome,arm,cytoband,start,end,call) 


# calculate segment level calls
segment_calls <- cbind(fData(corrected)[,c('chromosome','start','end')],called) %>% filter(!is.na(call)) %>%
    mutate(chromosome = paste0('chr',chromosome)) %>% 
    group_by(chromosome) %>%
    mutate(
        group = cumsum(call != dplyr::lag(call, default = dplyr::first(call)))
    ) %>%
    group_by(chromosome, group) %>%
    summarise(
        loc.start = min(start),
        loc.end   = max(end),
        num.mark  = n(),
        call = as.numeric(names(sort(table(call), decreasing = TRUE)[1])),
        .groups = "drop")

# Fetch arm level calls: fetch the proportion
arm_level_calls <-
    called_annotated %>% filter(!is.na(call)) %>% group_by(chromosome,arm) %>%
    summarise(CNA_status = median(call), CNA_proportion = mean(call != 0))


#-------------------------------------------------------------------------------
# 4.6 Calculate stats for export
#-------------------------------------------------------------------------------
CNA_stats <- data.frame(
    CNA_load  =  sum(segment_calls$call != 0) / nrow(segment_calls), # proportion of non-neutral segments
    
    CDKN2AB_status =  called_annotated %>% filter(chromosome == 'chr9',cytoband == 'p21.3', start > 22000000) %>% pull(call) %>% median(), #CDKN2AB status
    
    chr1p19q_status =
        ((arm_level_calls %>% filter(chromosome == 'chr1',arm == 'p') %>% .$CNA_status < 0) & (arm_level_calls %>% filter(chromosome == 'chr1',arm == 'p') %>% .$CNA_proportion > 0.7)) & ((arm_level_calls %>% filter(chromosome == 'chr19',arm == 'q') %>% .$CNA_status < 0) & (arm_level_calls %>% filter(chromosome == 'chr19',arm == 'q') %>% .$CNA_proportion > 0.7)),
    
    chr7chr10_status =
        all(arm_level_calls[arm_level_calls$chromosome == 'chr7',]$CNA_status > 0) && all(arm_level_calls[arm_level_calls$chromosome == 'chr7',]$CNA_proportion > 0.7) &
        all(arm_level_calls[arm_level_calls$chromosome == 'chr10',]$CNA_status < 0) && all(arm_level_calls[arm_level_calls$chromosome == 'chr7',]$CNA_proportion > 0.7) 
)
#-------------------------------------------------------------------------------
# 4.6 Run CNH
#-------------------------------------------------------------------------------
exportBins(corrected, file = Segments_igv_output, format = 'igv')
system(paste0('cd ',CNH_path, ' ; Rscript R/Main.R ', Segments_igv_output,' 0.2 ',CNH_results_output,' ',CNH_plot_output,' ', CNH_error_plot_output))


#-------------------------------------------------------------------------------
# 5.1 Write to file
#-------------------------------------------------------------------------------
# Save QDNAseq object
saveRDS(corrected,QDNAseq_output)

# Save segments
write.table(Segments, Segments_output, sep = '\t',quote = F, row.names = F)

# Save calls
write.table(called_annotated, Called_output, sep = '\t',quote = F, row.names = F)

# Save ACE results
write.table(ACE_results$minimadf, ACE_results_output, sep = '\t',quote = F, row.names = F)

# Save CNA stats
write.table(CNA_stats, CNA_stats_output, sep = '\t',quote = F, row.names = F)

