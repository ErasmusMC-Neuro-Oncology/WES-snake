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
# 1) Make robust for extracting sampleID and binsize
#
# History:
#  02-03-2026: File creation
#  29-05-2026: Add SigProfilerAssignment signatures
#  09-07-2026: Curate purities
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0.1  Load packages
#-------------------------------------------------------------------------------
suppressMessages(library(dplyr))

#-------------------------------------------------------------------------------
# 0.2 Parse command line arguments
#-------------------------------------------------------------------------------
if(exists("snakemake")){
    input_depth <- snakemake@input[["depth"]]
    input_CNA_stats <- snakemake@input[["CNA_stats"]]
    input_Segments <- snakemake@input[["Segments"]]
    input_CNH_results <- snakemake@input[["CNH_results"]]
    input_Purities <- snakemake@input[['Purities']]
    input_PureCN <- snakemake@input[["PureCN_purity"]]
    input_TMB <- snakemake@input[["TMB"]]
    input_variants <- snakemake@input[["variants"]]
    input_signatures <- snakemake@input[["signatures"]]
    output_SampleData <- snakemake@output[["SampleData"]]
}else{
    input_depth <- Sys.glob('/home/jurriaan/mnt/BIGR_home/SSLOWGRADE/output/WES/sarek/*/reports/mosdepth/*_tumor1/*_tumor1.md.mosdepth.summary.txt')
    input_CNA_stats <- Sys.glob('/home/jurriaan/mnt/BIGR_home/SSLOWGRADE/output/WES/QDNAseq/500kbp/*/data/CNA_stats.txt')
    input_Segments <- Sys.glob('/home/jurriaan/mnt/BIGR_home/SSLOWGRADE/output/WES/QDNAseq/500kbp/*/data/QDNAseq_Segments.txt')
    input_Purities<- Sys.glob('/home/jurriaan/mnt/BIGR_home/SSLOWGRADE/output/WES/PureCN_curated/500kbp/*/**_purity_bounds.csv')
    input_CNH_results <- Sys.glob('/home/jurriaan/mnt/BIGR_home/SSLOWGRADE/output/WES/CNH/500kbp/*/CNH_results.txt')
    input_PureCN <- Sys.glob('/home/jurriaan/mnt/BIGR_home/SSLOWGRADE/output/WES/PureCN_curated/500kbp/**/**_tumor1.csv')
    input_TMB <- Sys.glob('/home/jurriaan/mnt/BIGR_home/SSLOWGRADE/output/WES/PureCN_curated/500kbp/*/*_tumor1_mutation_burden.csv')
    input_variants <- Sys.glob('/home/jurriaan/mnt/BIGR_home/SSLOWGRADE/output/WES/PureCN_curated/500kbp/*/*_tumor1_variants.csv')
    input_signatures <- '/home/jurriaan/mnt/BIGR_home/SSLOWGRADE/output/WES/SigProfilerAssignment/Assignment_Solution/Activities/Assignment_Solution_Activities.txt'
    output_SampleData <- "/home/jurriaan/mnt/BIGR_home/SSLOWGRADE/output/WES/sampledata/SampleData_WES.txt"
}


#-------------------------------------------------------------------------------
# 1.1 Read data 
#-------------------------------------------------------------------------------
# Read datasets
depth <- tibble::tibble(sample = purrr::map(input_depth, ~strsplit(.x,'/')[[1]][9]), data = purrr::map(input_depth,read.delim)) %>% tidyr::unnest()

CNA_stats <- tibble::tibble(sample = purrr::map(input_CNA_stats, ~strsplit(.x,'/')[[1]][10]),binsize = purrr::map(input_CNA_stats, ~strsplit(.x,'/')[[1]][9]), data = purrr::map(input_CNA_stats,read.delim)) %>% tidyr::unnest()

CNH_results <- tibble::tibble(sample = purrr::map(input_CNH_results, ~strsplit(.x,'/')[[1]][10]),binsize = purrr::map(input_CNH_results, ~strsplit(.x,'/')[[1]][9]), data = purrr::map(input_CNH_results,~read.delim(.x, sep =' '))) %>% tidyr::unnest()

Purities <- tibble::tibble(sample = purrr::map(input_Purities, ~strsplit(.x,'/')[[1]][10]),binsize = purrr::map(input_Purities, ~strsplit(.x,'/')[[1]][9]), data = purrr::map(input_Purities,read.csv)) %>% tidyr::unnest()

PureCN <- tibble::tibble(sample = purrr::map(input_PureCN, ~strsplit(.x,'/')[[1]][10]),binsize = purrr::map(input_PureCN, ~strsplit(.x,'/')[[1]][9]), data = purrr::map(input_PureCN,~read.delim(.x, sep =','))) %>% tidyr::unnest() %>% rename(PureCN_ploidy = Ploidy, PureCN_purity = Purity)

TMB <- tibble::tibble(sample = purrr::map(input_TMB, ~strsplit(.x,'/')[[1]][10]),binsize = purrr::map(input_TMB, ~strsplit(.x,'/')[[1]][9]), data = purrr::map(input_TMB,~read.delim(.x, sep =','))) %>% tidyr::unnest()

variants <- tibble::tibble(sample = purrr::map(input_variants, ~strsplit(.x,'/')[[1]][10]),binsize = purrr::map(input_variants, ~strsplit(.x,'/')[[1]][9]), data = purrr::map(input_variants,~read.delim(.x, sep =','))) %>% tidyr::unnest()

signatures <- read.delim(input_signatures)

#signatures <- tibble::tibble(sample = purrr::map(input_signatures, ~strsplit(.x,'/')[[1]][9]),binsize = purrr::map(input_signatures, ~strsplit(.x,'/')[[1]][9]), data = purrr::map(input_signatures,~read.delim(.x, sep =','))) %>% tidyr::unnest()


#-------------------------------------------------------------------------------
# 1.2 Reformat and join data
#-------------------------------------------------------------------------------
# Calculate depth
target_depth <- depth %>%
    filter(grepl("_region", chrom)) %>%  
    group_by(sample) %>%
    summarise(mean_target_depth = sum(bases) / sum(length) )

# Fetch IDH mutation status and vaf
IDH_status <- variants %>%
  filter(ML.SOMATIC == TRUE) %>%
  group_by(sample, binsize) %>%
  summarise(
    IDHmt = case_when(
      any(gene.symbol == "IDH1") ~ "IDH1",
      any(gene.symbol == "IDH2") ~ "IDH2",
      TRUE ~ "IDHwt"
    ),
    IDHvaf_raw = {
      idh_vals <- AR[gene.symbol %in% c("IDH1", "IDH2")]
      if (length(idh_vals) > 0) max(idh_vals) else NA_real_
    }
  ) %>% ungroup()


mmr_signatures <- c('SBS6', 'SBS14', 'SBS15', 'SBS20', 'SBS21', 'SBS26', 'SBS44')
mmr_present    <- intersect(mmr_signatures, colnames(signatures))

# Keep absolute counts before normalizing
signatures_abs <- signatures %>%
    tibble::column_to_rownames('Samples')

# Normalize to relative contributions
signatures_rel <- signatures_abs %>%
    apply(1, function(x) x / sum(x)) %>%
    t() %>%
    as.data.frame()

# Build summary dataframe
signatures_summary <- signatures_rel %>%
    tibble::rownames_to_column(var = 'sample_short') %>%
    mutate(
        sample_short = gsub('_tumor1','',sample_short),
        TMZ_counts = signatures_abs[sample_short, 'SBS11'],
        MMR_counts     = rowSums(signatures_abs[sample_short, mmr_present, drop = FALSE]),
        MMR_rel     = rowSums(signatures_rel[sample_short, mmr_present, drop = FALSE]),
        MMR_active  = apply(signatures_abs[sample_short, mmr_present, drop = FALSE], 1,
                            function(x) paste(names(x[x > 0]), collapse = ','))) %>%
    left_join(CNA_stats %>% mutate(sample_short = purrr::map_chr(sample,~strsplit(.x, '\\.')[[1]][1])) %>% select(sample,sample_short)) %>%
    select(-sample_short)


# Join data
SampleData <- CNA_stats %>%
    left_join(target_depth) %>%
    left_join(CNH_results) %>%
    left_join(PureCN) %>%
    left_join(TMB) %>%
    left_join(signatures_summary) %>%
    left_join(IDH_status) %>%
    left_join(Purities)


# Rename and select columns
SampleData <- SampleData %>% rename(
                   CNH = Heterogeneity,
                   Nmutations = somatic.ontarget,
                   Codel_1p19 = chr1p19q_status,
                   TMB = somatic.rate.ontarget ) %>%
    select(sample,binsize,Sex,mean_target_depth,
           Nmutations, TMB, CNA_load,
           IDHmt, Codel_1p19, CDKN2AB_status,
           IDHvaf_raw,anchors_detail,IDH_purity,purity_1p19q,TP53_purity,anchor_consensus, ace_purity,ace_ploidy, PureCN_purity,PureCN_ploidy,
           TMZ_counts,contains('MMR'),
           contains('SBS'), Comment)


#-------------------------------------------------------------------------------
# 2.1 Write to file
#-------------------------------------------------------------------------------
write.table(SampleData, file =  output_SampleData, sep = '\t', quote = F, row.names = F)
