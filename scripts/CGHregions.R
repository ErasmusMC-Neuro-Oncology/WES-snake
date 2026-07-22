#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# CGHregions.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Run CGHregions on cellularity-based QDNAseq calls
#
# Author: Jurriaan Janssen (j.janssen.1@erasmusmc.nl)
#
# condaenv: CNA
# Usage: 
#
# TODO:
# 1) 
#
# History:
#  01-06-2026: File creation
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0.1  Load packages
#-------------------------------------------------------------------------------
if(!'BiocManager' %in% installed.packages()){install.packages('BiocManager',repos = "http://cran.us.r-project.org")}
if(!'CGHregions' %in% installed.packages()){BiocManager::install("CGHregions")}
suppressMessages(library(QDNAseq))
suppressMessages(library(CGHregions))
suppressMessages(library(dplyr))

#-------------------------------------------------------------------------------
# 0.2 Parse command line arguments
#-------------------------------------------------------------------------------
if(exists("snakemake")){
    input_QDNAseq <- snakemake@input[["QDNAseq"]]
    input_ACE_results <- snakemake@input[["ACE_results"]]
    input_cytobands <- snakemake@params[["cytobands"]]    

    output_Recalled <- snakemake@output[["Recalled"]]
    output_CGHregions <- snakemake@output[["CGHregions"]]
}else{
    input_QDNAseq <- Sys.glob('/home/jurriaan/mnt/BIGR_home/SSLOWGRADE/output/WES/QDNAseq/500kbp/*/data/QDNAseq_Segmented.Rds')
    input_ACE_results <- Sys.glob('/home/jurriaan/mnt/BIGR_home/SSLOWGRADE/output/WES/ACE/500kbp/*/ACE_fits.txt')
    input_cytobands <- '/home/jurriaan/mnt/BIGR_home/Resources/cytobands/hg38/cytoBand.txt'
   
}

#-------------------------------------------------------------------------------
# 1.1 Read data 
#-------------------------------------------------------------------------------
for(i in seq_along(input_QDNAseq)){
    object <- readRDS(input_QDNAseq[i])
    if(i == 1){
        QDNAseq_object <- object
    }else{
        QDNAseq_object <- Biobase::combine(QDNAseq_object,object)
    }
}


# Read ACE results
ACE_results <- tibble::tibble(sample = purrr::map(input_ACE_results, ~strsplit(.x,'/')[[1]][10]),binsize = purrr::map(input_ACE_results, ~strsplit(.x,'/')[[1]][9]), data = purrr::map(input_ACE_results,read.delim)) %>% tidyr::unnest() %>%
    group_by(sample) %>%
    filter(row_number() == 1)


#-------------------------------------------------------------------------------
# 2.1 Recall QDNAseq calls based on purity
#-------------------------------------------------------------------------------
purities <- ACE_results$cellularity[match(ACE_results$sample, Biobase::sampleNames(QDNAseq_object))]
recalled <- callBins(QDNAseq_object, nclass=3, cellularity=purities)


#-------------------------------------------------------------------------------
# 3.1 Run CGHregions
#-------------------------------------------------------------------------------
# Make CGH
CGH_recalled <- makeCgh(recalled)
# perform CGHregions
CGHregions_obj <- CGHregions(CGH_recalled)

#-------------------------------------------------------------------------------
# 3.2 Create CGHregions table
#-------------------------------------------------------------------------------
# read cytobands
cytobands <- read.delim(input_cytobands, header=F,stringsAsFactors=F)

addCytobands <- function(BED, cyto){
	colnames(cyto) <- c('chr', 'start', 'end', 'cytoband', 'bla')
	# convert 'chr1' to 1; and X into 23 and Y into 24
	new.chr <- c()
	for ( i in 1:nrow(cyto)){
	chr <- strsplit(cyto$chr[i], 'chr')[[1]][2]
	new.chr <- c(new.chr, chr)
	}
	y.ix <- which(new.chr=='Y')
	x.ix <- which(new.chr=='X')
	new.chr[y.ix] <- 24
	new.chr[x.ix] <- 23
	cyto$chr <- as.numeric(new.chr)

	## rename BED file columns:
	colnames(BED) <- c('chromo', 'start', 'end')

	##
	cytobands <- c()

	for ( i in 1:nrow(BED)){
		# get cytoband of region-start
		v1 <- BED$chromo[i] == cyto$chr
		v2 <- BED$start[i] > cyto$start
		v3 <- BED$start[i] < cyto$end
		cyto.start <- cyto$cytoband[which(v1&v2&v3)]

		# get cytoband of region-end
		v4 <- BED$chromo[i] == cyto$chr
		v5 <- BED$end[i] > cyto$start
		v6 <- BED$end[i] < cyto$end
		cyto.end <- cyto$cytoband[which(v4&v5&v6)]

		# add to vector with all cyto starts and ends
		region.cyto <- paste(cyto.start, cyto.end, sep='-')
		cytobands <- c(cytobands, region.cyto)
	}

	BED <- cbind(BED, cytobands)

}

regionsBED <- fData(CGHregions_obj)[,1:3]
regionsBED <- addCytobands(regionsBED, cytobands)

#-------------------------------------------------------------------------------
# 2.1 Write to file
#-------------------------------------------------------------------------------
# Write CGHregions
write.table(cbind(regionsBED, regions(CGHregions_obj)), file =  output_CGHregions, sep = '\t', quote = F, row.names = F)

# Write recalled object
saveRDS(recalled,output_Recalled)
