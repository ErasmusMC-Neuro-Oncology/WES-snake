#########################################################################/**
# @RdocFunction inputTextFile
#
# @title "Input text igv file input - return data frame of appropriate format"
#
# @synopsis
#
# \description{
#     @get "title".
# }
#
# \arguments{
#     \item{filename}{igv filename}
#     \item{type}{currently only supports igv}
# }



inputTextFile <- function(filename,type='igv'){
  input=read.table(filename,header=T,sep="\t",check.names=FALSE)
  input$feature=NULL
  colnames(input)[1:3] <- c("chrom","start","end")
  input = input[!(input$chrom %in% c("chrX","chrY","X","Y")),]

  return(input)
}



#########################################################################/**
# @RdocFunction inputTextFile        
#
# @title "Input chromosome borders"
#
# @synopsis
#
# \description{
#     @get "title".
# }
#
# \arguments{
#     \item{filename}{filename containing chromsome borders - only tested for human hg19}
# }


inputChromBorders <- function(filename="data/ChromInfoHG19.txt"){
  input=read.table(filename,stringsAsFactors=F)
  input=input[!(input$V1 %in% c("chrM","chrX","chrY")),]
  input =input[order(as.numeric(substring(input$V1,4))),]

  
  input$V3=cumsum(as.numeric(input$V2))
  return(input)
}

#########################################################################/**
# @RdocFunction inputQDNAseqObject
#
# @title "Input QDNAseq object"
#
# @synopsis
#
# \description{
#     @get "title".
# }
#
# \arguments{
#     \item{object}{QDNAseqobject containing segmented data}
#     \item{logTransform}{Boolean indicating whether data is log2 transformed}
#     \item{digits}{Number of digits to round segmented data - provided for compatibility with inputTextobject}
#     \item{chromosomeReplacements}{Indicating whether chromsomes should be replaced - by default replaces 23 by X, 24 by Y, and 25 by MT}
# }




inputQDNAseqObject <- function(object,  
    filter=TRUE, logTransform=TRUE, digits=3,
    chromosomeReplacements=c("23"="X", "24"="Y", "25"="MT"), ...) {


    if (inherits(object, "QDNAseqSignals")) {
        if (filter) {
            object <- object[QDNAseq:::binsToUse(object), ]
        }
        chromosome <- fData(object)$chromosome
        start <- fData(object)$start
        end <- fData(object)$end
        if (!"segmented" %in% assayDataElementNames(object))
            stop("Segments not found, please run segmentBins() first.")
        print(".") 
        dat <- assayDataElement(object, "segmented")
                   
        
        if (logTransform) {
            dat <- QDNAseq:::log2adhoc(dat)
        }
    }
    if (is.numeric(digits)) {
        dat <- round(dat, digits=digits)
    }
    oopts2 <- options(scipen=15)
    on.exit(options(scipen=oopts2), add=TRUE)
    out <- data.frame(chromosome=chromosome, start=start, end=end,
          dat, check.names=FALSE, stringsAsFactors=FALSE)
    
    return(out)
}

#########################################################################/**
# @RdocFunction inputQDNAseqObject
#
# @title "Input collection of bam files - returns QDNAseqObject with segmented data"
#
# @synopsis
#
# \description{
#     @get "title".
# }
#
# \arguments{
#     \item{dir}{directory containing bam files - QDNAseq will search directory recursively}
#     \item{binSize}{size of bins}
# }


inputBamFileDirectory <- function(dir,binSize){
  bins=QDNAseq:::getBinAnnotations(binSize=binSize)
  readCounts= QDNAseq:::binReadCounts(path=dir,bins=bins)
  readCountsFiltered <- QDNAseq:::applyFilters(readCounts,residual=TRUE, blacklist=TRUE)
  readCountsFiltered <- QDNAseq:::estimateCorrection(readCountsFiltered)
  copyNumbers <- QDNAseq:::correctBins(readCountsFiltered)
  copyNumbersNormalized <- QDNAseq:::normalizeBins(copyNumbers)
  copyNumbersSmooth <- QDNAseq:::smoothOutlierBins(copyNumbersNormalized)
  segmentedData <- QDNAseq:::segmentBins(copyNumbersSmooth,segmentStatistic='median',transformFun='sqrt')
  return(segmentedData)
}
