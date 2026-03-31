library(Biobase)

source("R/input.R")

args=commandArgs(trailingOnly=TRUE)




source("R/plot.R")
source("R/EstimateHeterogeneity.R")


input=args[1]
MinPurity=args[2]
Outputfile=args[3]
Plotfile=args[4]
ErrorPlotfile=args[5]

df=inputTextFile(filename=args[1])

samplenames = colnames(df)[4]
res=rep(0,3)

Heterogeneity = rep(0,length(samplenames))
Ploidy=rep(0,length(samplenames))
Purity=rep(0,length(samplenames))

PuritySweep = seq(MinPurity,1.001,0.01)
PloidySweep = seq(2,4,0.01)

ChromBorders=inputChromBorders(filename="data/ChromInfoHG19.txt")
tmp=c(1,ChromBorders$V3)
ChromBorders$V1=as.numeric(substring(ChromBorders$V1,4))


for(i in 4:ncol(df)){
  df$segments = df[,i]
  x=DataFrameToSegmentWeights(df)

  for(j in 1:nrow(x)){x$absEnd[j]=tmp[x$chrom[j]]+x$endCoordinates[j]}
  for(j in 1:nrow(x)){x$absStart[j]=tmp[x$chrom[j]]+x$startCoordinates[j]}


  r=EstimateHeterogeneity(x, Ploidy=PloidySweep,Purity = PuritySweep)
  res=ComputeMinimum(r=r,Ploidy=PloidySweep,Purity=PuritySweep)
  Heterogeneity[i-3]=res[[1]]
  Ploidy[i-3]=res[[3]]
  Purity[i-3]=res[[2]]
  plotSample(x,filename=Plotfile,reScale=TRUE,Purity=Purity[i-3],Ploidy=Ploidy[i-3],ChromBorders,main=colnames(df)[i],Heterogeneity=Heterogeneity[i-3])
  plotError(r=res[[4]],Ploidy=PloidySweep,Purity=PuritySweep,filename=ErrorPlotfile,Title=colnames(df)[i])
}

#write.table(format(data.frame(samplenames,Heterogeneity,Ploidy,Purity),digits=3,nsmall=4),"testCases/fitting_sc_Ploidy_curated.txt",quote=F)
write.table(format(data.frame(samplenames,Heterogeneity,Ploidy,Purity),digits=3,nsmall=4),file=Outputfile,quote=F)
