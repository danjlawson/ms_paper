
source("ancestryanomaly.R")

#############
## Setup the data
## We need a matrix containing POPULATIONS in the rows and CHROMOSOMES in the columns
## Each entry is the location of a FILE, which contains THREE VALUES:
## The position, the average painting value, and the summed painting value
## For each SNP in the chromosome
## NOTE That for our paper, ancestry was stored as an integer 0-9 and therefore we normalise these numbers by dividing by 9. If your data were generated differently, you want to set the "norm" parameter of "importRowAndColSumsFile" to (1/your maximum value), probably 1. This is called by "importSnpSums" below
## e.g.
## position(hg19),average_painting,summed_painting
## 78015180,0.0190625,457.5
## ...

chromosomes=as.character(1:22)
donorpops=c("Yamnaya","EHG","WHG","CHG","Farmer","African","EastAsian")
ancestryfilemat=sapply(chromosomes,function(chr){
    infiles=paste0("../ukb_steppe_rowandcolsums/",
                   donorpops,".",chr,".sum_ancestry_per_snp.csv")
    names(infiles)=donorpops
    infiles
})
colnames(ancestryfilemat)=as.character(chromosomes)

## Import from a linearised form
## NB We need to know the number of haplotypes that the data are made up from.
## Here is where we give details of the file format, e.g.
## > importeddata=importSnpSums(ancestryfilemat,nhaps=24000,norm=1,rev=FALSE)
## would treat the snpmeans and snpsums as being generated on a (0-1) scale, and also
## as if the SNPs were in ascending position order. (chromopainter reverses them, so by default we assume you want to put them back).
importeddata=importSnpSums(ancestryfilemat,nhaps=24000)

## processing & annotation
## Either rerun analysis, or just load it
importeddataa<-try(readRDS("../ukb_steppe_rowandcolsums/allchr_12kinds_annotated.RDS"))
if(class(importeddataa)=="try-error"){
    importeddataa=lapply(importeddata,annotatepvals_singlepop)
    saveRDS(importeddataa,file="../ukb_steppe_rowandcolsums/allchr_12kinds_annotated.RDS")
}


## Construction of a single painting matrix
mypainting=makeMeanPainting(importeddataa)
myp=whitened_pval(mypainting)
myrmse2=unwhitened_rmse(mypainting)
mysnps=which(myp[,"p"]<5e-8)

####### Reporting genome-wide significance
ancestrychisqscore=getdataasranges(mysnps,importeddataa[[1]]$snpmat,myp[,"p"])
write.table(ancestrychisqscore,file="ukb_ancestry_deviations_total.csv")

####### Reporting population-specific significance
snpsdf=snpranges_asdf(importeddataa)
write.table(snpsdf,file="ukb_ancestry_deviations_bypopulation.csv")


###########################

## Basic plotting
for(d in donorpops) {
    png(paste0("ukb_allchr_test",d,".png"),height=1000,width=1000)
    signifplot(importeddataa[[d]],d)
    dev.off()
}
###

#######

mypt=myp[,"p"]
pthresh=1e-20
mypt[mypt<pthresh]=pthresh

png("Ancestry_deviation.png",height=1000,width=1000)
par(mfrow=c(2,1))
plot(importeddataa[[1]]$snpmat$pos,myrmse2,
     frame.plot=F,xlab="Genome Position",
     type="l",ylab="RMSE")
text(getgenomemids(importeddataa[[1]]$snpmat),
     rep(0.01 + max(myrmse2),
         length(unique((importeddataa[[1]]$snpmat[,"chromo"])))),
     labels=unique(importeddataa[[1]]$snpmat[,"chromo"]),
     adj=c(0.5,1))
abline(v=importeddataa[[1]]$gaps,col="grey",lwd=2)
## p-value
plot(importeddataa[[1]]$snpmat$pos,-log10(mypt),
     frame.plot=F,xlab="Genome Position",
     type="l",ylab="RMSE")
text(getgenomemids(importeddataa[[1]]$snpmat),
     rep(0.01 + max(-log10(mypt)),
         length(unique((importeddataa[[1]]$snpmat[,"chromo"])))),
     labels=unique(importeddataa[[1]]$snpmat[,"chromo"]),
     adj=c(0.5,1))
points(importeddataa[[1]]$snpmat$pos[myp<5e-8],
       -log10(mypt[myp<5e-8]),col=2,pch=19,cex=0.2)
points(importeddataa[[1]]$snpmat$pos[myp<pthresh],
       -log10(mypt[myp<pthresh]),col=2)
abline(v=importeddataa[[1]]$gaps,col="grey",lwd=2)
abline(h=-log10(5e-8))
dev.off()


