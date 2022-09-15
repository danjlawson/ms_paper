###################################
## Functional enrichment analysis
library("magrittr")
library("gprofiler2")
## if (!require("BiocManager", quietly = TRUE))
##     install.packages("BiocManager")
## BiocManager::install("SNPlocs.Hsapiens.dbSNP144.GRCh37")
library("SNPlocs.Hsapiens.dbSNP144.GRCh37")

getsnpsranges<-function(x){
    ## Turn a list with chrom and snpranges as elemnts into a character string vector chr:r1-r2
    paste0(x$chrom,":",x$p0snpranges[,1],"-",x$p0snpranges[,2])
}
getsnps<-function(x,snpmat){
    ## Get snps in the range format that SNPlocs likes, chrom:pos-pos
    tsm=snpmat[x$snps,]
    unique(paste0(tsm[,"chromo"],":",tsm[,"pos0"],"-",tsm[,"pos0"]))
}
enrichmenttest=function(feature,snpmat,retlist=FALSE){
    ## Takes a "feature" (a list containing at least snps=vector of indices)
    ## and a snpmat with chromo and pos0 features
    ## and then:
    ## a) extracts those snps
    ## b) looks them up in GHCh37 to get rsids
    ## c) looks up which gene they are in using gprofiler::gconvert
    ## d) unique
    ## e) computes their enrichmentusing gprofiler::gost
    require("SNPlocs.Hsapiens.dbSNP144.GRCh37")
    require("gprofiler2")
    if(length(feature$snps)==0) return(NULL)
    rsid=feature %>% 
        getsnps(snpmat=snpmat) %>%  GRanges %>%
        snpsByOverlaps(x=SNPlocs.Hsapiens.dbSNP144.GRCh37)
    genes=rsid$RefSNP_id %>% gconvert
    gunique = genes$name %>% unique
    gost=gost(gunique)
    if((length(gost$result)==0 ) || (dim(gost$result)[1]==0)) return(NULL)
    result=gost$result[order(gost$result$p_value),,drop=FALSE]
    if(retlist){
        return(list(result=result,
         gost=gost,
         rsid=rsid,
         genes=gunique))
    }
    return(result)
}

## Just a convenience function to get the useful columns from a gost object
gostsummary<-function(gostres,
                      cols=c("p_value","term_size","query_size","intersection_size","term_name"))
{gostres[,cols]}

## Look up the enrichment properties of the rmse object
rmseenrichment=enrichmenttest(list(snps=which(myp<5e-8)),
                              snpmat=importeddataa[[1]]$snpmat,retlist=TRUE)
rmseenrichment_nomhc=
    enrichmenttest(list(snps=which((myp<5e-8)&(importeddataa[[1]]$snpmat$chromo!="6"))),
                                    snpmat=importeddataa[[1]]$snpmat,retlist=TRUE)
## Write them to disk
write.csv(data.frame(rmseenrichment$rsid)[,c(1,2,4)],file="rmse_snps.csv")
write.table(rmseenrichment$genes,file="rmseenrichmentgenes.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(rmseenrichment_nomhc$genes,file="rmseenrichmentgenes_nomhc.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
rmseenrichment$result %>% gostsummary %>% head(n=20)
## You can put the rmseenrichment$genes in at https://biit.cs.ut.ee/gprofiler/gost to get a lovely plot

## Do the same for all SNPs that are listed in any ancestry pval_low or high
anyexcesslist=lapply(importeddataa,function(x)c(x$snpwindowsfdr$pval_low$snps,x$snpwindowsfdr$pval_high$snps))
anyexcess=do.call("c",anyexcesslist) %>% unique
anyenrichment=enrichmenttest(list(snps=anyexcess),
                              snpmat=importeddataa[[1]]$snpmat,retlist=TRUE)
write.table(anyenrichment$genes,file="anyenrichmentgenes.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

importeddataa<-try(readRDS("../ukb_steppe_rowandcolsums/allchr_12kinds_annotated.RDS"))


## Make a data frame of all our enrichment results
allenrichment=lapply(importeddataa,function(x){
    list(low=cbind(enrichmenttest(x$snpwindowsfdr$pval_low,
                            snpmat=x$snpmat),"Direction"="low"),
         high=cbind(enrichmenttest(x$snpwindowsfdr$pval_high,
                            snpmat=x$snpmat),"Direction"="high"))
})
allenrichmentrbound=sapply(allenrichment,function(x){
    keep=sapply(x,dim)[2,]>1
    x=x[keep]
    do.call("rbind",x)
})
allenrichmentrbound=lapply(names(allenrichmentrbound),function(x){
    if(!is.null(allenrichmentrbound[[x]])) return(cbind(allenrichmentrbound[[x]],"Ancestry"=x))
    return(NULL)
})
allenrichmentdf=do.call("rbind",allenrichmentrbound)
write.table(data.frame(allenrichmentdf)[c(2:13,15:16)],file="allenrichment_bypopulation.tsv",sep="\t")

#########################
#########################
#########################
## Make a data frame of all our SNP-specific results
allsnpdata=
do.call("cbind",
lapply(names(importeddataa),function(n){
    res=importeddataa[[n]]$snpmat[,4:6]
    colnames(res)=paste0(n,"_",colnames(res))
    if(which(names(importeddataa)==n)==1) {
        res = cbind(importeddataa[[n]]$snpmat[,c(1,3)],res)
    }
    res
})
)
allsnpdata=cbind(allsnpdata,data.frame("raw_rmse"=myrmse,
                                       "whitened_rmse"=myt,
                                       "rmse_pval"=mypt))
head(allsnpdata)
write.csv(allsnpdata,file="allsnpdata_annotated.csv")
system("gzip -f allsnpdata_annotated.csv")

#########################
#########################
#########################
#########################

LDA<-read.csv("../fdrtool results/LDA_score_annotated.csv.gz")
rownames(LDA)=paste0(LDA$chr,":",LDA$pd)
rownames(importeddataa[[1]]$snpmat)=paste0(importeddataa[[1]]$snpmat$chromo,":",importeddataa[[1]]$snpmat$pos0)
LDA=LDA[rownames(importeddataa[[1]]$snpmat),]
rownames(importeddataa[[1]]$snpmat)=1:dim(importeddataa[[1]]$snpmat)[1]

lda_high_enrichment_0.2=enrichmenttest(list(snps=which(LDA$pval_high<0.2)),
                              snpmat=importeddataa[[1]]$snpmat,retlist=TRUE)
lda_high_enrichment_0.5=enrichmenttest(list(snps=which(LDA$pval_high<0.5)),
                              snpmat=importeddataa[[1]]$snpmat,retlist=TRUE)
lda_high_enrichment_0.9=enrichmenttest(list(snps=which(LDA$pval_high<0.9)),
                              snpmat=importeddataa[[1]]$snpmat,retlist=TRUE)
lda_low_enrichment=enrichmenttest(list(snps=which(LDA$pval_low<0.2)),
                              snpmat=importeddataa[[1]]$snpmat,retlist=TRUE)
write.table(lda_low_enrichment$genes,file="lda_low_enrichmentgenes.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(lda_high_enrichment$genes,file="lda_high_enrichmentgenes.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

## https://biit.cs.ut.ee/gplink/l/6e9INLUPTx

#########################
#########################
#########################
####### NEED TO GET THE WHOLE SET, THIS IS CHR6 ONLY
mhc=read.csv("MHC_snps.csv")
mhcmatch=sapply(mhc[,1],function(x)
    which.min((importeddataa$Farmer$snpmat$pos0-x)^2)
    )
mhcl=lapply(fullresa,function(x){
    x$snpmat[mhcmatch,]
})

rmseenrichment=enrichmenttest(list(snps=which(myp<5e-8)),
                              snpmat=importeddataa[[1]]$snpmat,retlist=TRUE)
