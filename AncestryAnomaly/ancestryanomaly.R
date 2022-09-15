################################
### READING THE NEW WAY
importRowAndColSumsFile<-function(ancestrypersnpfile,nhaps,norm=1/9,rev=TRUE){
    ## Read one file and return it in our format
    res=read.csv(ancestrypersnpfile)
    if(rev){
        index=dim(res)[1]:1
    }else{
        index=1:dim(res)[1]
    }
    list(snpnames=res[index,1],
         snpmeans=res[index,2]*norm,
         snpsums=res[index,3]*norm,
         nhaps=nhaps)
}

importSnpSums=function(ancestryfilemat,nhaps,...){
    ## Read a matrix of files and return it in our format
    fullres=lapply(rownames(ancestryfilemat),function(ancestry){
        chrdata=lapply(ancestryfilemat[ancestry,],importRowAndColSumsFile,nhaps=nhaps,...)
        return(list(chrdata=chrdata))
    })
    names(fullres)=rownames(ancestryfilemat)
    fullres
}

getSignifSNPs<-function(mysnps){
    ## Takes a vector containing indices
    ## And returns a matrix containing the start and end index
    ## of each contiguous region
  mysnps<-c(mysnps,tail(mysnps,1))
  ret<-matrix(0,ncol=2,nrow=0)
  tdiff<-diff(mysnps)
  while(length(tdiff)>0){
    tnext<-(which(tdiff>1))[1]
    if(is.na(tnext)) tnext<-length(tdiff)
    ret<-rbind(ret,as.numeric(c(mysnps[1],mysnps[tnext])))
    tdiff<-tdiff[-(1:(tnext))]
    mysnps<-mysnps[-(1:(tnext))]
  }
  ret
}

getdataasranges<-function(snps,snpmat,values,summary="min",
                          chromosome="chromo",position="pos0"){
    ## Take a set of snps identified by their row index
    ## and a snp matrix with chromosome info in the chromosome feature and position information in the position feature
    ## and extract it as a range, summarising the data by get(summary)
    mysnpranges=getSignifSNPs(snps)
    ret=do.call("rbind",
                lapply(1:dim(mysnpranges)[1],function(i){
                    x=mysnpranges[i,]
                    return(data.frame(
                        chromosome=snpmat[[chromosome]][x[1]],
                        left=snpmat[[position]][x[1]],
                        right=snpmat[[position]][x[2]],
                        summary=get(summary)(values[x[1]:x[2]]))
                        )
                })
    )
    colnames(ret)[4]=summary
    ret
}

snpranges_asdf_forpop<-function(x,who,what="fdr"){
    ## Return the SNP ranges for a particular population as a data frame
    xx=x[[paste0("snpwindows",what)]]
    lowdf=data.frame(chromo=xx$pval_low$chrom,
               left=xx$pval_low$p0snpwindow[,1],
               right=xx$pval_low$p0snpwindow[,2],
               minp=xx$pval_low$minpval,
               direction=rep("low",length(xx$pval_low$chrom)))
    highdf=data.frame(chromo=xx$pval_high$chrom,
               left=xx$pval_high$p0snpwindow[,1,drop=F],
               right=xx$pval_high$p0snpwindow[,2,drop=F],
               minp=xx$pval_high$minpval,
               direction=rep("high",length(xx$pval_high$chrom)))
    res=cbind(rbind(lowdf,highdf),
              population=rep(who,dim(lowdf)[1]+dim(highdf)[1]))
    res
}

snpranges_asdf<-function(x,what="fdr"){
    ## Return the SNP ranges for all populations as a data frame
    do.call("rbind",
            lapply(names(x),function(n)snpranges_asdf_forpop(x[[n]],n))
            )
}

getsnppos<-function(files){
    ## Return the SNP positions (header) for each file in files
    ret=lapply(files,function(x)
        as.character(read.table(x,nrows=1,row.names=1)[1,]))
    names(ret)=names(files)
    ret
}
getnsnps<-function(filemat){
    ## Count the number of SNPs for every file in a matrix of files
    ## Where each row is a different population with the same assumed number of SNPs
    ## Also accepts a single file
    if((class(filemat)=="character") & (length(filemat)==1)) return(length(getsnppos(filemat)[[1]]))
    snppos=getsnppos(filemat[1,])
    sapply(snppos,length)
}
getnhaps<-function(file){
    ## Count the number of haplotypes in a haplotype-formatted file
    ## This is simply the number of lines (excluding the end one)
   length(count.fields(file, sep = "\n"))-1
}
processFile<-function(file,nhaps=NULL,nsnps=NULL,scale=1,verbose=1){
    ## Read haplotypes format files using the fast cpp wrapper
    require("Rcpp")
    sourceCpp("rcppreadfile.cpp")
    if(all(is.null(nhaps))) {
        if(verbose>0) print(paste("Counting haplotypes..."))
        nhaps=getnhaps(file)
    }
    if(all(is.null(nsnps))) {
        if(verbose>0) print(paste("Counting SNPs..."))
        nsnps=getnsnps(file)
    }
    if(verbose>0) print(paste("Computing genome summaries..."))
    res=readfilecpp(file,nhaps,nsnps,scale,verbose)
    res$snpnames=as.numeric(res$snpnames)
    res$nsnps=nsnps
    res$nhaps=nhaps
    res
}
getmeanpainting<-function(filemat,nsnps=NULL,nhaps=NULL,scale=1/9,verbose=1,...){
    ## Run the data extraction on all files in the file matrix
    if(all(is.null(nhaps))) {
        if(verbose>0) print(paste("Counting haplotypes..."))
        nhaps=getnhaps(filemat[1,1])
    }
    if(all(is.null(nsnps))) {
        if(verbose>0) print(paste("Counting SNPs..."))
        nsnps=getnsnps(filemat)
    }
    rres=lapply(1:dim(filemat)[1],function(pop){
        chrres=lapply(1:dim(filemat)[2],function(chr){
            if(verbose>0) print(paste("Running chromosome",chr,
                                      "of",dim(filemat)[2],
                                      "with population",pop,"of",dim(filemat)[1]))
            processFile(filemat[pop,chr],
                        nhaps,nsnps[chr],
                        scale=scale,verbose=(verbose>1))
        })
        names(chrres)=colnames(filemat)
        list(chrdata=chrres)
    })
    names(rres)=rownames(filemat)
    return(rres)
}
oldreaddata=function(filemat){
    require("data.table")
    ## Suitable for testing with the old method from Nelson et al 2017
    alllist=list()
    alllist[[1]]=
        lapply(1:dim(filemat)[1],function(i){
            r=data.table::fread(filemat[i,1],header=T)
            rn=r[,1]
            r=as.matrix(r[,-1])
            r=t(r)
            colnames(r)=paste0(as.character(as.data.frame(rn)[,1]),c(".A",".B"))
            (r/9)[dim(r)[1]:1,]
        })
    names(alllist[[1]])=donorpops
    alllist
}

computesums = function(filepath,nsnps,nhaps,scale=1/9, sep=" ",verbose=T) {
    ## R approach: slow so this is for reference only
    sumsnps=rep(0,nsnps)
    sumhaps=rep(0,nhaps)
    con = file(filepath, "r")
    ## Skip the header
    myLine=scan(con,what="character",nlines=1,sep=sep,skip=0,quiet=TRUE)
    hapon=0
    while (length(
        myLine <- scan(con,what="character",nlines=1,sep=sep,skip=0,quiet=TRUE)
    ) > 0 ){
        hapon<-hapon + 1
        vals=as.numeric(myLine[-1])*scale
        sumhaps[hapon]=sum(vals)
        sumsnps=sumsnps + vals
    }
    close(con)
    return(list(sumsnps=rev(sumsnps),
                sumhaps=sumhaps,
                nsnps=nsnps,
                nhaps=nhaps))
} 
################################
################################
## Modelling
##################
## Poisson-Binomial appropriate for Hapmix/Mosaic
approxpars_poibin<-function(hapsums,nsnps){
    ## Takes a list of sums for each haplotype, and the corresponding number of snps
    ## Uses this to get a genome wide estimate of the ancestry distribution
    hapmat=do.call("rbind",hapsums)
    hapmeans=colSums(hapmat)/sum(nsnps)
    approxpars_poibin_simple(hapmeans)
}
approxpars_poibin_simple<-function(hapmeans){
    ## poibin has:
    ## mean = \sum_{i=1}^n p_i
    ## var = \sum_{i=1}^n p_i (1 - p_i) =  \sum_{i=1}^n p_i - p_i^2
    ## This matches that with the normal approximation to the poibin
    list(mu=sum(hapmeans),sigma=sqrt(sum(hapmeans*(1-hapmeans))))
}
approxpars_blockmedian_windows<-function(snpsums,nhaps,pos,windowsize=5e6){
    ## Gets a dataframe of the mean and sd within each window
    x=snpsums/nhaps
    posmin=pos[1]
    td=data.frame(tpos=pos-posmin+1,val=x)
    nwindows=ceiling(diff(range(pos))/windowsize)
    sums=rep(0,nwindows)
    sumsqs=rep(0,nwindows)
    counts=rep(0,nwindows)
    td$window=as.factor(ceiling((td$tpos)/windowsize))
    windowsds=sapply(levels(td$window),function(l){
        sd(td[td[,"window"]==l,"val"])
    })
    windowmeans=sapply(levels(td$window),function(l){
        mean(td[td[,"window"]==l,"val"])
    })
    data.frame(mean=windowmeans,sd=windowsds)
}

approxpars_blockmedian<-function(snpsums,nhaps,pos,windowsize=5e6){
    ## Block-Median
    ## Compute the mean and standard deviation of SNP paintings based on a median of windows approach
    allwindows=do.call("rbind",
            lapply(1:length(snpsums),function(chr){
        approxpars_blockmedian_windows(snpsums[[chr]],nhaps[chr],pos[[chr]],windowsize)
    }))
    return(list(mu=median(allwindows[,"mean"]),
                sigma=median(allwindows[,"sd"])))
}

###################
approxpoibin=function(tx,x){
    ## Just for testing
    ## Mean = sum_i p_i
    ## Var = sum_i p_i (1-p_i)
    mu=sum(x)
    sigma=sqrt(sum(x*(1-x)))
    data.frame(x=tx,y=dnorm(tx,mu,sigma))
}
###################
computepvals=function(x,pars,what="means"){
    ## Annotate a result by whether the value is high or low compared to the expected distribution
    x$pars=pars
    x$snpmeans=x$snpsums/x$nhaps
    if(what=="means"){
        x$pval_low=pnorm(x$snpmeans,x$pars$mu,x$pars$sigma,lower=T)
        x$pval_high=pnorm(x$snpmeans,x$pars$mu,x$pars$sigma,lower=F)
    }else{
        x$pval_low=pnorm(x$snpsums,x$pars$mu,x$pars$sigma,lower=T)
        x$pval_high=pnorm(x$snpsums,x$pars$mu,x$pars$sigma,lower=F)
    }
    x
}
gwsignificance<-function(signifmat,p=0.05){
    ## Computes the genome-wide significance assuming all SNPs are indepdendent
    apply(signifmat,2,function(x){
        p/length(na.omit(x))
    })
}
gweffective<-function(signifmat,pgw,lag.max=60){
    ## Computes the effective genome-wide significance threshold, accounting for correlation between snps
    sapply(1:dim(signifmat)[2],function(i){
        x=na.omit(signifmat[,i])
        ttacf<-acf(x,lag.max = lag.max,plot=FALSE)
        pgw[i]*(1+ttacf$acf[2])/(1-ttacf$acf[2])
    })
}
################
## PLOTTING THINGS
getgenomegaps<-function(snpmat){
    ## Identify where there are changes in the chromosome
    tt=as.numeric(as.factor(snpmat[,"chromo"]))
    cgaps<-which(diff(tt)>0)
    cgappos<-(snpmat[cgaps,"pos"]+snpmat[cgaps+1,2])/2
    cgappos
}
getgenomemids<-function(snpmat){
    ## Get the midpoint of each chromosome
    chromotextpos<-sapply(
        unique(snpmat[,"chromo"]),
        function(x){
            trange<-range(snpmat[snpmat[,"chromo"]==x,"pos"])
            mean(trange)
        })
    chromotextpos
}

getSnpRangesPval<-function(pvalsall,tthresh,snpwindow=500,offset=3){
    ## Extract extreme pvalue ranges from a list of pvals
    ## pvalsall is a matrix of pvalues (each column treated independently)
    ## tthresh is a vector of thresholds for each, for grouping
    ## snpwindow: SNPs that are within a certain window, but fall in different threshold bins, can be merged if they are closer than this
    ## offset: we ignore the first offset columbs of pvalsall
  nppops<-dim(pvalsall)[2]-offset
  ret<-list()
  if(length(tthresh)!=dim(pvalsall)[2]-offset) stop("Must provide 1 p-value threshold per population column")
  for(popon in 1:nppops){
    ret[[popon]]<-list()
    tw<-pvalsall[,popon+offset]<=tthresh[popon]
    tw[is.na(tw)]=FALSE
    if(any(tw)) {
      ret[[popon]][["snps"]]<-which(tw)
      ret[[popon]][["snpranges"]]<-getSignifSNPs(ret[[popon]][["snps"]])
      ret[[popon]][["chrom"]]<-apply(ret[[popon]][["snpranges"]],1,function(x){
        pvalsall[x[1],"chromo"]
      })
      ret[[popon]][["snpwindow"]]<-t(apply(
          ret[[popon]][["snpranges"]],1,
          function(x){
              tchrom<-pvalsall[x[1],"chromo"]
              c(max(x[1]-snpwindow,min(which(pvalsall[,"chromo"]==tchrom))),
                min(x[2]+snpwindow,max(which(pvalsall[,"chromo"]==tchrom))))
          }))
      ret[[popon]][["psnpranges"]]<-matrix(apply(ret[[popon]][["snpranges"]],2,function(x){pvalsall[x,"pos"]}),ncol=2)
      ret[[popon]][["psnpwindow"]]<-matrix(apply(ret[[popon]][["snpwindow"]],2,function(x){pvalsall[x,"pos"]}),ncol=2)
      ret[[popon]][["p0snpranges"]]<-matrix(apply(ret[[popon]][["snpranges"]],2,function(x){pvalsall[x,"pos0"]}),ncol=2)
      ret[[popon]][["p0snpwindow"]]<-matrix(apply(ret[[popon]][["snpwindow"]],2,function(x){pvalsall[x,"pos0"]}),ncol=2)
      ret[[popon]][["minpval"]]<-apply(ret[[popon]][["snpranges"]],1,function(x){
        min(pvalsall[x[1]:x[2],popon+offset])
      })
    }else{
        ret[[popon]][["p0snpwindow"]]<-ret[[popon]][["p0snpranges"]]<-ret[[popon]][["psnpwindow"]]<-ret[[popon]][["psnpranges"]]<-ret[[popon]][["snpwindow"]]<-ret[[popon]][["snpranges"]]<-matrix(nrow=0,ncol=2)
        ret[[popon]][["snps"]]<-ret[[popon]][["minpval"]]<-ret[[popon]][["chrom"]]<-numeric()
    }
  }
  names(ret)<-colnames(pvalsall)[-(1:offset)]
  ret
}

getlfdr<-function(signifmat,fdrmax=1-1e6,upperthreshold=1e-20){
    ## apply false discovery rate correction to the significance matrix.
    ## fdrmax is the significance threshold for reporting
    ## upperthreshold is an upper bound on extremely significant results, which cause numerical issues if not capped
  require("fdrtool")
  mykeep<-which(!is.na(signifmat[,1]))
  signifmat[signifmat<upperthreshold]=upperthreshold
  fdrvals<-apply(signifmat[mykeep,],2,function(x) {
    x[x<0]<-0
    x[x>1]<-1
    tsignif=data.frame(index=numeric(),lfdr=numeric(),z=numeric(),p=numeric())

    myz<- qnorm(x)
    tmp<-fdrtool(myz,plot=FALSE)
    default=list(signif=tsignif,lfdr=tmp$lfdr,lfdrmax=1,pmax=0,n50=0,fdrtool=tmp)
    tdull<-which(myz>0 & tmp$lfdr<1)
    tmine<-which(myz<0 & tmp$lfdr<1)
    if(length(tmine)==0) return(default)
    myfrac<-length(tmine)/(length(tdull)+length(tmine))

    n50<-floor(myfrac*as.numeric(floor((1-tmp$param[1,"eta0"])/2*length(tmp$pval))))

    if(n50==0) return(default)

    tmp$lfdr[myz>0]<-1
    torder<-order(myz,decreasing=F)
    tlfdr<-tmp$lfdr[torder[1:n50]]
    tsignif<-data.frame(index=torder[1:n50],lfdr=tmp$lfdr[torder[1:n50]],z=myz[torder[1:n50]],p=x[torder[1:n50]])
    tsignif<-tsignif[tsignif$lfdr<fdrmax,,drop=FALSE]
    n50<-dim(tsignif)[1]
    if(n50==0) return(default)

    list(signif=tsignif,lfdr=tmp$lfdr,lfdrmax=max(tsignif$lfdr),pmax=max(tsignif$p),n50=n50,fdrtool=tmp)
  })

  pthresh<-sapply(fdrvals,function(x){max(x$pmax)})
  n50<-sapply(fdrvals,function(x){x$n50})
  maxlfdr<-sapply(fdrvals,function(x){
    if(dim(x$signif)[1]==0)return(1);
    max(x$signif$lfdr) })
  minlfdr<-sapply(fdrvals,function(x){
    if(dim(x$signif)[1]==0)return(1);
    min(x$signif$lfdr)})
  list(fdrvals=fdrvals,thresh=pthresh,n50=n50,minlfdr=minlfdr,maxlfdr=maxlfdr)
}

annotatepvals_singlepop<-function(respoplist,fdrmax=0.2,
                                  chromosomegap=100000,
                                  snpwindow=200,pgw=0.05,
                                  lag.max=60,windowsize=5e6){
    ## Compute the parameters of the normal distribution fit
    ## The Poisson-Binomial approximation
    ## hapsums=lapply(respoplist$chrdata,function(x){
    ##     x$hapsums
    ## })
    ## nsnps=sapply(respoplist$chrdata,function(x){
    ##     x$nsnps
    ## })
    ## respoplist$pars=approxpars_full(hapsums,nsnps)
    ## The Block-median approximation
    snpsums=lapply(respoplist$chrdata,function(x){
        x$snpsums
    })
    nhaps=sapply(respoplist$chrdata,function(x){
        x$nhaps
    })
    pos=lapply(respoplist$chrdata,function(x){
        x$snpnames
    })
    respoplist$pars=approxpars_blockmedian(snpsums,nhaps,pos,windowsize)
    ## Create a chromosome-flattened matrix of significance values
    respoplist$chrdata=lapply(respoplist$chrdata,computepvals,pars=respoplist$pars)
    offsets = cumsum(c(0,chromosomegap +
                         sapply(respoplist$chrdata,function(x)max(x$snpnames))))
    respoplist$snpmat=do.call("rbind",
            lapply(1:length(respoplist$chrdata),function(i){
                data.frame(
                    chromo=names(respoplist$chrdata)[i],
                    pos=respoplist$chrdata[[i]]$snpnames+offsets[i],
                    pos0=respoplist$chrdata[[i]]$snpnames,
                    pval_low=respoplist$chrdata[[i]]$pval_low,
                    pval_high=respoplist$chrdata[[i]]$pval_high,
                    mean=respoplist$chrdata[[i]]$snpmeans)
            }))
    respoplist$gaps=getgenomegaps(respoplist$snpmat)
    ## Compute FDR thresholds
    respoplist$fdrdetails=getlfdr(respoplist$snpmat[,4:5],fdrmax=fdrmax)
    respoplist$threshfdr=respoplist$fdrdetails$thresh
    ## Compute GW thresholds
    respoplist$threshgw=gwsignificance(respoplist$snpmat[,4:5],pgw)
    respoplist$threshgwe=gweffective(respoplist$snpmat[,4:5],respoplist$threshgw,lag.max)
    ## Extract SNP windows
    respoplist$snpwindowsfdr<-getSnpRangesPval(respoplist$snpmat[,1:5],
                                            respoplist$threshfdr,
                                            snpwindow = snpwindow)
    respoplist$snpwindowsgw<-getSnpRangesPval(respoplist$snpmat[,1:5],
                                            respoplist$threshgw,
                                            snpwindow = snpwindow)
    respoplist$snpwindowsgwe<-getSnpRangesPval(respoplist$snpmat[,1:5],
                                            respoplist$threshgwe,
                                            snpwindow = snpwindow)
    respoplist
}

makeMeanPainting=function(importeddataa){
    ## Extract the mean painting over all chromosomes as a single matrix
    sapply(importeddataa,function(x) {  
        do.call("c",lapply(chromosomes,function(chron){
            x$chrdata[[chron]]$snpsums/x$chrdata[[chron]]$nhaps
        }))})
}

whitened_pval=function(mypainting){
    ## Combined p-value anomaly score for all populations
    ## Using Whitening
    mym=colMeans(mypainting)
    mymd=t(t(mypainting)-mym)
    mymdcov=cov(mymd)
    mymdcorr=cor(mymd)
    mymdinvcov=solve(mymdcov)
    mymdinvcovsvd=svd(mymdinvcov)
    myw=mymdinvcovsvd$u %*% diag(sqrt(mymdinvcovsvd$d))
    myz=mymd %*% myw
    myt=rowSums(myz^2)
    myp=pchisq(myt,dim(mypainting)[2],lower.tail=FALSE)
    data.frame(t=myt,p=myp)
}
unwhitened_rmse=function(mypainting,scaled=TRUE){
    ## Combined RMSE anomaly score for all populations
    ## Treating populations as independent (which they aren't)
    mym=colMeans(mypainting)
    mymd=t(t(mypainting)-mym)
    if(scaled){
        mysds=apply(mymd,2,sd)
        myrmse=sqrt(rowSums(t(t(mymd)/mysds)^2))
    }else{
        myrmse=sqrt(rowSums(mymd^2))
    }
    myrmse
}

