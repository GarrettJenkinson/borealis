# Main function that loads in data and produces all results
runBorealis <- function(inDir,
                suffix="_merged.cov.gz.CpG_report.merged_CpG_evidence.cov.gz",
                nThreads=8,minDepth=4,minSamps=5,timeout=10,laplaceSmooth=TRUE,
                chrs=c(paste0("chr",seq_len(22)),"chrX","chrY"),
                outprefix="borealis_",
                modelOutPrefix="CpG_model"){

    # get parallel backend ready
    if (nThreads>1){
        registerDoMC(nThreads)
    }

    # read in raw data and store matrix versions to disk
    BSobj <- loadBismarkData(inDir,suffix,chrs)
    chr <- as.vector(seqnames(BSobj), mode="character")
    pos <- start(BSobj)
    x <- as.array(getCoverage(BSobj, type="M"))
    n <- as.array(getCoverage(BSobj, type="Cov"))
    write.table(cbind(chr=chr,start=pos,end=pos+1,x),
                    file=paste0(modelOutPrefix,"_rawMethCount_",
                                paste(chrs,collapse="_"),".tsv"),
                    row.names=FALSE,quote=FALSE,sep="\t")
    write.table(cbind(chr=chr,start=pos,end=pos+1,n),
                    file=paste0(modelOutPrefix,"_rawTotalCount_",
                                paste(chrs,collapse="_"),".tsv"),
                    row.names=FALSE,quote=FALSE,sep="\t")

    # build the models and save the parameters to disk
    modelDF <- buildModels(BSobj, nThreads, minDepth=minDepth,
                            minSamps=minSamps,timeout=timeout,
                            laplaceSmooth=laplaceSmooth)
    write.csv(modelDF,
                file=paste0(modelOutPrefix,"_",paste(chrs,collapse="_"),".csv"),
                row.names=FALSE,quote=FALSE)

    # compute the final results and save them to disk
    result <- writeResults(BSobj,modelDF,outprefix,chrs)

    # return data for interactive use
    return(BSobj)
}

# Function to load data from bismark outputs
# assumes following pattern for full paths to bismark coverage gz files:
# ${inDir}/${sampleName}/${sampleName}${suffix}
loadBismarkData <- function(inDir,suffix,chrs){
    message("loading in Bismark data.")
    samples <- list.files(inDir)

    ind <- 1 # do this to stop Rcheck complaining
    dataList <- foreach(ind=seq_along(samples)) %dopar% {
        # read in bismark coverage file for this sample
        samp <- samples[ind]
        tmp <- read.table(gzfile(file.path(inDir,samp,paste0(samp,suffix))),
                            header=FALSE)

        # reformat to chr, pos, N, x columns
        tmp$V7 <- tmp$V5+tmp$V6
        tmp <- tmp[,c(1,2,7,5)]
        colnames(tmp) <- c("chr","pos","N","X")

        # filter to relevant chromosomes
        tmp[tmp$chr %in% chrs,]
    }
    BSobj <- DSS::makeBSseqData(dataList,samples)
}

# Function to build the beta-binomial models for each CpG using the data from
# the entire cohort.
buildModels <- function(BSobj,nThreads,minDepth=4,minSamps=5,
                        timeout=10,laplaceSmooth=TRUE) {

    ## grab counts
    n <- as.array(getCoverage(BSobj, type="Cov"))
    keepInd <- rowSums(n>=minDepth)>=minSamps
    x <- as.array(getCoverage(BSobj, type="M"))[keepInd,]
    n <- n[keepInd,]
    chr <- as.vector(BSobj$chr,mode="character")[keepInd]
    pos <- as.array(getCoverage(BSobj, type="Cov"))[keepInd]
    rm(list="BSobj")

    niter <- nrow(x)
    message('Looping over ',niter,
            ' positions in model building. It will take some time.')
    ind <- 1 # do this to stop Rcheck complaining
    df <- foreach(ind=seq_len(nThreads), .combine=rbind) %dopar% {
        # Subset the full data down to just this processor's data
        currInds <- seq(ind,niter,by=nThreads)
        chrLoc <- chr[currInds]
        posLoc <- pos[currInds]
        x <- x[currInds,]
        n <- n[currInds,]

        # loop along and store results in dfInd,
        # which will be merged by foreach across processors
        dfInd <- data.frame(chr=chrLoc,pos=posLoc,mu=NA_real_,theta=NA_real_)
        for(ind2 in seq_along(currInds)){
            pred <- tryCatch({
                # get this cpg's data and fit model
                x1 <- x[ind2,]
                n1 <- n[ind2,]
                pred <- fitGamlss(x1,n1,minDepth,laplaceSmooth,timeout)
            }, error=function(cond) {return(NULL)})
            if(!is.null(pred)){
                dfInd$mu[ind2] <- pred$mu[1]
                dfInd$theta[ind2] <- pred$sigma[1]
            }
        }
        dfInd
    }
    message('modeling complete. Storing results.')

    ## create result and sort by position
    df <- na.omit(df)
    df <- df[with(df, order(chr,pos)), ]
}

# Helper function for model building at each CpG site
fitGamlss <- function(x1,n1,minDepth,laplaceSmooth,timeout){
    # filter samples with low depth
    x1 <- x1[n1 > minDepth]
    n1 <- n1[n1 > minDepth]

    if(laplaceSmooth){
        # Laplace smoothing before beta binomial modeling
        df <- data.frame(x=x1+1,n=n1+2)
    }
    else{
        df <- data.frame(x=x1,n=n1)
    }

    # do gamlss fit, using timeout and silencing messages
    fit <- R.utils::withTimeout( { purrr::quietly(gamlss::gamlss)(
                            cbind(x,n-x)~1,data=df,
                            family = BB(mu.link = "logit",sigma.link = "log")
                        )$result}, timeout=timeout, onTimeout="silent" )
    if (!is.null(fit)){ # fit will be null if it timed out
        pred <- predictAll(fit,data=df[1,,drop=FALSE])
    } else{
        pred <- NULL
    }
    return(pred)
}

# Helper function to write the results for each individual in cohort after model
# has been built.
writeResults <- function(BSobj,modelDF,outprefix,chrs,minObsDepth=10){
    # pull out data elements
    chr <- as.vector(seqnames(BSobj), mode="character")
    pos <- start(BSobj)
    x <- as.array(getCoverage(BSobj, type="M"))
    n <- as.array(getCoverage(BSobj, type="Cov"))
    rm(list="BSobj")

    message("Writing results for each sample to file.")
    result <- foreach(rep=seq_len(ncol(x)), .combine=c,
                        .errorhandling="remove",.multicombine=TRUE) %dopar% {
        df <- data.frame(chr=chr,pos=pos,x=as.integer(x[,rep]),
                            n=as.integer(n[,rep]))

        df <- df[df$n>=minObsDepth,]
        df <- na.omit(plyr::join(df,modelDF,by=c("chr","pos"),type="inner"))

        df <- computePvalsAndEffSize(df)

        # write sorted table
        write.table(df[order(df$pVal),c(1,2,2,3:length(df))],
                        file=paste0(outprefix,colnames(x)[rep],"_",
                        paste(chrs,collapse="_"),"_DMLs.tsv"),
                        quote=FALSE,row.names=FALSE,sep="\t")

        # return to the foreach a dummy value
        retVal <- 0
        retVal
    }
}

# Helper function to produce p-values and effect sizes from a model given the
# count data in one sample.
computePvalsAndEffSize <- function(df){
    # Compute p-values and effect size statistics
    df$leftTailProb <- pmax(0,gamlss.dist::pBB(df$x,mu=df$mu,
                                                sigma=df$theta,bd=df$n))
    df$rightTailProb <- pmax(0,gamlss.dist::pBB((df$x-1)*(df$x>1),mu=df$mu,
                                                sigma=df$theta,
                                                bd=df$n,lower.tail=FALSE))

    df$pVal <- pmin(1.0,2*pmin(df$rightTailProb,df$leftTailProb))
    df$isHypo <- NA #df$rightTailProb > df$leftTailProb
    df[["rightTailProb"]] <- NULL
    df[["leftTailProb"]] <- NULL
    # will not be what everyone wants/needs, particularly if split on chr
    #df$pAdj <- p.adjust(df$pVal,method="BH")
    df$effSize <- (df$x/df$n) - df$mu
    df$isHypo[df$effSize < -0.1] <- TRUE
    df$isHypo[df$effSize > 0.1] <- FALSE
    return(df)
}

# Main function to run a single new sample against existing model from cohort.
runSingleNewSample <- function(inFile,outFile,minObsDepth=10,
                                modelFile="CpG_model.csv"){
    # read in sample data
    df <- read.table(gzfile(inFile),header=FALSE)
    df$V7 <- df$V5+df$V6
    df <- df[,c(1,2,7,5)]
    colnames(df) <- c("chr","pos","n","x")

    # read in models
    modelDF <- read.csv(modelFile)

    # Do modeling
    df <- df[df$n>=minObsDepth,]
    df <- na.omit(plyr::join(df,modelDF,by=c("chr","pos"),type="inner"))

    df <- computePvalsAndEffSize(df)

    # write sorted table
    df <- df[order(df$pVal),c(1,2,2,3:length(df))]
    if(!is.null(outFile)){
        write.table(df,file=outFile,quote=FALSE,row.names=FALSE,sep="\t")
    }
    return(df)
}

plotCpGsite <- function(cpgSites, sampleOfInterest=NA,
                        modelFile="CpG_model.csv",
                        methCountFile="CpG_model_rawMethCount.tsv",
                        totalCountFile="CpG_model_rawTotalCount.tsv"){
    # read in data
    modelDF <- read.csv(modelFile)
    methDF  <- read.table(methCountFile, header=TRUE)
    methDF$chr <- as.character(methDF$chr)
    totDF   <- read.table(totalCountFile, header=TRUE)
    totDF$chr <- as.character(methDF$chr)
    samps <- names(methDF)[c(-1,-2,-3)]

    # make a plot per cpg
    plots <- list()
    for(cpg in cpgSites){
        # filter down to this locus
        my_chr <- strsplit(cpg,":")[[1]][1]
        my_start <- as.numeric(strsplit(cpg,":")[[1]][2])
        x <- methDF %>% dplyr::filter(.data$chr==my_chr,
                                        .data$start==my_start) %>%
                dplyr::select(-c(1,2,3))
        n <- totDF %>%  dplyr::filter(.data$chr==my_chr,
                                        .data$start==my_start) %>%
                dplyr::select(-c(1,2,3))
        model <- modelDF %>% dplyr::filter(.data$chr==my_chr,
                                            .data$pos==my_start) %>%
                dplyr::select(-c(1,2)) %>%
                dplyr::mutate(alpha=.data$mu/(.data$theta+1e-6),
                                beta=(1-.data$mu)/(.data$theta+1e-6))

        # make binom confidence intervals for each samp's raw data
        df <- data.frame(samples=samps,lowerBnd=NA,est=NA,upperBnd=NA,
                            sampleOfInterest=(samps %in% sampleOfInterest))
        for(ind in seq_along(samps)){
            res <- binom.test(as.numeric(x[ind][[1]]),as.numeric(n[ind][[1]]))
            df[ind,"lowerBnd"] <- res$conf.int[1]
            df[ind,"upperBnd"] <- res$conf.int[2]
            df[ind,"est"] <- res$estimate
        }
        # make plot and store in list to be returned
        plots[[cpg]] <- generatePlot(df,model,sampleOfInterest)
    }
    return(plots)
}

generatePlot <- function(df,model,sampleOfInterest){
    # Make plot of raw data confidence intervals and point estimates
    if(!is.na(sampleOfInterest)){
        p1 <- ggplot2::ggplot(df,aes(y=.data$samples,
                                        color=.data$sampleOfInterest)) +
            ggplot2::geom_errorbar(aes(xmin=.data$lowerBnd,
                                        xmax=.data$upperBnd)) +
            ggplot2::geom_point(aes(x=.data$est),size=2) +
            ggplot2::xlim(c(0,1)) +
            ggplot2::xlab("Percent Methylation")+
            ggplot2::theme(legend.position="none")
    } else{
        p1 <- ggplot2::ggplot(df,aes(y=.data$samples)) +
            ggplot2::geom_errorbar(aes(xmin=.data$lowerBnd,
                                        xmax=.data$upperBnd))+
            ggplot2::geom_point(aes(x=.data$est),size=2) +
            ggplot2::xlim(c(0,1))+
            ggplot2::xlab("Percent Methylation")
    }

    # Make plot of the estimated population beta distribution
    betaDF <- data.frame(x=seq(0,1,length=1000), y=dbeta(seq(0,1,length=1000),
                                                        model$alpha,model$beta))
    p2 <- ggplot2::ggplot(betaDF,aes(x=.data$x,y=.data$y)) +
            ggplot2::geom_line(color="blue") +
            ggplot2::expand_limits(y=0) +
            ggplot2::xlab("") + ggplot2::ylab("Model")

    # return a stacked plot
    return(plot_grid(p2,p1,ncol=1,rel_heights = c(1,2),align="v"))
}

