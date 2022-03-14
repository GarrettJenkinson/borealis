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
    write.table(cbind(as.vector(BSobj$chr,mode="character"),
                        BSobj$pos,BSobj$pos,BSobj$x),
                    file=paste0(modelOutPrefix,"_rawMethCount.tsv"),
                    row.names=FALSE,quote=FALSE,sep="\t")
    write.table(cbind(as.vector(BSobj$chr,mode="character"),
                        BSobj$pos,BSobj$pos,BSobj$n),
                    file=paste0(modelOutPrefix,"_rawTotalCount.tsv"),
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

    # return data and results for interactive use
    return(list(BSobj=BSobj,result=result))
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
    BSobj <- list(  chr=seqnames(BSobj),
                    pos=start(BSobj),
                    x=as.array(getCoverage(BSobj, type="M")),
                    n=as.array(getCoverage(BSobj, type="Cov")))
}

# Function to build the beta-binomial models for each CpG using the data from
# the entire cohort.
buildModels <- function(BSobj,nThreads,minDepth=4,minSamps=5,
                        timeout=10,laplaceSmooth=TRUE) {
    ## grab counts
    keepInd <- rowSums(BSobj$n>=minDepth)>=minSamps
    x <- BSobj$x[keepInd,]
    n <- BSobj$n[keepInd,]
    chr <- as.vector(BSobj$chr,mode="character")[keepInd]
    pos <- BSobj$pos[keepInd]
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

    return(df)
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

    nreps1 <- ncol(BSobj$x)
    repNames <- colnames(BSobj$x)

    message("Writing results for each sample to file.")
    result <- foreach(rep=seq_len(nreps1), .combine=c,
                        .errorhandling = "remove",.multicombine=TRUE) %dopar% {
        df <- data.frame(chr=as.vector(BSobj$chr,mode="character"),
                            pos=BSobj$pos,
                            x = as.integer(BSobj$x[,rep]),
                            n=as.integer(BSobj$n[,rep]))

        df <- df[df$n>=minObsDepth,]
        df <- na.omit(plyr::join(df,modelDF,by=c("chr","pos"),type="inner"))

        df <- computePvalsAndEffSize(df)

        # write sorted table
        write.table(df[order(df$pVal),c(1,2,2,3:length(df))],
                        file=paste0(outprefix,repNames[rep],"_",
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
