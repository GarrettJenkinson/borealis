\name{runBorealis}
\alias{runBorealis}
\title{Run the full borealis pipeline}
\usage{
runBorealis(inDir,
            suffix ="_merged.cov.gz.CpG_report.merged_CpG_evidence.cov.gz",
            nThreads = 8, minDepth = 4, minSamps = 5, timeout = 10,
            laplaceSmooth = TRUE,
            chrs = c(paste0("chr",seq_len(22)), "chrX", "chrY"),
            outprefix = "borealis_", modelOutPrefix = "CpG_model")
}
\arguments{

    \item{inDir}{\code{character(1)} Directory path to bismark results. NOTE:
    this assumes following pattern for full paths to bismark coverage gz files:
    \code{${inDir}/${sampleName}/${sampleName}${suffix}}}

    \item{suffix}{(optional) \code{character(1)} File suffix for the bismark
    coverage files.}

    \item{nThreads}{(optional) \code{numeric(1)} Number of compute threads to be
    used in multithreading computations.}

    \item{minDepth}{(optional) \code{numeric(1)} The minimum depth of coverage
    for sample to go into modeling.}

    \item{minSamps}{(optional) \code{numeric(1)} The minimum number of samples
    with \code{minDepth} coverage required to build a model at a given CpG
    site.}

    \item{timeout}{(optional) \code{numeric(1)} The maximum time in seconds to
    spend trying to build a model at a given CpG site (if it takes longer, we
    skip the site).}

    \item{laplaceSmooth}{(optional) \code{logical(1)} Whether or not to do
    Laplace (i.e., add one) smoothing on the counts.}

    \item{chrs}{(optional) A character vector listing the chromosomes to be
    loaded.}

    \item{outprefix}{(optional) \code{character(1)} The sample output file
    prefix (can include a full file path if current working directory is not
    desired output location).}

    \item{modelOutPrefix}{(optional) \code{character(1)} The cohort modeling
    output file prefix (can include a full file path if current working
    directory is not desired output location).}
}
\description{
Run the full borealis pipeline. It will load in bismark data and save out to
disk matrix-based methylation and total count files, then it will build the
beta-binomial statistical models for the cohort at each CpG site and save the
parameters of this model to disk, and finally provide outlier p-values and
summary statistics for each sample in the cohort at each CpG site.
}
\examples{
extdata <- system.file("extdata","bismark", package="borealis")
outdir <- tempdir()
results <- runBorealis(extdata,nThreads=2,chrs="chr14",suffix=".gz",
                        outprefix = file.path(outdir,"borealis_"),
                        modelOutPrefix = file.path(outdir,"CpG_model"))
}
\value{
Returns an object of "BSseq" class with raw dataset loaded and used for
modeling purposes.
}
