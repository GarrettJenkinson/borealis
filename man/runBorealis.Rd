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

    \item{suffix}{\code{character(1)} File suffix for the bismark coverage
    files.}

    \item{nThreads}{\code{numeric(1)} Number of compute threads to be used in
    multithreading computations.}

    \item{minDepth}{\code{numeric(1)} The minimum depth of coverage for sample
    to go into modeling.}

    \item{minSamps}{\code{numeric(1)} The minimum number of samples with
    \code{minDepth} coverage required to build a model at a given CpG site.}

    \item{timeout}{\code{numeric(1)} The maximum time in seconds to spend trying
    to build a model at a given CpG site (if it takes longer, we skip the
    site).}

    \item{laplaceSmooth}{\code{logical(1)} Whether or not to do Laplace (i.e.,
    add one) smoothing on the counts.}

    \item{chrs}{A character vector listing the chromosomes to be loaded.}

    \item{outprefix}{\code{character(1)} The sample output file prefix (can
    include a full file path if current working directory is not desired output
    location).}

    \item{modelOutPrefix}{\code{character(1)} The cohort modeling output file
    prefix (can include a full file path if current working directory is not
    desired output location).}
}
\description{
Run the full borealis pipeline.
}
\examples{
extdata <- system.file("extdata","bismark", package="borealis")
results <- runBorealis(extdata,nThreads=1,chrs="chr1")
file.remove(c('CpG_model_chr1.csv', 'CpG_model_rawMethCount.tsv',
    'CpG_model_rawTotalCount.tsv', 'borealis_pat1_chr1_DMLs.tsv',
    'borealis_pat2_chr1_DMLs.tsv', 'borealis_pat3_chr1_DMLs.tsv',
    'borealis_pat4_chr1_DMLs.tsv', 'borealis_pat5_chr1_DMLs.tsv',
    'borealis_pat6_chr1_DMLs.tsv', 'borealis_pat7_chr1_DMLs.tsv',
    'borealis_pat8_chr1_DMLs.tsv'))
}
\value{
Returns a list with raw data and modeling results. The list has elements
\code{BSobj} and \code{result}, which have the dataset and modeling results,
respectively.
}
