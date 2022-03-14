\name{loadBismarkData}
\alias{loadBismarkData}
\title{Load in bismark results}
\usage{
loadBismarkData(inDir, suffix, chrs)
}
\arguments{

    \item{inDir}{\code{character(1)} Directory path to bismark results.}

    \item{suffix}{\code{character(1)} File suffix for the bismark coverage
    files.}

    \item{chrs}{A character vector listing the chromosomes to be loaded.}

}
\description{
Load in bismark results.
}
\examples{
extdata <- system.file("extdata","bismark", package="borealis")
loadBismarkData(extdata, "_merged.cov.gz.CpG_report.merged_CpG_evidence.cov.gz",
                'chr1')
}
\value{
Returns an object of "SeqCountSet" class.
}
