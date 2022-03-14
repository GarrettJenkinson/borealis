\name{runSingleNewSample}
\alias{runSingleNewSample}
\title{Run a single new sample after modeling complete}
\usage{
runSingleNewSample(inFile, outFile, minObsDepth=10, modelFile="CpG_model.csv")
}
\arguments{

    \item{inFile}{\code{character(1)} File name (including full path if not
    current working directory) to the bismark coverage file.}

    \item{outFile}{\code{character(1)} File name (including full path if not
    current working directory) for the sample's modeling outputs. If \code{NULL}
    is provided, no outputs will be written to disk.}

    \item{minObsDepth}{\code{numeric(1)} Minimum depth of coverage in this
    sample for a modeling output/p-value to be produced at a given CpG.}

    \item{modelFile}{\code{character(1)} File name (including full path if not
    current working directory) for the model files (built by running
    \code{runBorealis} function).}

}
\description{
Run a single new sample after modeling complete.
}
\examples{
extdata <- system.file("extdata", package="borealis")
runSingleNewSample(file.path(extdata,'bismark','pat1',
                'pat1_merged.cov.gz.CpG_report.merged_CpG_evidence.cov.gz'),
                'out.txt',modelFile=file.path(extdata,'CpG_model_chr1.csv'))
file.remove('out.txt')
}
\value{
Returns a \code{data.frame} with modeling results.
}
