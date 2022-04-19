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

    \item{minObsDepth}{(optional) \code{numeric(1)} Minimum depth of coverage in
    this sample for a modeling output/p-value to be produced at a given CpG.}

    \item{modelFile}{(optional) \code{character(1)} File name (including full
    path if not current working directory) for the model files (built by running
    \code{runBorealis} function).}

}
\description{
Run a single new sample after modeling using \code{runBorealis} has already
been completed in a cohort of samples. It will not rebuild the models and only
predict using previously estimated model specified by \code{modelFile}.
}
\examples{
extdata <- system.file("extdata", package="borealis")
outdir <- tempdir()
df <- runSingleNewSample(file.path(extdata,'bismark','patient_72',
                            'patient_72.gz'),file.path(outdir,'output.txt'),
                            modelFile=file.path(extdata,'CpG_model_chr14.csv'))
}
\value{
Returns a \code{GRanges} object with modeling results.
}
