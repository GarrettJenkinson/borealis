\name{plotCpGsite}
\alias{plotCpGsite}
\title{Generate a plot of the model and raw data at one or more CpG sites}
\usage{
plotCpGsite(cpgSites, sampleOfInterest=NA, modelFile="CpG_model.csv",
                methCountFile="CpG_model_rawMethCount.tsv",
                totalCountFile="CpG_model_rawTotalCount.tsv")
}
\arguments{

    \item{cpgSites}{A character vector of CpG sites specified as "chr1:71732"
    representing the chromosome and start position of the CpG site. A separate
    plot will be generated for each site specified.}

    \item{sampleOfInterest}{(optional) \code{character(1)} Name of sample of
    interest which will be colored differently than the rest of the samples in
    the cohort. If \code{NA} then all samples will be plotted with same color.}

    \item{modelFile}{\code{character(1)} The mode file (including full path if
    not current working directory) with beta-binomial parameter estimates
    produced by \code{runBorealis}.}

    \item{methCountFile}{\code{character(1)} File name (including full path if
    not current working directory) for the methylated count file produced by
    \code{runBorealis}.}

    \item{totalCountFile}{\code{character(1)} File name (including full path if
    not current working directory) for the total count file produced by
    \code{runBorealis}.}

}
\description{
Generate plots of model and results. The top panel of the plot will be the beta
distribution in the beta-binomial model estimated for the cohort. The bottom
panel will be the 95 percent confidence intervals around the percent methylation
in each sample at that CpG site.
}
\examples{
extdata <- system.file("extdata", package="borealis")
plots <- plotCpGsite(c("chr1:71732","chr1:71788"),
                sampleOfInterest="pat1",
                modelFile=file.path(extdata,'CpG_model_chr1.csv'),
                methCountFile=file.path(extdata,'CpG_model_rawMethCount.tsv'),
                totalCountFile=file.path(extdata,"CpG_model_rawTotalCount.tsv"))

}
\value{
Returns a \code{list} with each element indexed by the provided \code{cpgSites}
and storing a ggplot/cowplot object.
}
