test_borealis <- function(){
    extdata <- system.file("extdata", package="borealis")
    df <- runSingleNewSample(paste0(extdata,'/bismark/pat1/',
        'pat1_merged.cov.gz.CpG_report.merged_CpG_evidence.cov.gz'),
        NULL,modelFile=paste0(extdata,'/CpG_model_chr1.csv'))
    checkTrue(df$pVal[1]<0.05)
}
