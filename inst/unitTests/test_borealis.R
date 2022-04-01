test_borealis <- function(){
    extdata <- system.file("extdata", package="borealis")
    df <- runSingleNewSample(file.path(extdata,'bismark','patient_72',
                                        'patient_72.gz'),
                        NULL,modelFile=file.path(extdata,'CpG_model_chr14.csv'))
    pVal <- df$pVal[df$start=="24780288"]
    checkTrue(pVal<0.05)
}
