test_borealis <- function(){
    extdata <- system.file("extdata", package="borealis")
    gr <- runSingleNewSample(file.path(extdata,'bismark','patient_72',
                                        'patient_72.gz'),
                        NULL,modelFile=file.path(extdata,'CpG_model_chr14.csv'))
    pVal <- gr[start(gr)==24780288]$pVal
    checkTrue(pVal<0.05)
}
