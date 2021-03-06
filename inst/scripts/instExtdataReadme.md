The data utilized in the vignette and included in the "inst/extdata/bismark"
directory are outputs from the Bismark Bisulfite Mapper
for Bisulfite next-generation sequencing data available at:

https://www.bioinformatics.babraham.ac.uk/projects/bismark/

Paired end sequencing data (FASTQ format) were used as input to Bismark and
aligned to human genome version hg19. Comprehensive instructions on how to run
Bismark are provided in its user guide at:

https://rawgit.com/FelixKrueger/Bismark/master/Docs/Bismark_User_Guide.html

Assuming we have paired-end sequencing data for Patient 7 in FASTQ format, and the
forward and reverse read files are named Patient_7_1.fq and Patient_7_2.fq,
then a command similar to the following could be run to generate Bismark
alignments:

1) Call bismark 

    bismark \
    --bowtie2 \
    --path_to_bowtie2 /your/pathTo/bowtie2/ \
    --output_dir Patient_7 \
    --genome /path/to/genome/hg19 \
    -1 Patient_7_1.fq -2 Patient_7_2.fq


Then post-alignment extract the required outputs for downstream work:

2) Call methylation extractor:
 
    bismark_methylation_extractor \
    --gzip \
    --no_overlap \
    --bedGraph \
    --counts \
    --ignore 4 \
    --ignore_3prime_r2 4 \
    --buffer_size 50G \
    --comprehensive \
    --merge_non_CpG \
    --output Patient_7 \
    --genome_folder /path/to/genome/hg19 \
    Patient_7.bam

3) Call coverage to cytosine conversion:

    coverage2cytosine \
    --merge_CpG \
    --gzip \
    --genome_folder /path/to/genome/hg19 \
    --output Patient_7_merged.cov.gz \
    --dir Patient_7 \
    Patient_7/Patient_7_bismark_bt2_pe.bismark.cov.gz

Running these series of commands for multiple samples would produce outputs
similar to those in the BOREALIS package inst/extdata/bismark folder. Note the
outputs were subsetted to a small region of chromosome 14 for the purposes of
creating a lightweight vignette.

The Bismark output used as direct input to BOREALIS is generated by default with
suffix "_merged.cov.gz.CpG_report.merged_CpG_evidence.cov.gz". This is the
default suffix specified in the runBorealis() method which is used at the outset
of the vignette.

Following the examples above, you will generate outputs per patient, with each
patient's Bismark data in a separate, named directory. The parent directory
containing these patient-specific drectories is supplied to runBorealis() as
explained in the vignette.

The data included in the "inst/extdata" directory (prefixed "Cpg_model_*") are
model building files output by runBorealis() and are provided to enable running
the example code provide with the BOREALIS package, e.g.

    runBorealis(extdata,nThreads=2, chrs="chr14", suffix=".gz",
    outprefix=file.path(extdata,"vignette_borealis_"),
    modelOutPrefix=file.path("my/output_dir","vignette_CpG_model"))


Thie information here, combined with the vignette, should provide you with all
you need to fully replicate your own analysis. Good luck!
