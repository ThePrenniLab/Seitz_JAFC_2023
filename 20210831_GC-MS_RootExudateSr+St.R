
library(xcms)
library(csu.pmf.tools)
library(RAMClustR)

###############################################
## Section A start - XCMS processing
###############################################
## we first choose our project directory, then define our experimental parameters
wd  <-  choose.dir()  # wd <- getwd()
setwd(wd)
dir.create('datasets')
ExpDes  <-  pmfDefineExperiment() 
save(ExpDes, file = "datasets/ExpDes.Rdata")
# load("ExpDes.Rdata")
# load("datasets/ExpDes.Rdata")


# run XCMS - GCMS
#parameters set as run_pmfWorkflow for GC-MS of RootExudates from covercrops
xset  <-  pmfxcms(
  cores=4, 
  minpw=1.5,
  maxpw=15, 
  ExpDes = ExpDes, 
  filetype = "cdf",
  bwpre=3, 
  bwpost=1, 
  minTIC=10, 
  outTIC=TRUE, 
  outPCA=TRUE, 
  snthresh=5, 
  minfrac=0.05,
  seqskip=0, 
  regroup=FALSE
  )
# load("datasets/xcmsObject.Rdata")

## Section A End - XCMS



###############################################
## Section B Start - Build RC obect, remove features
## found at levels in blank samples similar to QC.
###############################################

## build empty ramclustObj
RC <- rc.get.xcms.data(
  xcmsObj = xset,
  ExpDes = ExpDes,
  mzdec = 1
)

save(RC, file = "datasets/RCobject_0_postXCMS.Rdata")

load("datasets/RCobject_0_postXCMS.Rdata")

fix.sample.names<-read.csv(file='fix.sample.names.csv', header=TRUE, sep=",")
RC$phenoData$sample.names<-fix.sample.names$name

#  ExpDes$design
## turn sample names into factor names
#RC <- rc.expand.sample.names(
#  ramclustObj = RC, 
#  delim = "-",
#  factor.names = c(
#    "ID", "Genotype", 
#    "Type", "Extraction",
#    "Extract", "Rep")
#)

## remove 'system peak' features
##DO NOT run unless you have blanks
RC <- rc.feature.filter.blanks(
  ramclustObj = RC,
  qc.tag = "QC",
  blank.tag = "Blank",
  sn = 3
)
##53.8% of features move forward
##Features which failed to demonstrate signal intensity of at least 3 fold greater in QC samples
##than in blanks were removed from the feature dataset. 2165 of 4691 features were removed.

## replace missing values with low level noise (optionally also zero values)

RC <- rc.feature.replace.na(
  ramclustObj = RC,
  replace.int = 0.1,
  replace.noise = 0.1,
  replace.zero = TRUE,
  samp.max.missing = 0.8
)
##replaced 58 of 171768 total features ( 0 % )
## section B end 


###############################################
## Section C: Feature normalization to correct for analytical drift
###############################################

## create 'qc' plots prenormalization for comparison
RC <- rc.qc(
  ramclustObj = RC, 
  outfile.basename = "preNorm",
  qc.tag = 'QC',
  remove.qc = FALSE
)
## then restore 'qc' samples back to dataset for normalization
#RC <- rc.restore.qc.samples(
#  ramclustObj = RC
#)



## normalize feature values to QC samples by run order and batch effect

RC <- rc.feature.normalize.qc(
  ramclustObj = RC, 
  order = order(RC$phenoData$filenames),
  batch = rep(1, nrow(RC$phenoData)),
  qc.tag = "QC",
  rsq.cut = 0.1,
  p.cut = 0.05
)
##Features were normalized by linearly regressing run order versus qc feature intensities to 
##account for instrument signal intensity drift. Only features with a regression pvalue less 
##than 0.05 and an r-squared greater than 0.1 were corrected.  Of 2526 features, 1281 were 
##corrected for run order effects.

## normalize feature values to total ion signal
## can account for both sensitivity drift and sample 
## concentration differences (i.e. urine dilution)
RC <- rc.normalize.tic(
  ramclustObj = RC
)

## rerun rc.qc to visualize QC variance post-normalization. 
RC <- rc.qc(
  ramclustObj = RC, 
  outfile.basename = "post_QCNorm"
)

save(RC, file = "datasets/RCobject_1_preClustering.Rdata")

# load("datasets/RCobject_1_preClustering.Rdata")

## EVALUATE OUTPUT BEFORE PROCEEDING
## IF DATA QUALITY LOOKS INSUFFICIENT
## FIX BEFORE CLUSTERING AND ANNOTATING

## Section C: End


###############################################
## Section D: Cluster Features
###############################################

# cluster Features into compounds. 
#set certain parameter base don file like xcms above
RC <- rc.ramclustr(
  ramclustObj = RC, 
  st = 10,
  sr = 100,
  maxt = NULL,
  deepSplit = FALSE,
  blocksize = 2000,
  mult = 5,
  hmax = 0.9,
  collapse = TRUE,
  minModuleSize = 3,
  linkage = "average",
  cor.method = "pearson",
  rt.only.low.n = TRUE,
  fftempdir = NULL
)
#RAMClust has condensed 2526 features into 80 spectra collapsing feature into spectral signal intensities 

save(RC, file = "datasets/RCobject_2_postClustering.Rdata")

exportDataset(ramclustObj=RC, which.data="SpecAbund")

RC <- rc.qc(
  ramclustObj = RC, 
  outfile.basename = "post_Clustering",
  qc.tag = 'QC',
  remove.qc = TRUE
)


# load("datasets/RCobject_2_postClustering.Rdata")

# Export spectra
rc.export.msp.rc(ramclustObj = RC, one.file = TRUE)

exportDataset(ramclustObj=RC, which.data="SpecAbund")

save(RC, file = "datasets/RCobject_2_postClusteringandQCremoval.Rdata")

## Section D: End


#####################

# BREAK in SCRIPT:
# ____ run RAMSearch manually

## import MSFinder results
RC <- impRamSearch(ramclustObj = RC)

RC <- annotate(
  ramclustObj = RC,
  standardize.names = FALSE,
  min.msms.score = 6.5,
  database.priority = NULL,
  citation.score = TRUE,
  find.inchikey = FALSE,
  taxonomy.inchi = NULL,
  reset = TRUE)


## Error with this command <- Error in `[.data.frame`(tmp, , 1:5) : undefined columns selected
## use pubchem rest API to retrieve a great deal of data on each annoation
# RC <- rc.cmpd.get.pubchem(ramclustObj = RC)

# The below didn't give an error but likely didn't work based upon the fact that the pubchem command didn't work
## use classyFire to retreive compounds classification for each annotation
# RC <- getClassyFire(ramclustObj = RC)

RC <- degolm(ramclustObj = RC)

save(RC, file = "datasets/RCobject_3_postAnnotation.Rdata")

## Section E: End



###############################################
## Section F: Statistical Analysis
###############################################

# Examine factors to guide analysis
des <- getData(ramclustObj = RC)[[1]]
head(des, n = 3)



#  PCA
# factors were misnamed: sn=sn, species = label, tissue = species, date = tissue, label=date, trt = trt, use = use
RC <- pmfpca(ramclustObj = RC,  
             which.factors = c("date","trt","label"), 
             npc = "auto", num.factors = NULL)


#  ANOVA
RC <- pmfanova(ramclustObj = RC, anova.call = "label", subset=c("date", "Leaf", "use", "Change over time"))
RC <- pmfanova(ramclustObj = RC, anova.call = "trt*date", subset=c( "use", "Metabolomics at Harvest", "label", "March8"))
RC <- pmfanova(ramclustObj = RC, anova.call = "trt*date")

save(RC, file = "datasets/RCobject_4_postStats.Rdata")


## Section F: End



###############################################
## Section G: Reporting and file export
###############################################

## Export SpecAbund Dataset - contains quantitative signal intensity values
exportDataset(ramclustObj=RC, which.data="SpecAbund")

## Export Annotation Summary file (spectra directory)
annotation.summary(ramclustObj=RC)

## Export methods summary for all R based post-XCMS steps
write.methods(ramclustObj = RC)

## print sample names to console, return summarized values in RC object
RC <- reportSampleNames(ramclustObj = RC)

## zip files for sharing
make.zip.files(do.raw = FALSE)

## write a text file with descriptions of output files.  
write.file.summary()

## Section G: End


