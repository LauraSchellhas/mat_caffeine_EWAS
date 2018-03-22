# Load in all the data and perform probe QC
# You can do this using your usual pipeline, for example, this is how I do it for ALSPAC:

# Set initial parameters for analysis (TP=time point)
Phenofile = "alspac_pheno"
CellData = "gervinandlyle"
WD = "/panfs/panasas01/sscm/gs8094/EWAS/pat_bmi/alspac"

#Set the working directory
setwd(WD)

#Load description of samples
load("/panfs/panasas01/dedicated-mrcieu/studies/latest/alspac/epigenetic/methylation/450k/aries/released/2016-05-03/data/samplesheet/data.Robj")
samplesheet<-subset(samplesheet, time_point=="cord")
qletB<-samplesheet$ALN[which(samplesheet$QLET=="B")] #find alns for multiple pregnancies
samplesheet.cord<-samplesheet[-which(samplesheet$ALN %in% qletB),] #remove multiple pregnancies

#Load phenotype data
complete.pheno<-readRDS(paste0("pheno/",Phenofile,".rds"))
#Add Sample_Name to Pheno
complete.pheno.cord<-merge(complete.pheno,samplesheet.cord[,c("ALN","Sample_Name")],by.x="aln",by.y="ALN")

#Load the methylation data
load("/panfs/panasas01/dedicated-mrcieu/studies/latest/alspac/epigenetic/methylation/450k/aries/released/2016-05-03/data/betas/data.Robj")
meth.matrix.cord <- norm.beta.random[,samplesheet.cord$Sample_Name]
rm(norm.beta.random)

#Load detection P-values (used to filter out probes with a high detection P-value)
load("/panfs/panasas01/dedicated-mrcieu/studies/latest/alspac/epigenetic/methylation/450k/aries/released/2016-05-03/data/detection_p_values/data.Robj")
pvals.cord <- detp[,samplesheet.cord$Sample_Name] 
rm(detp)

#Load annotation data (used to filter out control probes)
require(meffil)
annotation <- meffil.get.features("450k")

#Filter methylation data (remove control probes and probes with high detection P-values)
pvalue.over.0.05 <- pvals.cord > 0.05
count.over.0.05 <- rowSums(sign(pvalue.over.0.05))
probes.to.exclude.p <- rownames(pvals.cord)[which(count.over.0.05 > ncol(pvals.cord)*0.05)]
snps.and.controls <- as.character(annotation$name[-grep("cg|ch", annotation$name)])
annotation1<- annotation[-which(annotation$name %in% c(snps.and.controls,probes.to.exclude.p)),]
meth.matrix.cord <- subset(meth.matrix.cord, row.names(meth.matrix.cord) %in% annotation1$name)

#Load cell.counts
cell.counts.cord<-read.table(paste0("/panfs/panasas01/dedicated-mrcieu/studies/latest/alspac/epigenetic/methylation/450k/aries/released/2016-05-03/data/derived/cellcounts/cord/",CellData,"/data.txt"),header=T)
colnames(cell.counts.cord)[1]<-"Sample_Name"
complete.pheno.cord <- merge(complete.pheno.cord,cell.counts.cord,by="Sample_Name",all=FALSE)

#make sure all variable names are correct for this analysis
colnames(complete.pheno.cord) <- tolower(colnames(complete.pheno.cord))
colnames(complete.pheno.cord)[1]<-"sample.id"
colnames(complete.pheno.cord)[6] <- "mat.ses"

###################################################
#      MATERNAL CAFFEINE ANALYSIS PLAN CODE       # 
#                    GEMMA SHARP                  #
#              PACKAGES AND FUNCTIONS             #
###################################################

setwd("results")

# Load required packages (if these are not already installed, you will have to install them as the first step)
library(sva)
library(tableone)
library(matrixStats)
library(limma)

# Setup the necessary functions
## Function to remove outliers using the IQR*3 (Tukey) method
IQR.removal <- function(meth.matrix){
  rowIQR <- rowIQRs(meth.matrix, na.rm = T)
  row2575 <- rowQuantiles(meth.matrix, probs = c(0.25, 0.75), na.rm = T)
  maskL <- meth.matrix < row2575[,1] - 3 * rowIQR 
  maskU <- meth.matrix > row2575[,2] + 3 * rowIQR 
  meth.matrix[maskL] <- NA
  meth.matrix[maskU] <- NA
  meth.matrix
}
## Function to generate surrogate variables and merge them with the phenotype data (used to adjust for batch)
SVA.generate <- function(meth.matrix, pheno.data, variable.of.interest, model.covariates,n.sv){
  intersecting.samples <- intersect(pheno.data$sample.id,colnames(meth.matrix))
  pheno.data <- na.omit(pheno.data[which(pheno.data$sample.id %in% intersecting.samples),unique(c("sample.id",variable.of.interest,model.covariates))])
  meth.matrix <- meth.matrix[,match(pheno.data$sample.id,colnames(meth.matrix))]
  k = which(is.na(meth.matrix), arr.ind=TRUE)
  meth.matrix[k] = rowMedians(meth.matrix, na.rm=TRUE)[k[,1]]
  mod = model.matrix(reformulate(paste0("pheno.data$",colnames(pheno.data[-1]))))
  mod0 = mod[,-grep(paste0(variable.of.interest,collapse="|"),colnames(mod))]
  sva.ret = sva(meth.matrix, mod=mod, mod0=mod0, n.sv=n.sv)
  SVs = as.data.frame(sva.ret$sv)
  colnames(SVs) <-paste0("sv",1:ncol(SVs))
  cbind(pheno.data,SVs)
}
## Function to run EWAS
ewas.function <-  function(meth.matrix, pheno.data, variable.of.interest){   
  meth.matrix <- meth.matrix[,match(pheno.data$sample.id,colnames(meth.matrix))]
  model.covariates <- colnames(pheno.data)[-which(colnames(pheno.data) %in% c(variable.of.interest,"sample.id"))]
  des = model.matrix(reformulate(paste0("pheno.data$",c(variable.of.interest,model.covariates))))
  fit = lmFit(meth.matrix, des)
  fit.ebayes = eBayes(fit)
  n = rowSums(!is.na(meth.matrix))
  se = (sqrt(fit.ebayes$s2.post) * fit.ebayes$stdev.unscaled[,grep(paste0(variable.of.interest,collapse="|"),colnames(fit.ebayes$stdev.unscaled))])
  res = data.frame(n=n,
                   coef=fit.ebayes$coefficient[,grep(paste0(variable.of.interest,collapse="|"),colnames(fit.ebayes$coefficient))],
                   se=se,
                   p=fit.ebayes$p.value[,grep(paste0(variable.of.interest,collapse="|"),colnames(fit.ebayes$p.value))])
  res
}

###################################################
# PACE PATERNAL BMI PHASE ONE ANALYSIS PLAN CODE  # 
#                    GEMMA SHARP                  #
#                    CORD BLOOD                   #
###################################################

# Set initial parameters
study <- "ALSPAC" #change to your study identifier
timepoint <- "birth" #change depending on the age of the children with methylation samples. Can be "birth", "early_childhood", "late_childhood", "adolescence" or "adult"
cell.names <- tolower(names(cell.counts.cord))[-1]
traits.and.covariates <- c("mat.caff", "mat.bmi","sex","mat.ses","mat.age","mat.smoking","parity","gest.age")
covariates <- c("mat.bmi","mat.ses","mat.age","mat.smoking","parity", cell.names)

# Check phenotype data
pheno <- complete.pheno.cord
for(i in 1:length(c("sample.id",traits.and.covariates,cell.names))) {
  print(ifelse(c("sample.id",traits.and.covariates,cell.names)[i] %in% colnames(pheno)==FALSE,
               paste("CAUTION: the variable called",c("sample.id",traits.and.covariates,cell.names)[i],"is missing from pheno"),
               paste("variable called",c("sample.id",traits.and.covariates,cell.names)[i],"is present in pheno")))
}

# Load methylation data (meth)
meth <- meth.matrix.cord

# IQR*3 method to remove outliers (if this has not already been applied to your data)
log.iqr <- data.frame(cpgs = row.names(meth),NAs.before.IQR3 = rowSums(is.na(meth)))
meth <- IQR.removal(meth)
log.iqr$NAs.after.IQR3 <- rowSums(is.na(meth))
save(log.iqr, file=paste0(study,".matcaff.logIQR.",timepoint,".Rdata"))

# Generate surrogate variables for technical batch and merge with pheno data to create the phenotype dataframes for the mutually adjusted models
pheno.covs <- SVA.generate(meth, pheno, variable.of.interest = "mat.caff", model.covariates = c(covariates,"sex"),n.sv=20)
pheno.gest <- SVA.generate(meth, pheno, variable.of.interest = "mat.caff", model.covariates = c(covariates,"gest.age","sex"),n.sv=20)
pheno.min <- SVA.generate(meth, pheno, variable.of.interest = "mat.caff", model.covariates = cell.names,n.sv=20)

#Create the phenotype dataframes for the sex stratified EWASs
pheno.covs.boys.only <- pheno.covs[which(pheno.covs$sex == 0),]
pheno.covs.girls.only <- pheno.covs[which(pheno.covs$sex == 1),]

# Summarise pheno data and save summaries as .csv files
min.tableone <- as.data.frame(print(CreateTableOne(data=pheno.min[,-1])),stringsAsFactors=FALSE)
covs.tableone <- as.data.frame(print(CreateTableOne(data=pheno.covs[,-1],factorVars=c("mat.smoking","mat.ses","parity","sex"))),stringsAsFactors=FALSE)
gest.tableone <- as.data.frame(print(CreateTableOne(data=pheno.gest[,-1],factorVars=c("mat.smoking","mat.ses","parity","sex"))),stringsAsFactors=FALSE)

write.csv(min.tableone,file=paste0(study,".matcaff.min.summary.",timepoint,".csv"))
write.csv(covs.tableone,file=paste0(study,".matcaff.covs.summary.",timepoint,".csv"))
write.csv(gest.tableone,file=paste0(study,".matcaff.covs.summary.",timepoint,".csv"))

# Test associations between maternal caffeine and cell types
cells <- pheno.min[,which(colnames(pheno.min) %in% cell.names)]
cells.res <- t(apply(cells,2,function(x) summary(lm(x ~ pheno.min$mat.caff))$coef[2,]))
write.csv(cells.res,file=paste0(study,".matcaff.cells.res.summary.",timepoint,".csv"))

# Run each EWAS
ewas.res.min <- ewas.function(meth, pheno.min, variable.of.interest = "mat.caff")
ewas.res.covs <- ewas.function(meth, pheno.covs[,!colnames(pheno.covs) == "sex"], variable.of.interest = "mat.caff")
ewas.res.gest <- ewas.function(meth, pheno.gest[,!colnames(pheno.gest) == "sex"], variable.of.interest = "mat.caff")
ewas.res.covs.boys.only <- ewas.function(meth, pheno.covs.boys.only[,!colnames(pheno.covs.boys.only) == "sex"], variable.of.interest = "mat.caff")
ewas.res.covs.girls.only <- ewas.function(meth, pheno.covs.girls.only[,!colnames(pheno.covs.girls.only) == "sex"], variable.of.interest = "mat.caff")

# Save EWAS results as an Rdata file
save(list=intersect(ls(),
                    c("ewas.res.min",
                      "ewas.res.covs",
                      "ewas.res.gest",
                      "ewas.res.covs.boys.only",
                      "ewas.res.covs.girls.only")),
     file=paste0(study,".matcaff.ewasresults.",timepoint,".Rdata"))
     
