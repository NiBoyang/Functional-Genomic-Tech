# At the top of each section, it was stated that which tutorial the codes were
# mainly from, specific modification would be comment on that line.
# invisible() funtion was used to suppress the system output, make the output more tidy

#===============================================================================
#import required libraries======================================================
#===============================================================================
library(limma)
library(affy)
library(gcrma)
library(annotate)
library(mouse4302.db)# load chip-specific annotation
library(scatterplot3d)
#BiocManager::install("arrayQualityMetrics")
library(arrayQualityMetrics)
options(warn=-1)# suppress warning

#===============================================================================
#QC of Raw Data=================================================================
#===============================================================================
##load the files for the workflow, this part of code is mainly from tutorial 2
# load the targets file and get affyBatch from it
targets <- read.table("targets.txt",header=T,as.is=T)
rawData <- ReadAffy(filenames=targets$Filename, sampleNames=targets$Name)

# create a directory for QC results
invisible(dir.create("QC_Results"))

# boxplot of raw data
png(filename = "./QC_Results/boxplot_raw.png")
par(mar=c(3,7,2,1))
# add xlabel, title, and change to horizontal
boxplot(rawData, cex.axis = 1, las = 1,
        xlab = "Expression",
        col = targets$Colour, 
        main = "Distribution of raw intensity values",
        horizontal = T)
invisible(dev.off())
print("boxplot of raw data created -> ./QC_Results/boxplot_raw.png")

print("Running comprehensive QC analysis for the raw data!")

# modified from https://www.bioconductor.org/packages/release/bioc/vignettes/arrayQualityMetrics/inst/doc/arrayQualityMetrics.R
# comprehensive QC of the raw data
arrayQualityMetrics(expressionset = rawData,
                    outdir = "./QC_Results/Report_for_raw",
                    force = TRUE,
                    do.logtransform = TRUE)
print("Completed")
print("Result can be found in './QC_Results/Report_for_raw/index.html'")

#===============================================================================
#QC of Normalized Data==========================================================
#===============================================================================
# using gcRMA to normalize the raw data
print("Starting normalization by gcRMA")

# set fast to False, referred to https://www.bioconductor.org/packages/devel/bioc/manuals/gcrma/man/gcrma.pdf
eset <- gcrma(rawData, fast=F)

# boxplot of normalized data
png("./QC_Results/boxplot_normalized.png")
par(mar=c(3,7,2,1))
boxplot(eset, cex.axis = 1, las = 1,
        xlab = "Expression",
        col = targets$Colour, 
        main = "Distribution of normalized intensity values",
        horizontal = T)
invisible(dev.off())
print("boxplot of raw data created -> ./QC_Results/boxplot_normalized.png")

# comprehensive QC for the normalized data
print("Running comprehensive QC analysis for the normalized data!")
arrayQualityMetrics(expressionset = eset,
                    outdir = "./QC_Results/Report_for_normalized",
                    force = TRUE)
print("Completed")
print("Result can be found in './QC_Results/Report_for_normalized/index.html'")

# get expression value matrix
exprsvals <- exprs(eset)

# 3d-PCA plot of the normalized data, from tutorial 3
pca <- prcomp(t(exprsvals), scale = T) 
png("./QC_Results/3d-pca.png")
s3d <- scatterplot3d(pca$x[,1:3], pch = 19, color = targets$Colour, 
                     main = "3D-PCA of Normalized Dataset")
s3d.coords <- s3d$xyz.convert(pca$x[,1:3])
text(s3d.coords$x, s3d.coords$y, labels = colnames(exprsvals),pos = 3,offset = 0.5)
invisible(dev.off())
print("The 3D-PCA plot of normalized data created -> ./QC_Results/3d-pca.png")

# mva plot, codes from tutorial 3
png("./QC_Results/MVA_plot.png")
mva.pairs(exprsvals, cex = 1, main = "MVA plot of normalized data")
invisible(dev.off())
#===============================================================================
#Fold Change Analysis===========================================================
#===============================================================================
## perform fold filtering, this part of code is mainly modified from tutorial 3
# gcRMA outputs log2 data while MAS5 outputs linear data
# To convert from log¦
exprsvals10 <- 2^exprsvals

# get the mean of infected samples
c2c12.inf <- apply(exprsvals10[,c("C2C12_Inf.1",
                                  "C2C12_Inf.2",
                                  "C2C12_Inf.3")],1,mean)

# get the mean of control samples
c2c12.ctrl <- apply(exprsvals10[,c("C2C12_Ctrl.1",
                                   "C2C12_Ctrl.2",
                                   "C2C12_Ctrl.3")],1,mean)

# calculate the fold
c2c12.inf_to_ctrl <- c2c12.inf / c2c12.ctrl

# combine all the data including log2-transformation
all.data = cbind(c2c12.ctrl, c2c12.inf, c2c12.inf_to_ctrl, abs(1 - c2c12.inf_to_ctrl))
colnames(all.data)[4] <- "abs(1-FC)"
# write the combined data to file
write.table(all.data, file="Fold_change.txt", quote=F, sep = '\t', col.names = NA)
print("Fold Change Result created -> ./Fold_change.txt")

#===============================================================================
#Statistical Analysis===========================================================
#===============================================================================
## Beginning statistical analysis, code mainly modified from tutorial 4
#build an annotation table
ID <- featureNames(eset)
Symbol <- getSYMBOL(ID, "mouse4302.db")
Name <- as.character(lookUp(ID, "mouse4302.db", "GENENAME"))
tmp <- data.frame(ID=ID, Symbol=Symbol, Name=Name, stringsAsFactors=F)
tmp[tmp=="NA"] <- NA #fix padding with NA characters 

#assign as feature data of the current Eset
fData(eset) <- tmp

## Statistical analysis using Limma, code mainly modified from tutorial 4
# create design matrix
design <- model.matrix(~-1+factor(c(1,1,1,2,2,2)))
colnames(design) <- c("C2C12_Inf", "C2C12_Ctrl")

# create constratmatrix, where Limma comparisons to make
contrastmatrix <- makeContrasts(C2C12_Inf-C2C12_Ctrl,
                                levels=design)

# issue these commands to fit the model
# and make the contrasts
fit <- lmFit(eset, design)
fit2 <- contrasts.fit(fit, contrastmatrix)

# this last part essentially moderates the t-statistic using 
# the borrowed variance approach described in class
fit2 <- eBayes(fit2)

# get the results
myresults <-topTable(fit2,coef=1, adjust="fdr", number=nrow(eset))
write.table(myresults,"differential.txt")
print("Statistical Differential Analysis Result created -> ./differential.txt")

#===============================================================================
#Enrichment Analysis============================================================
#===============================================================================
## Carry out Functional Enrichment analysis, code mainly modified from tutorial 5
# read annotation database
Mm.H <- readRDS("Mm.h.all.v7.1.entrez.rds") 

## Process annotation for functional enrichment

# Here we select from the annotation a number of keys with the primary key being PROBEID
# due to the use of arrayQualityMetrics, use AnnotationDbi::select instead
# ref: https://support.bioconductor.org/p/101902/#101903
res <- AnnotationDbi::select(mouse4302.db, keys = rownames(eset), 
                             columns = c("ENTREZID", "ENSEMBL","SYMBOL"), 
                             keytype="PROBEID")

#find the index of each row of the expression set in the #annotation object res
idx <- match(rownames(eset), res$PROBEID)

#Use the index to set the phenotypic data in the ExpressionSet
fData(eset) <- res[idx, ]

#Find all rows that don't have an EntrezID and remove then
eset_t<-eset[is.na(fData(eset)$ENTREZID)==0,]


## Functional Enrichment Analysis, codes referred to Tutorial 5

#convert to indexes
H.indices <- ids2indices(Mm.H,fData(eset_t)$ENTREZID)

# run mroast for enrichmen analysis
results <-camera(eset_t,index=H.indices,design=design,contrast=contrastmatrix[,1],adjust.method = "BH")

# write the enrichment analysis results to file
write.table(results,"enrichment.txt",sep="\t")

print("Functional Enrichment Analysis result created -> ./enrichment.txt")
