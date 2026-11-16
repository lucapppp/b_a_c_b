###################################################
### chunk number 1: 
###################################################
setwd("~/BCB ESERCITAZIONI/Esercitazione_5 - E2B - 17nov2026 - ESERCITAZIONE R LIMMA/WD_R")
# MCF-7 human breast cancer cells treated Dioxin, (GSE7765)
# Let’s analyse gene expression changes without and with the treatment 
x <-read.table("./Data/Example-matrix.txt", skip=66, header=TRUE, sep="\t", row.names=1)
# HU133A and B single-channel arrays were used to quantify differences in RNA after treatment with 100 nM dioxin for 16 h. 


###################################################
### chunk number 2: 
###################################################
dim(x)      # Dimensions
names(x)    # Vector of strings to name the columns
attach(x)   # This makes x the default data frame so that we don't have to use 
# x$ to access columns, E.G.: 
x$GSM188013[1:5]
GSM188013[1:5]


###################################################
### chunk number 3: 
###################################################
summary(x)
round(apply(x,2, summary),3)  # Column-wise summary statistics of x


###################################################
### chunk number 4: 
###################################################
logX <- log2(x)
round(apply(logX,2, summary),3)  # Column-wise summary statistics of log2x


###################################################
### chunk number 5: hist.array.0
###################################################
# histograms comparing the expression values of an array (e.g. #1) with its logged version

opt <- par(mfrow=c(1,2)) #one row and two columns
hist(x[,1])  
hist(logX[,1])


###################################################
### chunk number 6: hist.array.1
###################################################
opt <- par(mfrow=c(1,1))
hist(logX[,1], main="Expression values",
    xlab="Expression (Log2 transformed)",
    ylab="# of spots",
    col=4)  


###################################################
### chunk number 7: MAPlot.1
###################################################
# scatter plot and MA plot
opt<-par(mfcol=c(2,1),mar = c(5.1, 4.1, 4.1, 2.1)) # par(mar = c(bottom, left, top, right)
plot(X1<-x[,1], Y1<-x[,2]) #compare pre-treatment sample with after-treatment sample   
abline(0,1, col="yellow")
A.1<-0.5*log2(X1*Y1) #average intensity
M.1 <- log2(X1/Y1) #log2ratio
plot(A.1, M.1)
abline (h=0,col="yellow")



###################################################
### chunk number 8: boxplot.1
###################################################
opt <- par(mfcol=c(1,2)) 
boxplot(x, col=c(2,3,2,3,2,3),main="EV for 3 control and 3 treated samples",
    ylab="Expression", las=2, cex.axis=0.7, cex.main=0.9) # las=2 all axis labels vertical (1 for horizontal)
boxplot(logX, col=c(2,3,2,3,2,3), main="EV for 3 control and 3 treated samples",
    ylab="Expression (Log2)", las=2, cex.axis=0.7, cex.main=0.9)
abline(0,0, col="black")   



###################################################
### chunk number 9: toPdf
###################################################
pdf("diagnostics.pdf")
opt <- par(mfcol=c(1,2)) 
boxplot(x, col=c(2,3,2,3,2,3),main="EV for 3 control and 3 treated patients",
    ylab="Expression", las=2, cex.axis=0.5, cex.main=0.9)
boxplot(logX, col=c(2,3,2,3,2,3), main="EV for 3 control and 3 treated patients",
    ylab="Expression", las=2, cex.axis=0.5, cex.main=0.9)
abline(0,0, col="black")   
dev.off()


###################################################
### chunk number 10: labels
###################################################
# Looking for differentially expressed genes:
# first, we label treated samples and control samples. 
descriptions<-read.table("./Data/Example-matrix.txt", skip=29, header=F, sep="\t", row.names=1, nrows=1)
descriptions
classes<-c(0,1,0,1,0,1)


###################################################
### chunk number 11: compute expression ratios
###################################################
# Then, we compute the treated/control expression ratios as the difference between 
# the treated and the control log-transformed mean expression values
get.er<-function(data){
  return(mean(data[classes==1])-mean(data[classes==0]))}

ers=apply(logX,1,get.er) # apply over the rows, i.e. gene probes
#View(ers)
plot(ers ,ylab = 'Mean Log2ratio', xlab = "Gene probes")
# over/under expressed genes (at least of a factor 4)
plot(which(abs(ers)>2),ers[abs(ers)>2],ylab = 'Mean Log2ratio', xlab = "Gene probes") # a factor 4 of difference (2 in log2ratio)
length(ers[abs(ers)>2])
names(ers[abs(ers)>2])

#BiocManager::install("hgu133plus2.db")
library(hgu133plus2.db)
probes <- names(ers[abs(ers)>2])
mapped_genes <- unique(mapIds(hgu133plus2.db, keys=probes, column="SYMBOL", keytype="PROBEID"))


###################################################
### chunk number 12: Student's t-test
###################################################
      ttest=function(x){
       tt=t.test(x[classes==0],x[classes==1])
       return(c(tt$statistic,tt$p.value))}
# it computes significance of the difference between the means of the two groups

###################################################
### chunk number 13: Student's t-test
###################################################
# compare these data using the Student's t-test
ans=apply(logX,1,ttest)
row.names(ans)<-c('ttest', 'pval')
ts<- ans[1,] # tscores: differences between each pair of means
pvals<-ans[2,]
res<-t(ans)
#View(res)

View(final<-res[which(pvals<0.05),])
#View(final<-res[which(res[,'pval']<0.05),])
dim(final)[1]
rownames(final)
mapped_genes_tt <- unique(mapIds(hgu133plus2.db, keys=rownames(final), column="SYMBOL", keytype="PROBEID"))


###################################################
### chunk number 14: hist
###################################################
opt <- par(mfcol=c(1,1)) 
hist(ts)


###################################################
### chunk number 15: QQ plot
###################################################
# qqnorm is a generic function which produces a normal QQ plot of the values in y.
# qqline adds a line which passes through the probe lower and upper quartiles
opt <- par(mfcol=c(1,1)) 
qqnorm(ts)
qqline(ts)
abline(0,1, col="yellow")


###################################################
### chunk number 16: howManyGenes
###################################################
for (i in c(0.05, 0.01, 0.001, 0.0001, 0.00001)){
  print(paste("Probes with p-values smaller than", i, length(which(pvals < i))))
}

###################################################
### chunk number 17: whichGenes
###################################################
which(pvals<0.001) #positions and names for a named vector
View(pvals[pvals<0.001]) #extract cases of interest

#-------------------------------------------------------------------------------------------------------------

# EX.1
# Two-channel microarrays preprocessing

# For this tutorial we will consider the ApoA1 case-study
# .	ApoA1 (apolipoprotein A1) plays an important role in HDL (High Density Lipoprotein) metabolism (“good” cholesterol)
# .	In this experiment we want to compare 8 ApoA1 knock-out mice with 8 normal C57BL/6 ("Black six") mice.
# .	Target/test mRNA has been labelled with Cy5 (red-fluorescent) dye
#.	The common control -wild type reference- has been obtained by pooling RNA extracted from the 8 control mice
# .	Common control reference mRNA was labelled with Cy3 (yellow-fluorescent) dye
 

rm(list = ls())
# Load limma library
library(limma)
#limmaUsersGuide()
# Extract the ApoA1 folder and change your current working directory with the one containing the data 
setwd("Data/ApoAI")

# Load the metadata target file (ApoAITargets.txt)
targets <- readTargets("ApoAITargets.txt")

# Explore the targets data
targets

# Load microarray data
RG <- read.maimages(targets, source="spot") #source is different for different MA types

# The RG object contains R/G channel intensity data and metadata.
# The source option is needed to specify the source file type. 
# Since we have spot files, this option must be set to "spot".
# Unlike most microarray image analysis software, SPOT does not 
# store probe annotation together with the output (.spot file), 
# so a GAL (GenePix Array List) file is used to store gene name information. 
# We can add this information to the RG object using:
RG$genes <- readGAL("ApoAI.gal")

# Background correction
RG_b <- backgroundCorrect(RG, method="normexp")

# Within array normalization
opt <- par(mfrow=c(2,1))
MA_w <- normalizeWithinArrays(RG_b, method="none")
plot(MA_w[['A']], MA_w[['M']])
abline (h=0,col="yellow")

MA_w <- normalizeWithinArrays(RG_b, method="loess") #locally weighted regression model normalization
plot(MA_w[['A']], MA_w[['M']])
abline (h=0,col="yellow") #please zoom in the plot

# Plot the intensity density 
x11(); plotDensities(MA_w) #x11 is used to open an interactive graphics device

#Between array normalization/summarization
MA_q <- normalizeBetweenArrays(MA_w, method="quantile") # default for single-channel MAs

MA_aq <- normalizeBetweenArrays(MA_w, method="Aquantile") # default for two-channel MAs
# Quantile normalization was explored by Yang and Thorne (2003) for two-color cDNA arrays. 
# method="quantile" ensures that the intensities have the same empirical distribution 
# across arrays and across channels. 
# method="Aquantile" ensures that the A-values (average intensities) have the same empirical 
# distribution across arrays leaving the M-values (log-ratios) unchanged. 
# Plot the intensity density histogram
x11(); plotDensities(MA_q)

x11(); plotDensities(MA_aq)

# Generate an intensity boxplot for the tested cases (Cy5 dye)
# Set color vectors for each condition (blue
# for C57BL/6 and yellow/gold for ApoAI-/-)
colors <- MA_q$targets$Cy5
colors
colors[colors=="wild type"] <- "skyblue"
colors[colors=="ApoAI KO"] <- "gold"
colors

# Boxplot of the M values (one box for each array)
x11(); boxplot(MA_q$M ~ col(MA_q$M),   #col returns a matrix of integers indicating their column number
  names=rownames(MA_q$targets),
  col=colors, xlab="Sample",
  ylab="M-values (Log2ratios)")

#-------------------------------------------------------------------------------------------------------------

#EX.2

# Single-channel microarray preprocessing before proceeding with differential analysis

rm(list = ls())
library(affy) # BiocManager::install("affy")
library(limma)
#BiocManager::install("RSQLite")
#BiocManager::install("hgu133plus2cdf")
library(hgu133plus2cdf)


# For this tutorial we will consider the GSE29316 case-study
# We have three biological replicates of SHH stimulated (Sonic Hedgehog factor - involved 
#in organogenesis and tumorigenesis) primary colon myofibroblast cells CCD-18Co,
# and three replicates of unstimulated primary colon myofibroblast cells CCD-18Co.
# We will use these samples to analyze gene expression, but, differently from the previous case, 
# we work on one-channel microarray data

setwd("./../SHH") 
getwd()
targets <- readTargets("targets.txt")

#Now, we can produce an affy batch object that, similarly to the RG object for the dual channel data,
#will contain raw signal intensity for each probe and metadata of our single-channel microarray data.
aBatch <- ReadAffy(filenames=targets$filename,
                     sampleNames=targets$sample)

# It is possible to skip all the passages reported in the previous example by using the RMA algorithm. 
#This algorithm will perform correction and normalization automatically:
#1. Probe specific correction of the PM probes using a model based on observed intensity being the sum
#   of signal and noise
#2. Normalization of corrected PM probes using quantile normalization (Bolstad et al., 2003)
#3. Calculation of Expression measure using median polish.

# Build an expression set object from the affy batch by applying the RMA algorithm
eset <- rma(aBatch)


#Density histograms before and after the RMA processing

par(mfrow=c(2,1))
hist(aBatch)

library(affyPLM)
hist(eset, xlab='log intensity')

#hist(eset@assayData[["exprs"]])
#hist( exprs(eset), freq=F)

#PDF file of density histograms before and after the RMA processing

pdf("hist.pdf") 
par(mfrow=c(2,1))
hist(aBatch)
hist(eset, xlab='log intensity')
dev.off()

library(lattice)
df <- data.frame(Expression=as.vector(exprs(eset)),
                 Sample=sampleNames(eset)[col(eset)])
densityplot(~Expression|Sample, df, plot.points=F, ylab="Density")



#Pairwise MA-plots histogram before the RMA processing
x11(); MAplot(aBatch, pairs=T,
                plot.method="smoothScatter",
                main="M vs. A plot (aBatch)",
                cex=0.9)

#Pairwise MA-plots histogram after  the RMA processing
library(affyPLM)
x11(); MAplot(eset, pairs=T,
                plot.method="smoothScatter",
                main="M vs. A plot (eset)",
                cex=0.9)

# The experiment design type is: case-control
# We only have two conditions: SHH-treated and untreated fibroblasts
# Define the factors (i.e. the conditions we want to consider)
F <- factor(paste(targets$condition))
F

# Estimate fold changes and standard errors by fitting a linear model 
# for each gene, using the design matrix


#Build the design matrix
design.matrix <- model.matrix(~0+F)
colnames(design.matrix) <- levels(F)
# levels(F) are the values of the categorical variable F.

#Fitting the linear model
fit0 <- lmFit(eset, design.matrix) #Microarray linear model fit-class

#Build the contrasts 
cont.matrix <- makeContrasts(H="untreated - SHH", levels=design.matrix)
#This function expresses contrasts between a set of parameters as a numeric matrix. 
#The parameters are usually the coefficients from a linear model fit, so the matrix
#specifies which comparisons between the coefficients are to be extracted from the fit.

            
#Fitting contrasts
fit <- contrasts.fit(fit0, cont.matrix)
#the function re-orientates the fitted model object from the coefficients of the original
#design matrix to any set of contrasts of the original coefficients. The coefficients, 
#unscaled standard deviations and correlation matrix are re-calculated in terms of the contrasts.

#Apply the empirical Bayes method for smoothing of standard errors: 

fit_f <- eBayes(fit)
#These functions are used to rank genes in order of evidence for differential expression. 
#they compute statistics of differential expression using an empirical Bayes method to squeeze 
#the gene-wise residual variances towards a common value (standard errors that are much larger 
# or smaller than those from other genes are squeezed towards the average standard error

#Table of the top-ranked genes from the linear model fit
View(topTable(fit_f, coef=1, number=nrow(fit_f), sort.by="p"))

#We can also filter genes by their adjusted p-value (i.e. corrected for multiple testing - Benjamini Hochberg by defauly)
R <- topTable(fit_f, coef=1, number=nrow(fit_f))
L <- rownames(R[R$adj.P.Val < .01,])
length(rownames(R)) # or dim(R)[1]
length(L)
#L is the list of DEGs that passed the
# threshold: 804/54675 (about 1.47% of the total).
length(L)*100/dim(R)[1] 
L

#Print every probe in L vector in the limmaList.txt file
#Please, remember that this is a list of probe IDs. 
#It must be converted into a list of gene IDs (or gene symbols)
# Convert probe IDs into gene symbols
gene_symbols <- mapIds(
  hgu133plus2.db,
  keys = as.character(L),
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)
write(gene_symbols,"limmaList.txt", append=T, ncolumns=1)



