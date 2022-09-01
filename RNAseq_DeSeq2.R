#######################################

#__Date:__ April
#__Author:__ Austin Rosen
#__Script:__ RNAseq_DeSeq2.R
#__Project:__ To analyze RNA-seq in a project that looks at defense priming
#__Requires:__ 
# 
# + DESeq2 (1.20.0)   
# + corrplot (0.84)
# + RColorBrewer (1.1.2)
# + pheatmap (1.0.10)

######################################

######### FOR FIRST TIME USE ONLY ##############
######### After use, comment this section ##############

# Intall required packages/libraries:

#Install Bioconductor R version 3.5 or less
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")

#OR 

# Install bioconductor R version 3.6 or greater
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.10")

## When prompted to install (a/s) all or some of the dependent packages, move your cursor down to the console and type "a" for all

# Install DESeq2:
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

# Install 'apeglm'
BiocManager::install("apeglm")

# install corrplot:
install.packages("corrplot")

#Install pretty heatmap - pheatmap: https://cran.r-project.org/web/packages/pheatmap/index.html
install.packages("pheatmap")

#Install RColorBrewer: https://cran.r-project.org/web/packages/RColorBrewer/index.html
install.packages("RColorBrewer")

######### DONE WITH: FIRST TIME USE ONLY SECTION ##############
################################


###########  LOAD PACKAGE  #####################

# Do this every time.
# Load packages:
library(DESeq2)
library(corrplot)
library(RColorBrewer)
library(pheatmap)

#################################################




###########  READ IN THE DATA  #####################

# import the counts data
getwd()

#Set this to your working directory:
# You may need to set this to your own working directory to your scripts directory:
setwd("~/RNAseq_project/DeSeq2/02_scripts")
#setwd("/Users/erinnishimura/Dropbox/Courses/2019_DSCI512/05_courseContent/PROJ08_DESeq2/02_scripts")
getwd()
countsData <- read.table(file = "../01_input/Arabid_feature_counts_parameters.txt", header = FALSE, row.names = 1, skip = 2) # 

# :!: EXERCISE

# Explore the countsData object using head and dim and class

head(countsData)
dim(countsData)
class(countsData)

# Read in the metadata
metadata1 <- read.table(file = "../01_input/metadata_Arabid.txt", header = FALSE) # import the data
metadata1


# Organize the metadata file
colnames(metadata1) <- c("fasta", "plant", "sam", "name", "treatment", "rep")
metadata1

# Exercise: Looking at the metadata1 object and countsData objects (try head, dim, str, class, colnames):


# Organize the countsData file.
# Notice that the countsData file doesn't have any column headers:
head(countsData)

# Let's give countsData some columns names. The first names will be... chr', 'start', etc...
# The last names will be names for each sample. We can pull those names from metadata1:
as.vector(metadata1$name)



# Name countsData columns headers:
colnames(countsData) <- c("chr", "start", "stop", "strand", "length", as.vector(metadata1$name))

# :!: EXERCISE: Now look at the top of countsData again usign head():
head(countsData)





################### COUNT MATRIX INPUT ###################

# In this section we will prepare our input data for analysis.
# In the instruction for DESeq2, it states: "We read in a count matrix, which we will name cts, and the sample information table, which we will name coldata."

# OK, our task will be to generate a table called "cts" out of the countsData table.
# Subset the countsData 
head(countsData)
dim(countsData)
head(countsData[,6:23])


# Save just the subset as an object called cts:
cts <- as.matrix(countsData[,6:23])
head(cts)

# Yay, we made cts

# Next we need to make an information called coltable. We can make this out of the metadata table.

class(metadata1)
# Reorganize the metadata table so the names column are now row headers
metadata1
rownames(metadata1)<- metadata1$name
metadata1

coldata <- metadata1[,c("treatment", "rep")]
coldata$treatment <- as.factor(coldata$treatment)
coldata$rep <- as.factor(coldata$rep)
rownames(coldata)
colnames(cts)

# Yay! Now we have coldata! This is a new metadata object where we have just selected the type of information that is critical for deseq2 to use.

# One thing we need to explicitly check. The rownames of coldata need to exactly match the colnames of cts.
#Check that the names match  --> Should be TRUE
all(rownames(coldata) == colnames(cts))


# Next we will create an ddsHTSeq object out of cts and coldata:
# This will set a base design for your experiment:
# Load all the _counts.txt files and to attach them to the metadata.

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ treatment + rep)


################# PRE-FILTERING -- FILTER FOR PRESENT GENES: ################# 
# Not necessary, but helps keep things fast. 
# Exclude all samples that have less than 5 reads:
keep <- rowSums(counts(dds)) >= 0
dds <- dds[keep,]

# Exercise: How many did we exclude?
dim(dds)

str(dds)
coldata
################### NOTE ON FACTOR LEVELS ###################
# Organize the categories based on what makes sense:
dds$condition <- factor(dds$treatment, levels = c("C28","C28P3", "C28P3L7", "C28P3L7T3", "C35", "C35P3"))




# PERFORM DESEQ2 analysis:
####### This is where the magic happens! ########### 
# This will transform dds into a specialized object with many more fields filled in.
#dds <- DESeq(dds)

dds <- DESeq(dds, fitType = "local")

# Exercise: Check the output (class, str, dim, plotDispEsts) by adding the object dds inside the parentheses:
class(dds)
str(dds)
dim(dds)
plotDispEsts(dds)


# Here is a demonstration of the size Factor scaling that was calculated (sizeFactor):
dds$sizeFactor


# Demo: Access the normalized counts using counts(x, normalized = TRUE)
# Demo: Access the raw count info using counts(x, normalized = FALSE)
# Both of these datasets are good supplemental tables for papers
head(counts(dds, normalized = TRUE))
head(counts(dds, normalized = FALSE))



############## DIFFERENTIAL EXPRESSION ANALYSIS #####################

# calculate the statistically significant differences between E.coli and B. subtilis diets

# comparison A
res_C35vsC28 <- results(dds,
                          lfc = 0.5,
                          contrast=c("treatment","C35","C28"))

summary(res_C35vsC28)

# comparison B
res_C28P3vsC28 <- results(dds,
                        lfc = 0.5,
                        contrast=c("treatment","C28P3","C28"))

summary(res_C28P3vsC28)

# comparison C
res_C35P3vsC28 <- results(dds,
                          lfc = 0.5,
                          contrast=c("treatment","C35P3","C28"))

summary(res_C35P3vsC28)

# comparison D
res_C35P3vsC35 <- results(dds,
                          lfc = 0.5,
                          contrast=c("treatment","C35P3","C35"))

summary(res_C35P3vsC35)

# comparison E
res_C35P3vsC28P3 <- results(dds,
                          lfc = 0.5,
                          contrast=c("treatment","C35P3","C28P3"))

summary(res_C35P3vsC28P3)

# comparison F
res_C28P3L7vsC28 <- results(dds,
                          lfc = 0.5,
                          contrast=c("treatment","C28P3L7","C28"))

summary(res_C28P3L7vsC28)

# comparison G
res_C28P3L7vsC35 <- results(dds,
                            lfc = 0.5,
                            contrast=c("treatment","C28P3L7","C35"))

summary(res_C28P3L7vsC35)

# comparison H
res_C28P3L7T3vsC28 <- results(dds,
                            lfc = 0.5,
                            contrast=c("treatment","C28P3L7T3","C28"))

summary(res_C28P3L7T3vsC28)

# comparison I
res_C28P3L7T3vsC28P3 <- results(dds,
                              lfc = 0.5,
                              contrast=c("treatment","C28P3L7T3","C28P3"))

summary(res_C28P3L7T3vsC28P3)

# comparison J
res_C28P3L7T3vsC28P3L7 <- results(dds,
                                lfc = 0.5,
                                contrast=c("treatment","C28P3L7T3","C28P3L7"))

summary(res_C28P3L7T3vsC28P3L7)

# comparison K
res_C28P3L7T3vsC35P3 <- results(dds,
                                  lfc = 0.5,
                                  contrast=c("treatment","C28P3L7T3","C35P3"))

summary(res_C28P3L7T3vsC35P3)

#####################################
#Perform PCA
### Transform counts for data visualization
rld <- rlog(dds, blind=TRUE)

### Plot PCA 
plotPCA(rld, intgroup="treatment")







############## PERFORM LOG FOLD CHANGE SHRINKAGE FOR VISUALIZATION #####################

# An input requirement of the lfcShrink function is a coef term. This is pulled from the resultsNames of dds:


resultsNames(dds)

resLFC_EcolVBsubt <- lfcShrink(dds, 
                               coef="food_Ecoli_vs_Bsubtilis", res = res_EcolVBsubt, type='apeglm')



summary(resLFC_EcolVBsubt)


# Exercise 2: Inspect the results object called res_AcVuntreated and resLFC_AcVuntreated (head, dim, summary, str):


##################  Exploring and exporting results ##################  

##### KNOWN GENES:

# Check known genes to make sure everything is working as predicted. Check out :WBGene00018393 # aka msra-1

plotCounts(dds, gene=which(rownames(resLFC_EcolVBsubt) == "WBGene00018393"), intgroup="treatment")
plotCounts(dds, gene=which(rownames(resLFC_EcolVBsubt) == "WBGene00018393"), intgroup=c("treatment", "time"))

##### MA PLOTS:

# Plot the the defaul MA-plot and the shrunken MA-plot:

par(mfrow=c(1,1))
plotMA(res_EcolVBsubt, main="Acetic Acid versus untreated\nunshrunken", ylim = c(-7,7), 
       ylab = "log fold change (ratio of normalized E.coli / B.subtilis)",
       xlab = "means of normalized counts")

plotMA(resLFC_EcolVBsubt, main="Acetic Acid versus untreated\nshrunken", ylim = c(-7,7), 
       ylab = "log fold change (ratio of normalized E.coli / B.subtilis)",
       xlab = "means of normalized counts")

#Identify genes on the plot
#  Step1 -> execute idx code line below. 
#  Step2 -> Click on a dot in the plot. 
#  Step3 -> To finish, click on "finish" in the upper right hand corner of the plot
idx <- identify(resLFC_EcolVBsubt$baseMean, resLFC_EcolVBsubt$log2FoldChange)

#  Step4 -> click here to see what you got!
rownames(resLFC_EcolVBsubt)[idx]

# Cool thing to do -> look up these WBGene #'s on wormbase.org, then find the genes on the UCSC browser containing browser shots of the data:
# wormbase.org
# http://genome.ucsc.edu/s/Erin%20Osborne/GomezOrte_191210

# Cool thing to do.. plot one of your genes: Change "WBGene00010961" to "WBGeneYOURGENE"
plotCounts(dds, gene=which(rownames(resLFC_EcolVBsubt) == "WBGene00010961"), intgroup="treatment")


# Save the plots


pdf("../03_output/MAplots.pdf", height = 6, width = 10)

par(mfrow=c(1,2)) # This changes how many plot panels you can generate
plotMA(res_EcolVBsubt, main="Acetic Acid versus untreated\nunshrunken", ylim = c(-7,7), 
       ylab = "log fold change (ratio of normalized E.coli / B.subtilis)",
       xlab = "means of normalized counts")

plotMA(resLFC_EcolVBsubt, main="Acetic Acid versus untreated\nshrunken", ylim = c(-7,7), 
       ylab = "log fold change (ratio of normalized E.coli / B.subtilis)",
       xlab = "means of normalized counts")


dev.off() 
dev.off() #Note: sometimes you need 2x dev.off() lines of code to really truly escape out of the .pdf printer

##### CORRELATION MATRICES:

#Take r-stabilized log transformations of all the normalized count data. This will help with the problem that the data is noisy and it will help with the problem that the data is spread across a wide range of values.
rld <- rlog(dds, blind=FALSE)  #Take the r-stabilized log transformed data:

# Calculate the distances between each sample
sampleDists <- dist(t(assay(rld))) # calculate distance matrices:


sampleDistMatrix <- as.matrix(sampleDists) #convert from data.frame -> matrix
rownames(sampleDistMatrix) <- colnames(rld) # Add some labels
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) #Pick some pretty colors

# Draw the heatmap
par(mfrow=c(1,1))
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors, width = 4, height = 4) # Plot the heatmap


# Save the CORRELATION MATRIX as a .pdf file


pdf("../03_output/corr_matrix_plots.pdf", height = 6, width = 7)
par(mfrow=c(1,1))
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors) # Plot the heatmap

dev.off() # Do you see "null device"
dev.off() #-> sometimes you need a second dev.off until you see "null device"



##### VOLCANO PLOTS:###################

resultsNames(dds)
# Volcano plots are nice ways of displaying the fold change against the p-value.
res_volcano_plot <- lfcShrink(dds,coef="food_Ecoli_vs_Bsubtilis", type='apeglm') # Calculate padj without a lower lfc limit

significantLFC <- subset(resLFC_EcolVBsubt, padj < 0.05) # Identify significantly changing genes
significant_points_to_plot <-res_volcano_plot[which(rownames(res_volcano_plot) %in% rownames(significantLFC)),] 

# We will set the top limit of the plot as 100. I'll need to find all the points that exceed that measure:
maxedout <- subset(res_volcano_plot, padj < 10e-100)

#Draw plot:
par(mfrow=c(1,1)) # one plot only 

# Draw the plot
plot(res_volcano_plot$log2FoldChange, -log10(res_volcano_plot$padj),
     main="Volcano plot", xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)", pch=20, cex = 0.4, ylim = c(0, 100), xlim = c(-6,6), col = "#00000030")

# Add points
points(maxedout$log2FoldChange, rep(102, length(maxedout$log2FoldChange)), 
       pch=17, cex = 0.4, ylim = c(0, 100), col = "red")

points(significant_points_to_plot$log2FoldChange, -log10(significant_points_to_plot$padj),
       pch=20, cex = 0.4, ylim = c(0, 100), col = "red")

# Add lines
abline(v=0, col = "blue")
abline(v=0.5, col = "blue", lty = "dashed")
abline(v=-0.5, col = "blue", lty = "dashed")
abline(h=-log10(0.005), col = "blue", lty = "dashed")



## Save the VOLCANO PLOT as a .pdf
pdf("../03_output/volcanoplot.pdf", width = 6, height = 6)

par(mfrow=c(1,1))
plot(res_volcano_plot$log2FoldChange, -log10(res_volcano_plot$padj),
     main="Volcano plot", xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)",
     pch=20, cex = 0.4, ylim = c(0, 100), xlim = c(-4.3, 4.3),col = "#00000030")

points(maxedout$log2FoldChange, rep(102, length(maxedout$log2FoldChange)), 
       pch=17, cex = 0.4, col = "red")

points(significant_points_to_plot$log2FoldChange, -log10(significant_points_to_plot$padj),
       pch=20, cex = 0.4, col = "red")

abline(v=0, col = "blue")
abline(v=0.5, col = "blue", lty = "dashed")
abline(v=-0.5, col = "blue", lty = "dashed")
abline(h=-log10(0.005), col = "blue", lty = "dashed")

dev.off()
dev.off()



############## GETTING LISTS OF SIGNIFICANTLY CHANGING GENES #####################

# We will use the set of significantly changing genes from our variance shrunken analysis:
#Remember we calculated this above as:
#res_EcolVBsubt <- results(dds,
#                          lfc = 0.5,
#                          contrast=c("food", "Ecoli", "Bsubtilis"))

#resultsNames(dds)
#resLFC_AcVuntreated <- lfcShrink(dds, 
#                                 coef="treatment_aceticAcidTreated_vs_untreated", res = res_AcVuntreated)


# Check the results table:
summary(resLFC_EcolVBsubt)
head(resLFC_EcolVBsubt)


# Select the significant subset of genes that are up-regulated in acetic acid
Up_in_e.coli <- subset(resLFC_EcolVBsubt, padj < 0.05 & log2FoldChange > 0.5)
Up_in_e.coli <- Up_in_e.coli[order(Up_in_e.coli$padj),] #order them

head(Up_in_e.coli) # Check them
dim(Up_in_e.coli)

# Select the significant subset of genes that are down-regulated in acetic acid
Down_in_e.coli <- subset(resLFC_EcolVBsubt, padj < 0.05 & log2FoldChange < -0.5)
Down_in_e.coli <- Down_in_e.coli[order(Down_in_e.coli$padj),]

head(Down_in_e.coli)
dim(Down_in_e.coli)

# Save these lists to output files:
write(rownames(Up_in_e.coli), file = "../03_output/Genes Up in E.coli.txt", sep = "\n")
write(rownames(Down_in_e.coli), file = "../03_output/Genes Down in E.coli.txt", sep = "\n")
write(rownames(resLFC_EcolVBsubt), file = "../03_output/all_present_genes.txt", sep = "\n") # sometimes you need a control/background set of genes, too.


############## MAKE HIERARCHICALLY CLUSTERED HEATMAPS OF ALL CHANGING GENES #####################

# Let's loosen our restrictions on significance to all genes with any log fold change and adjusted p-values less than 0.1 (both are default)

# Get acetic acid v. untreated differentially expressed genes:
res_EcolVBsubt_2 <- results(dds,
                            contrast=c("food", "Ecoli", "Bsubtilis"))

# Get 20 v 15 C differentially expressed genes:
res_20v15_2 <- results(dds,
                       contrast=c("temp", "20", "15"))

# Get 25 v 15 C differentially expressed genes:
res_25v15_2 <- results(dds,
                       contrast=c("temp", "25", "15"))

# Get 25 v 20 C differentially expressed genes:
res_25v20_2 <- results(dds,
                       contrast=c("temp", "25", "20"))

#Subset each results table for just the differentially expressed genes:
sign_food <- subset(res_EcolVBsubt_2, padj < 0.05)
dim(subset(res_EcolVBsubt_2, padj < 0.05))

sign_20v15 <- subset(res_20v15_2, padj < 0.05)
sign_25v15 <- subset(res_25v15_2, padj < 0.05)
sign_25v20 <- subset(res_25v20_2, padj < 0.05)

#Determine how many genes were captured and merge them:
changing_genes <- rbind(sign_food, sign_20v15, sign_25v15, sign_25v20)

dim(changing_genes)
length(unique(rownames(changing_genes)))
# Get all r-stabilized log values of all the data:
rdl_all <- assay(rlog(dds, blind=FALSE))

#Subset just the changing genes:
changing_lrt_rdl <- subset(rdl_all, rownames(rdl_all) %in% rownames(changing_genes))
dim(changing_lrt_rdl)
head(changing_lrt_rdl)

# Make sure it is in matrix form:
class(changing_lrt_rdl)

# Draw a heat map
# Scale by row, listed below as scale = "row". This is really important. It sets the mean of every row to 0 and the standard deviation of every row to 1:
p <- pheatmap(changing_lrt_rdl, 
              scale="row", 
              color = colorRampPalette(c("blue", "white", "red"), space = "Lab")(100),
              cluster_rows=TRUE, 
              cluster_cols=FALSE, 
              clustering_distance_rows = "euclidean", 
              clustering_method = "complete",
              show_rownames = FALSE)

p
help(pheatmap)
# Tip: If euclidean doesn't look good, try "correlation" other clustering distances
# Tip: If complete doesn't look good, try other methods
# Save the clustered heatmap plot as a .pdf:

pdf("../03_output/clustered_genes.pdf", width = 6, height = 8)

p

dev.off()
dev.off()


# Get versions
sessionInfo()
#######################################

#__Date:__ December
#__Author:__ Erin Osborne Nishimura
#__Script:__ 191210_GomezOrte.R
#__Project:__ To analyze RNA-seq in a project that compares yeast grown in rich media v. acetic acid media at different time points
#__Requires:__ 
# 
# + DESeq2 (1.20.0)   
# + corrplot (0.84)
# + RColorBrewer (1.1.2)
# + pheatmap (1.0.10)

######################################

######### FOR FIRST TIME USE ONLY ##############
######### After use, comment this section ##############

# Intall required packages/libraries:

#Install Bioconductor R version 3.5 or less
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")

#OR 

# Install bioconductor R version 3.6 or greater
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.10")

## When prompted to install (a/s) all or some of the dependent packages, move your cursor down to the console and type "a" for all

# Install DESeq2:
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

# Install 'apeglm'
BiocManager::install("apeglm")

# install corrplot:
install.packages("corrplot")

#Install pretty heatmap - pheatmap: https://cran.r-project.org/web/packages/pheatmap/index.html
install.packages("pheatmap")

#Install RColorBrewer: https://cran.r-project.org/web/packages/RColorBrewer/index.html
install.packages("RColorBrewer")

######### DONE WITH: FIRST TIME USE ONLY SECTION ##############
################################


###########  LOAD PACKAGE  #####################

# Do this every time.
# Load packages:
library(DESeq2)
library(corrplot)
library(RColorBrewer)
library(pheatmap)

#################################################




###########  READ IN THE DATA  #####################

# import the counts data
getwd()

#Set this to your working directory:
# You may need to set this to your own working directory to your scripts directory:
setwd("/set/to/your/scripts/directory")
#setwd("/Users/erinnishimura/Dropbox/Courses/2019_DSCI512/05_courseContent/PROJ08_DESeq2/02_scripts")
getwd()
countsData <- read.table(file = "../01_input/counts.txt", header = FALSE, row.names = 1, skip = 2) # 

# :!: EXERCISE

# Explore the countsData object using head and dim and class

head(countsData)
dim(countsData)
class(countsData)

# Read in the metadata
metadata1 <- read.table(file = "../01_input/metadata_gomezOrte.txt", header = FALSE) # import the data
metadata1


# Organize the metadata file
colnames(metadata1) <- c("fasta1", "fasta2", "names1", "names2", "food", "temp", "rep")
metadata1

# Exercise: Looking at the metadata1 object and countsData objects (try head, dim, str, class, colnames):


# Organize the countsData file.
# Notice that the countsData file doesn't have any column headers:
head(countsData)

# Let's give countsData some columns names. The first names will be... chr', 'start', etc...
# The last names will be names for each sample. We can pull those names from metadata1:
as.vector(metadata1$names2)



# Name countsData columns headers:
colnames(countsData) <- c("chr", "start", "stop", "strand", "length", as.vector(metadata1$names2))

# :!: EXERCISE: Now look at the top of countsData again usign head():





################### COUNT MATRIX INPUT ###################

# In this section we will prepare our input data for analysis.
# In the instruction for DESeq2, it states: "We read in a count matrix, which we will name cts, and the sample information table, which we will name coldata."

# OK, our task will be to generate a table called "cts" out of the countsData table.
# Subset the countsData 
head(countsData)
dim(countsData)
head(countsData[,6:23])


# Save just the subset as an object called cts:
cts <- as.matrix(countsData[,6:23])
head(cts)

# Yay, we made cts

# Next we need to make an information called coltable. We can make this out of the metadata table.

class(metadata1)
# Reorganize the metadata table so the names2 column are now row headers
metadata1
rownames(metadata1)<- metadata1$names2
metadata1

coldata <- metadata1[,c("food", "temp", "rep")]
coldata$temp <- as.factor(coldata$temp)
coldata$rep <- as.factor(coldata$rep)
rownames(coldata)
colnames(cts)

# Yay! Now we have coldata! This is a new metadata object where we have just selected the type of information that is critical for deseq2 to use.

# One thing we need to explicitly check. The rownames of coldata need to exactly match the colnames of cts.
#Check that the names match  --> Should be TRUE
all(rownames(coldata) == colnames(cts))


# Next we will create an ddsHTSeq object out of cts and coldata:
# This will set a base design for your experiment:
# Load all the _counts.txt files and to attach them to the metadata.

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ food + temp)


################# PRE-FILTERING -- FILTER FOR PRESENT GENES: ################# 
# Not necessary, but helps keep things fast. 
# Exclude all samples that have 0 reads:
keep <- rowSums(counts(dds)) >= 1
dds <- dds[keep,]

# Exercise: How many did we exclude?
dim(dds)

str(dds)
coldata
################### NOTE ON FACTOR LEVELS ###################
# Organize the categories based on what makes sense:
dds$treatment <- factor(dds$food, levels = c("Ecoli","Bsubtilis"))
dds$time <- factor(dds$temp, levels = c("15", "20", "25"))



# PERFORM DESEQ2 analysis:
####### This is where the magic happens! ########### 
# This will transform dds into a specialized object with many more fields filled in.
dds <- DESeq(dds)



# Exercise: Check the output (class, str, dim, plotDispEsts) by adding the object dds inside the parentheses:
class()
str()
dim()
plotDispEsts()


# Here is a demonstration of the size Factor scaling that was calculated (sizeFactor):
dds$sizeFactor


# Demo: Access the normalized counts using counts(x, normalized = TRUE)
# Demo: Access the raw count info using counts(x, normalized = FALSE)
# Both of these datasets are good supplemental tables for papers
head(counts(dds, normalized = TRUE))
head(counts(dds, normalized = FALSE))



############## DIFFERENTIAL EXPRESSION ANALYSIS #####################

# calculate the statistically significant differences between E.coli and B. subtilis diets
res_EcolVBsubt <- results(dds,
                          lfc = 0.5,
                          contrast=c("food", "Ecoli", "Bsubtilis"))

summary(res_EcolVBsubt)

############## PERFORM LOG FOLD CHANGE SHRINKAGE FOR VISUALIZATION #####################

# An input requirement of the lfcShrink function is a coef term. This is pulled from the resultsNames of dds:


resultsNames(dds)

resLFC_EcolVBsubt <- lfcShrink(dds, 
                               coef="food_Ecoli_vs_Bsubtilis", res = res_EcolVBsubt, type='apeglm')



summary(resLFC_EcolVBsubt)


# Exercise 2: Inspect the results object called res_AcVuntreated and resLFC_AcVuntreated (head, dim, summary, str):


##################  Exploring and exporting results ##################  

##### KNOWN GENES:

# Check known genes to make sure everything is working as predicted. Check out :WBGene00018393 # aka msra-1

plotCounts(dds, gene=which(rownames(resLFC_EcolVBsubt) == "WBGene00018393"), intgroup="treatment")
plotCounts(dds, gene=which(rownames(resLFC_EcolVBsubt) == "WBGene00018393"), intgroup=c("treatment", "time"))

##### MA PLOTS:

# Plot the the defaul MA-plot and the shrunken MA-plot:

par(mfrow=c(1,1))
plotMA(res_EcolVBsubt, main="Acetic Acid versus untreated\nunshrunken", ylim = c(-7,7), 
       ylab = "log fold change (ratio of normalized E.coli / B.subtilis)",
       xlab = "means of normalized counts")

plotMA(resLFC_EcolVBsubt, main="Acetic Acid versus untreated\nshrunken", ylim = c(-7,7), 
       ylab = "log fold change (ratio of normalized E.coli / B.subtilis)",
       xlab = "means of normalized counts")

#Identify genes on the plot
#  Step1 -> execute idx code line below. 
#  Step2 -> Click on a dot in the plot. 
#  Step3 -> To finish, click on "finish" in the upper right hand corner of the plot
idx <- identify(resLFC_EcolVBsubt$baseMean, resLFC_EcolVBsubt$log2FoldChange)

#  Step4 -> click here to see what you got!
rownames(resLFC_EcolVBsubt)[idx]

# Cool thing to do -> look up these WBGene #'s on wormbase.org, then find the genes on the UCSC browser containing browser shots of the data:
# wormbase.org
# http://genome.ucsc.edu/s/Erin%20Osborne/GomezOrte_191210

# Cool thing to do.. plot one of your genes: Change "WBGene00010961" to "WBGeneYOURGENE"
plotCounts(dds, gene=which(rownames(resLFC_EcolVBsubt) == "WBGene00010961"), intgroup="treatment")


# Save the plots


pdf("../03_output/MAplots.pdf", height = 6, width = 10)

par(mfrow=c(1,2)) # This changes how many plot panels you can generate
plotMA(res_EcolVBsubt, main="Acetic Acid versus untreated\nunshrunken", ylim = c(-7,7), 
       ylab = "log fold change (ratio of normalized E.coli / B.subtilis)",
       xlab = "means of normalized counts")

plotMA(resLFC_EcolVBsubt, main="Acetic Acid versus untreated\nshrunken", ylim = c(-7,7), 
       ylab = "log fold change (ratio of normalized E.coli / B.subtilis)",
       xlab = "means of normalized counts")


dev.off() 
dev.off() #Note: sometimes you need 2x dev.off() lines of code to really truly escape out of the .pdf printer

##### CORRELATION MATRICES:

#Take r-stabilized log transformations of all the normalized count data. This will help with the problem that the data is noisy and it will help with the problem that the data is spread across a wide range of values.
rld <- rlog(dds, blind=FALSE)  #Take the r-stabilized log transformed data:

# Calculate the distances between each sample
sampleDists <- dist(t(assay(rld))) # calculate distance matrices:


sampleDistMatrix <- as.matrix(sampleDists) #convert from data.frame -> matrix
rownames(sampleDistMatrix) <- colnames(rld) # Add some labels
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255) #Pick some pretty colors

# Draw the heatmap
par(mfrow=c(1,1))
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors, width = 4, height = 4) # Plot the heatmap


# Save the CORRELATION MATRIX as a .pdf file


pdf("../03_output/corr_matrix_plots.pdf", height = 6, width = 7)
par(mfrow=c(1,1))
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors) # Plot the heatmap

dev.off() # Do you see "null device"
dev.off() #-> sometimes you need a second dev.off until you see "null device"



##### VOLCANO PLOTS:###################

resultsNames(dds)
# Volcano plots are nice ways of displaying the fold change against the p-value.
res_volcano_plot <- lfcShrink(dds,coef="food_Ecoli_vs_Bsubtilis", type='apeglm') # Calculate padj without a lower lfc limit

significantLFC <- subset(resLFC_EcolVBsubt, padj < 0.05) # Identify significantly changing genes
significant_points_to_plot <-res_volcano_plot[which(rownames(res_volcano_plot) %in% rownames(significantLFC)),] 

# We will set the top limit of the plot as 100. I'll need to find all the points that exceed that measure:
maxedout <- subset(res_volcano_plot, padj < 10e-100)

#Draw plot:
par(mfrow=c(1,1)) # one plot only 

# Draw the plot
plot(res_volcano_plot$log2FoldChange, -log10(res_volcano_plot$padj),
     main="Volcano plot", xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)", pch=20, cex = 0.4, ylim = c(0, 100), xlim = c(-6,6), col = "#00000030")

# Add points
points(maxedout$log2FoldChange, rep(102, length(maxedout$log2FoldChange)), 
       pch=17, cex = 0.4, ylim = c(0, 100), col = "red")

points(significant_points_to_plot$log2FoldChange, -log10(significant_points_to_plot$padj),
       pch=20, cex = 0.4, ylim = c(0, 100), col = "red")

# Add lines
abline(v=0, col = "blue")
abline(v=0.5, col = "blue", lty = "dashed")
abline(v=-0.5, col = "blue", lty = "dashed")
abline(h=-log10(0.005), col = "blue", lty = "dashed")



## Save the VOLCANO PLOT as a .pdf
pdf("../03_output/volcanoplot.pdf", width = 6, height = 6)

par(mfrow=c(1,1))
plot(res_volcano_plot$log2FoldChange, -log10(res_volcano_plot$padj),
     main="Volcano plot", xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)",
     pch=20, cex = 0.4, ylim = c(0, 100), xlim = c(-4.3, 4.3),col = "#00000030")

points(maxedout$log2FoldChange, rep(102, length(maxedout$log2FoldChange)), 
       pch=17, cex = 0.4, col = "red")

points(significant_points_to_plot$log2FoldChange, -log10(significant_points_to_plot$padj),
       pch=20, cex = 0.4, col = "red")

abline(v=0, col = "blue")
abline(v=0.5, col = "blue", lty = "dashed")
abline(v=-0.5, col = "blue", lty = "dashed")
abline(h=-log10(0.005), col = "blue", lty = "dashed")

dev.off()
dev.off()



############## GETTING LISTS OF SIGNIFICANTLY CHANGING GENES #####################

# We will use the set of significantly changing genes from our variance shrunken analysis:
#Remember we calculated this above as:
#res_EcolVBsubt <- results(dds,
#                          lfc = 0.5,
#                          contrast=c("food", "Ecoli", "Bsubtilis"))

#resultsNames(dds)
#resLFC_AcVuntreated <- lfcShrink(dds, 
#                                 coef="treatment_aceticAcidTreated_vs_untreated", res = res_AcVuntreated)


# Check the results table:
summary(resLFC_EcolVBsubt)
head(resLFC_EcolVBsubt)


# Select the significant subset of genes that are up-regulated in acetic acid
Up_in_e.coli <- subset(resLFC_EcolVBsubt, padj < 0.05 & log2FoldChange > 0.5)
Up_in_e.coli <- Up_in_e.coli[order(Up_in_e.coli$padj),] #order them

head(Up_in_e.coli) # Check them
dim(Up_in_e.coli)

# Select the significant subset of genes that are down-regulated in acetic acid
Down_in_e.coli <- subset(resLFC_EcolVBsubt, padj < 0.05 & log2FoldChange < -0.5)
Down_in_e.coli <- Down_in_e.coli[order(Down_in_e.coli$padj),]

head(Down_in_e.coli)
dim(Down_in_e.coli)

# Save these lists to output files:
write(rownames(Up_in_e.coli), file = "../03_output/Genes Up in E.coli.txt", sep = "\n")
write(rownames(Down_in_e.coli), file = "../03_output/Genes Down in E.coli.txt", sep = "\n")
write(rownames(resLFC_EcolVBsubt), file = "../03_output/all_present_genes.txt", sep = "\n") # sometimes you need a control/background set of genes, too.


############## MAKE HIERARCHICALLY CLUSTERED HEATMAPS OF ALL CHANGING GENES #####################

# Let's loosen our restrictions on significance to all genes with any log fold change and adjusted p-values less than 0.1 (both are default)

# Get acetic acid v. untreated differentially expressed genes:
res_EcolVBsubt_2 <- results(dds,
                            contrast=c("food", "Ecoli", "Bsubtilis"))

# Get 20 v 15 C differentially expressed genes:
res_20v15_2 <- results(dds,
                       contrast=c("temp", "20", "15"))

# Get 25 v 15 C differentially expressed genes:
res_25v15_2 <- results(dds,
                       contrast=c("temp", "25", "15"))

# Get 25 v 20 C differentially expressed genes:
res_25v20_2 <- results(dds,
                       contrast=c("temp", "25", "20"))

#Subset each results table for just the differentially expressed genes:
sign_food <- subset(res_EcolVBsubt_2, padj < 0.05)
dim(subset(res_EcolVBsubt_2, padj < 0.05))

sign_20v15 <- subset(res_20v15_2, padj < 0.05)
sign_25v15 <- subset(res_25v15_2, padj < 0.05)
sign_25v20 <- subset(res_25v20_2, padj < 0.05)

#Determine how many genes were captured and merge them:
changing_genes <- rbind(sign_food, sign_20v15, sign_25v15, sign_25v20)

dim(changing_genes)
length(unique(rownames(changing_genes)))
# Get all r-stabilized log values of all the data:
rdl_all <- assay(rlog(dds, blind=FALSE))

#Subset just the changing genes:
changing_lrt_rdl <- subset(rdl_all, rownames(rdl_all) %in% rownames(changing_genes))
dim(changing_lrt_rdl)
head(changing_lrt_rdl)

# Make sure it is in matrix form:
class(changing_lrt_rdl)

# Draw a heat map
# Scale by row, listed below as scale = "row". This is really important. It sets the mean of every row to 0 and the standard deviation of every row to 1:
p <- pheatmap(changing_lrt_rdl, 
              scale="row", 
              color = colorRampPalette(c("blue", "white", "red"), space = "Lab")(100),
              cluster_rows=TRUE, 
              cluster_cols=FALSE, 
              clustering_distance_rows = "euclidean", 
              clustering_method = "complete",
              show_rownames = FALSE)

p
help(pheatmap)
# Tip: If euclidean doesn't look good, try "correlation" other clustering distances
# Tip: If complete doesn't look good, try other methods
# Save the clustered heatmap plot as a .pdf:

pdf("../03_output/clustered_genes.pdf", width = 6, height = 8)

p

dev.off()
dev.off()


# Get versions
sessionInfo()
