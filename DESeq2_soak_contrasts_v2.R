#Credit: https://lashlock.github.io/compbio/R_presentation.html

rm(list=ls()) 
# Remove stuff from the R environment
# Set working directory to HTSeq_Outputs_alignments

myFiles <- list.files(pattern=glob2rx("*.htseq.txt"), 
                      all.files=T, 
                      full.names=T)
# Create a vector containing the names of your files
# Raw read counts are fine for DESeq2, but not for WGCNA

myFiles 
# myFiles is the list of files in your working directory which 
# match the pattern *.htseq.txt.

myData <- lapply(myFiles, 
                 read.table, colClasses=c("character","numeric"))
# Use list apply to apply the native R function read.table() to myFiles.

names(myData) <- c("x025x",
                   "x026x","x027x","x028x","x029x","x030x","x031x","x032x",
                   "x033x","x034x","x035x","x036x",
                   "x037x","x038x","x039x",
                   "x040x","x041x","x042x","x043x","x044x","x045x","x046x",
                   "x047x","x048x","x049x","x050x","x051x","x052x","x053x",
                   "x054x","x055x","x056x","x057x","x058x","x059x","x060x",
                   "x061x","x062x","x063x","x064x","x065x","x066x","x067x",
                   "x068x","x069x","x070x","x071x","x072x","x073x","x074x",
                   "x075x","x076x","x077x","x078x","x079x","x080x","x081x",
                   "x082x","x083x","x084x")

summary(myData) 
# Check out information about the list we created.
nrow(myData$x025x) 
# There should be 27,438 genes

head(myData$x025x) 
# Look at first ten lines of the first file

tail(myData$x025x) 
# Look at the last ten lines of the first file. 

# Only apply this function if there is any metadata at the end of the files.
# If there is metadata, it must be removed, or it will be treated as
# gene counts.
rmLines <- function(x){ 
  x <- x[1:27433,]
}
# This is a really quick and dirty solution. This is a function to select 
# the first 42,189 lines.


myData <- lapply(myData, rmLines) 
# Apply the rmLines() function to each object in the myData list.

nrow(myData$x025x) 
# Check to see how many lines there are in each element of myData.
# There will be fewer lines if you remove the metadata.
tail(myData$x025x)

# Generate combined data frame for counts data
counts <- data.frame(row.names=myData$x025x[,1],
                     x025x=myData$x025x[,2], 
                     x026x=myData$x026x[,2],
                     x027x=myData$x027x[,2],
                     x028x=myData$x028x[,2],
                     x029x=myData$x029x[,2],
                     x030x=myData$x030x[,2], 
                     x031x=myData$x031x[,2], 
                     x032x=myData$x032x[,2], 
                     x033x=myData$x033x[,2],
                     x034x=myData$x034x[,2], 
                     x035x=myData$x035x[,2], 
                     x036x=myData$x036x[,2],
                     x037x=myData$x037x[,2],
                     x038x=myData$x038x[,2],
                     x039x=myData$x039x[,2],
                     x040x=myData$x040x[,2],
                     x041x=myData$x041x[,2], 
                     x042x=myData$x042x[,2], 
                     x043x=myData$x043x[,2], 
                     x044x=myData$x044x[,2],
                     x045x=myData$x045x[,2], 
                     x046x=myData$x046x[,2], 
                     x047x=myData$x047x[,2], 
                     x048x=myData$x048x[,2],
                     x049x=myData$x049x[,2],
                     x050x=myData$x050x[,2],
                     x051x=myData$x051x[,2], 
                     x052x=myData$x052x[,2], 
                     x053x=myData$x053x[,2], 
                     x054x=myData$x054x[,2],
                     x055x=myData$x055x[,2],
                     x056x=myData$x056x[,2], 
                     x057x=myData$x057x[,2],
                     x058x=myData$x058x[,2],
                     x059x=myData$x059x[,2],
                     x060x=myData$x060x[,2],
                     x061x=myData$x061x[,2], 
                     x062x=myData$x062x[,2], 
                     x063x=myData$x063x[,2], 
                     x064x=myData$x064x[,2],
                     x065x=myData$x065x[,2],
                     x066x=myData$x066x[,2], 
                     x067x=myData$x067x[,2],
                     x068x=myData$x068x[,2],
                     x069x=myData$x069x[,2],
                     x070x=myData$x070x[,2],
                     x071x=myData$x071x[,2], 
                     x072x=myData$x072x[,2], 
                     x073x=myData$x073x[,2], 
                     x074x=myData$x074x[,2],
                     x075x=myData$x075x[,2], 
                     x076x=myData$x076x[,2],
                     x077x=myData$x077x[,2],
                     x078x=myData$x078x[,2],
                     x079x=myData$x079x[,2],
                     x080x=myData$x080x[,2], 
                     x081x=myData$x081x[,2], 
                     x082x=myData$x082x[,2], 
                     x083x=myData$x083x[,2],
                     x084x=myData$x084x[,2]
)

head(counts)

geneTotals <- rowSums(counts) 
# Evaluate the sum of each row and save 
# to a new vector called geneTotals.

countsNonZero <- counts[geneTotals>0,] 
# Subset the rows where the geneTotal is greater than 0.

countsNonZero

nrow(countsNonZero) # See how many genes are left in the analysis.
# 24378 genes remain in the full set of samples
# 23689 genes remain in the developmental subset
# 23350 genes remain in the soaking subset

# Redefine the names.
names(myData) <- c(rep("TZ27_0hr", 3),
                   rep("TZ27_3hr", 3),
                   rep("TZ27_6hr", 3),
                   rep("TZ27_12hr", 3),
                   rep("TZ27_18hr", 3),
                   rep("TZ37_0hr", 3),
                   rep("TZ37_3hr", 3),
                   rep("TZ37_6hr", 3),
                   rep("TZ37_12hr", 3),
                   rep("TZ37_18hr", 3),
                   rep("PI_0hr", 3),
                   rep("PI_3hr", 3),
                   rep("PI_6hr", 3),
                   rep("PI_12hr", 3),
                   rep("PI_18hr", 3),
                   rep("Erv_0hr", 3),
                   rep("Erv_3hr", 3),
                   rep("Erv_6hr", 3),
                   rep("Erv_12hr", 3),
                   rep("Erv_18hr", 3))



names(myData)
# Distinguish between fast and slow cooking genotypes at different soaking times.
# The annotated_genes_2 subset did not have this step.
# The annotated_genes_4 subset did have this step.  Do not use it.
treatments <- as.factor(c(rep("TZ27_0hr", 3),
                          rep("TZ27_3hr", 3),
                          rep("TZ27_6hr", 3),
                          rep("TZ27_12hr", 3),
                          rep("TZ27_18hr", 3),
                          rep("TZ37_0hr", 3),
                          rep("TZ37_3hr", 3),
                          rep("TZ37_6hr", 3),
                          rep("TZ37_12hr", 3),
                          rep("TZ37_18hr", 3),
                          rep("PI_0hr", 3),
                          rep("PI_3hr", 3),
                          rep("PI_6hr", 3),
                          rep("PI_12hr", 3),
                          rep("PI_18hr", 3),
                          rep("Erv_0hr", 3),
                          rep("Erv_3hr", 3),
                          rep("Erv_6hr", 3),
                          rep("Erv_12hr", 3),
                          rep("Erv_18hr", 3)))

treatments

library("BiocManager") # load the needed packages
library("limma")
library("edgeR")
library("DESeq2")

colData <- DataFrame(treatments) 

dds <- DESeqDataSetFromMatrix(countsNonZero, 
                              colData=colData, 
                              design=formula(~treatments)) 
# Create the DESeq object from the counts matrix we made that has 
# no non-expressed genes, uses the new colData object we made to 
# label the tissues, and the same experimental design using treatments 
# as the nominal categories.

dds <- DESeq(dds)

dds_pv <- results(dds)

dds_pv <- dds_pv[!(is.na(dds_pv$padj)),] 

#######################################################

# Create a DGEList object for edgeR
dge_1 <- DGEList(counts=countsNonZero, 
                 group=treatments, 
                 remove.zeros=T) 
# We won't be using the dge_1 object in our analysis because this 
# function removes all of the non-expressed genes from the subset 
# before plotting the data, but we will use just for making the plot.

# Generate MDS plots
plotMDS(dge_1,main="Log Fold Change in Gene Expression between Treatments", 
        labels=treatments,
        #gene.selection = "pairwise",
        gene.selection = "common",
        col=c(rep("#A46424", 15), 
              rep("#ECA524", 15),
              rep("#666B25", 15),
              rep("#E4CD24", 15)))

plotMDSlibrary("DESeq2")

colData_1 <- DataFrame(treatments) 
# 'Treatments' contains our column labels.  It is organized as factors.

dds_1 <- DESeqDataSetFromMatrix(countsNonZero, 
                                colData = colData_1,
                                design=formula(~treatments)) 

# Create the DESeq object from the counts matrix we made that 
# has no non-expressed genes, uses the new colData object we made 
# to label the tissues, and the same experimental design using 
# treatments as the nominal categories.

dds_1 <- DESeq(dds_1)

#Make PCA Plot
library(ggplot2)

tdds <- DESeqTransform(dds)
pcaData<-plotPCA(tdds, intgroup = "treatments", returnData = TRUE)
percentVar <- round(100*attr(pcaData, "percentVar"))
ggplot(pcaData, aes(PC1, PC2, color=treatments))+
  geom_point(size=4)+
  labs(title="Principle Component Analysis of All Replicates\n")+
  theme(plot.title = element_text(hjust = 0.5))+
  guides(color=guide_legend("Treatments"))+
  xlab(paste0("PC1: ",percentVar[1],"% variance"))+
  ylab(paste0("PC2: ",percentVar[2],"% variance"))+
  geom_text(aes(label=treatments), color="black", hjust=.5, 
            vjust=-.8, size=2.5)

# Save PCA data for future analysis
PCA<-plotPCA(tdds_1, intgroup = "treatments", returnData = TRUE)
write.csv(PCA, file = "PCA.csv")

plotDispEsts(dds_1) 
# Check the trend of normalized counts 
# to see how it corresponds to the dispersion estimate.

treatmentsbeans <- dds_1@rowRanges@elementMetadata@listData$treatments

geneDispersion <- dds_1@rowRanges@elementMetadata@listData$dispGeneEst

plot(treatmentsbeans, # See how fold change corresponds to dispersion.
     geneDispersion, 
     xlab="Fold Change <Treatment #1> Expression Relative to <Treatment #2>", 
     ylab="Gene Dispersion Estimate", pch=19, cex=0.4)
# Compare the data to see if the data is strangly dispersed or not 
# (the dispersal plots should lie as close to zero along the y-axis as possible
# and as close to zero along the x-axis as possible.

install.packages("gplots")
library("gplots")

rdds_1 <- rlog(dds_1) 
# Apply regularized logarithmic transformation to the gene counts.

# Using the assay() function in DESeq you can 
# extract the normalized count matrix, then use 
# t() to transpose the matrix and calculate 
# distance between rows of the matrix with dist()
distsRL_1 <- dist(t(assay(rdds_1))) 

mat_1 <- as.matrix(distsRL_1)

rownames(mat_1)<- colnames(mat_1) <- with(colData(dds_1), 
                                          treatments)

View(dds_1)
# Heatmap for soaking samples
heatmap.2(mat_1, 
          main = "Heat Map of All Replicates",
          key.xlab = "Dissimilarity Score",
          trace="none",
          scale = "none",
          keysize = 1.2,
          breaks = seq(0,100,1),
          #col=redblue(100), 
          dendrogram="none", 
          #density.info="none",
          margins=c(10,10),
          cexRow = 1,
          cexCol = 1,
          Rowv = FALSE,
          Colv = FALSE)

# Pearson's correlation coefficients
library("Hmisc")
cor_matrix_pearson <-rcorr(dds@assays@data$counts[dds_pv$padj < 0.05,], 
                           type=c("pearson"))
print(cor_matrix_pearson)
cor_matrix_spearman <-rcorr(dds@assays@data$counts[dds_pv$padj < 0.05,], 
                            type=c("spearman"))
print(cor_matrix_spearman)
# We used Spearman's correlation coefficients due to the data having a
# right-tailed skew.
# The standard deviation of "corResult" will sometimes be zero when
# comparing replicates.

# Make a clustered heat map of the differentially expressed genes
dds_pv_1 <- results(dds_1)

dds_pv_1 <- dds_pv_1[!(is.na(dds_pv_1$padj)),] 
# Remove rows with "NA" from the adjusted P-value column, 
# which DESeq inputs when it detects extreme outliers in samples 
# using Cook's distance calculation. These may have a significant 
# effect on the model, or, perhaps the model we used (~treatments) 
# does not account for some property of the data...
View(dds_pv_1$padj < 0.05)
View(dds_1@assays@data$counts[dds_pv_1$padj < 0.05,])
heatmap(dds_1@assays@data$counts[dds_pv_1$padj < 0.05,], 
        labCol=c("TZ27_0hr", "TZ27_0hr", "TZ27_0hr",
                 "TZ27_3hr", "TZ27_3hr", "TZ27_3hr",
                 "TZ27_6hr", "TZ27_6hr", "TZ27_6hr",
                 "TZ27_12hr", "TZ27_12hr", "TZ27_12hr",
                 "TZ27_18hr", "TZ27_18hr", "TZ27_18hr", 
                 "TZ37_0hr", "TZ37_0hr", "TZ37_0hr",
                 "TZ37_3hr", "TZ37_3hr", "TZ37_3hr",
                 "TZ37_6hr", "TZ37_6hr", "TZ37_6hr",
                 "TZ37_12hr", "TZ37_12hr", "TZ37_12hr",
                 "TZ37_18hr", "TZ37_18hr", "TZ37_18hr",
                 "PI_0hr", "PI_0hr", "PI_0hr",
                 "PI_3hr", "PI_3hr", "PI_3hr",
                 "PI_6hr", "PI_6hr", "PI_6hr",
                 "PI_12hr", "PI_12hr", "PI_12hr",
                 "PI_18hr", "PI_18hr", "PI_18hr", 
                 "Ervilha_0hr", "Ervilha_0hr", "Ervilha_0hr",
                 "Ervilha_3hr", "Ervilha_3hr", "Ervilha_3hr",
                 "Ervilha_6hr", "Ervilha_6hr", "Ervilha_6hr",
                 "Ervilha_12hr", "Ervilha_12hr", "Ervilha_12hr",
                 "Ervilha_18hr", "Ervilha_18hr", "Ervilha_18hr"),
        Colv = NA, Rowv = NA,
        cexRow=0.5, cexCol=1.3, margins=c(8,8),
        main = "Gene Expression Heatmap",
        col=heat.colors(5))
legend(x = "left", legend=c("min", "mean", "max"),
       fill=heat.colors(5))
# We can also plot the matrix of counts to see genes with 
# similar expression values in the conditions. Here we are 
# just plotting the heatmap of counts, not covariance between 
# genes across tissue type.

#######################################################

# Create an excel spreadsheet that contains all of the differentially 
# expressed genes with an adjusted p-value that is less than 0.05.
sigGenes_1 <- as.data.frame(dds_pv[dds_pv$padj < 0.05,])
write.csv(sigGenes_1, file = "sigGenes_4.csv")
# Modify the spreadsheet in Excel.
sigGenes_4<-read.csv(file = "sigGenes_4.csv")
colnames(sigGenes_4)<-c("V2", "Avg # of Transcripts", "log(Difference in Expression)", 
                        "Standard Error", "Wald's Stat", "p-Value", "Adjusted p-Value")
#View(sigGenes_4)

# Modify the annotation file
#View(Pvulgaris_442_v2.1.annotation_info_dedup)
df<-Pvulgaris_442_v2.1.annotation_info_dedup
df$V1<-NULL
df$V3<-NULL
df$V4<-NULL
df$V5<-NULL
df$V6<-NULL
df$V7<-NULL
df$V8<-NULL
df$V9<-NULL
df$V10<-NULL
df$V11<-NULL
df$V12<-NULL
View(df)

# Some steps were performed prior to creating the annotated data set:
# 1) Deleted the first column that only contained numbers.
# 2) Removed the ending from all of the gene names (i.e. ".v2.1").
# 3) Saved the file manually as a text file.
# 4) Imported the data file manaully at a base text file.

# Create the list of annotated, differentially expressed genes
annotated_sigGenes_1<-merge(df, sigGenes_4)
View(annotated_sigGenes_1)
write.table(annotated_sigGenes_1, file = "annotated_sigGenes_soak_4.txt")
write.csv(annotated_sigGenes_1, file = "annotated_sigGenes_soak_4.csv")
