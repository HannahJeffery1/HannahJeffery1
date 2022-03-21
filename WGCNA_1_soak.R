#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================

rm(list=ls())

# Display the current working directory
# Load the WGCNA package
install.packages("WGCNA")

BiocManager::install("Biobase")
BiocManager::install("GO.db")
BiocManager::install("impute")
BiocManager::install("preprocessCore")
BiocManager::valid("GO.db")
BiocManager::valid("impute")
BiocManager::valid("preprocessCore")

library(WGCNA);
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);


#WGCNA_file = read.csv("WGCNA_norm5.6.csv");
#WGCNA_file = read.csv("WGCNA_norm10.6_soak.csv");
#WGCNA_file = read.csv("soakDEGs_2.csv");

library(readxl)
#WGCNA_file = read_excel("soakDEGs_2_yellow.xlsx");
#WGCNA_file = read_excel("soakDEGs_2_brown.xlsx");
WGCNA_file = read.csv("WGCNA_norm10.6_soak.csv");
#WGCNA_file = read.csv("WGCNA_norm10.6_soak_yellow.csv");
#WGCNA_file = read.csv("WGCNA_norm10.6_soak_brown.csv");
# Take a quick look at what is in the data set:
dim(WGCNA_file);
names(WGCNA_file);
View(WGCNA_file)

#=====================================================================================
#
#  Code chunk 2
#

#=====================================================================================

##### This step was already performed.  The WGCNA_norm file is the product of this data processing step.

# Create a table of normalized gene counts using the element from "treatments" in DESeq2.contrasts
# and the bamfiles.txt file in Samples_new_org > HiSat2_polya_bam_outputs.

library("limma")
library("edgeR")
BiocManager::install("Rsubread")
library("Rsubread")

samps <- factor(treatments)
mat<-model.matrix(~0 + samps)
y <- voom(samps, mat, plot = T) # Check the count data for heteroskedasticity
targets<- readTargets(file = "bamfiles.txt")

fc <- featureCounts(files=targets$FeatureCountsFile,
                    annot.ext="Pvulgaris_442_v2.1.gene.gtf",
                    isGTFAnnotationFile=TRUE,GTF.featureType="exon",
                    GTF.attrType = "gene_id",
                    chrAliases=NULL,useMetaFeatures=TRUE,allowMultiOverlap=FALSE,
                    readExtension5=0,
                    readExtension3=0,read2pos=NULL,countMultiMappingReads=FALSE,
                    minMQS=0,
                    ignoreDup=FALSE,strandSpecific=0,checkFragLength=FALSE,
                    minFragLength=50,
                    maxFragLength=600,countChimericFragments=TRUE,nthreads=4,
                    reportReads=NULL)
x <- DGEList(counts=fc$counts, genes=fc$annotation[,c("GeneID","Length")])

# Filtering parameter
isexpr <- rowSums(cpm(x) > 10) >= 6
x <- x[isexpr,]
x_rpkm <- rpkm(x,x$genes$Length)
write.table(x_rpkm,file="RPKM_Normalized_Gene_Counts10.6.txt")
write.csv(RPKM_Normalized_Gene_Counts10.6, file = "RPKM_Normalized_Gene_Counts10.6.csv")
# Average the columns by hand in Excel.

#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================

View(WGCNA_file)
datExprbean = as.data.frame(t(WGCNA_file[, -c(1:3)]))
View(datExprbean)

gsgbean = goodSamplesGenes(datExprbean, verbose = 3);
# The package automatically deducts genes that are not often expressed.
gsgbean$allOK
View(gsgbean)

#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================
View(datExprbean)

if (!gsgbean$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsgbean$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExprbean)[!gsgbean$goodGenes], collapse = ", ")));
  if (sum(!gsgbean$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExprbean)[!gsgbean$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExprbean = datExprbean[gsgbean$goodSamples, gsgbean$goodGenes]
}
gsgbean = goodSamplesGenes(datExprbean, verbose = 3);
gsgbean$allOK

#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


beansampleTree = hclust(dist(datExprbean), method = "average");
#View(beansampleTree)
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(beansampleTree, main = "Sample Clustering to Detect Outliers (>60,000)", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)


#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================


# Plot a line to show the cut
abline(h = 60000, col = "red");
# Determine cluster under the line
clust = cutreeStatic(beansampleTree, cutHeight = 60000, minSize = 100)
table(clust)
# clust 1 contains the samples we want to keep.
#keepSamples = (clust==1)
#datExprbean = datExprbean[keepSamples, ]
#beannGenes = ncol(datExprbean)
#beannSamples = nrow(datExprbean)


#=====================================================================================
#
#  Code chunk 7
#
#=====================================================================================

beantraitData = read.csv("traits_file_soak_tog.csv");
#beantraitData = read_excel("traits_file_soak_sep.xlsx")
#beantraitData = read_excel("traits_file_soak_yellow.xlsx")
#beantraitData = read_excel("traits_file_soak_brown.xlsx")

dim(beantraitData)
names(beantraitData)
View(beantraitData)

# Form a data frame analogous to expression data that will hold the phenotypic traits.
View(datExprbean)
beanSamples = rownames(datExprbean);
beantraitRows = match(beanSamples, beantraitData$Sample);
beandatTraits = beantraitData[beantraitRows, -1];
beandatTraits$Soak_time_hr<-as.numeric(beandatTraits$Soak_time_hr)
beandatTraits$Cook_time_min<-as.numeric(beandatTraits$Cook_time_min)
#beandatTraits$SoakT_Y<-as.numeric(beandatTraits$SoakT_Y)
#beandatTraits$SoakT_B<-as.numeric(beandatTraits$SoakT_B)
#beandatTraits$CookT_Y<-as.numeric(beandatTraits$CookT_Y)
#beandatTraits$CookT_B<-as.numeric(beandatTraits$CookT_B)
#rownames(beandatTraits) = beantraitData[beantraitRows, 1];

View(beanSamples)
View(beantraitRows)
View(beandatTraits)
rownames(beandatTraits)

collectGarbage();


#=====================================================================================
#
#  Code chunk 8
#
#=====================================================================================


# Re-cluster samples
beansampleTree2 = hclust(dist(datExprbean), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
beantraitColors = numbers2colors(beandatTraits, signed = FALSE);
#beantraitColors = numbers2colors(beandatTraits, signed = TRUE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(beansampleTree2, beantraitColors,
                    groupLabels = names(beandatTraits), 
                    main = "Sample Dendrogram and Trait Heatmap")


#=====================================================================================
#
#  Code chunk 9
#
#=====================================================================================

#save(datExprbean, beandatTraits, file = "Bean_cooking_data_WGCNA_5.6.RData")
save(datExprbean, beandatTraits, file = "Bean_cooking_data_WGCNA_10.6.RData")
#save(datExprbean, beandatTraits, file = "Bean_cooking_data_WGCNA_10.6_yellow.RData")
#save(datExprbean, beandatTraits, file = "Bean_cooking_data_WGCNA_10.6_brown.RData")

