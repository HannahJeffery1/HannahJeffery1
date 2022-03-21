#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================

rm(list=ls())

library(WGCNA)
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the data saved in the first part
#beanlnames = load(file = "Bean_cooking_data_WGCNA_5.6.RData");
beanlnames = load(file = "Bean_cooking_data_WGCNA_10.6.RData");
#beanlnames = load(file = "Bean_cooking_data_WGCNA_10.6_yellow.RData");
#beanlnames = load(file = "Bean_cooking_data_WGCNA_10.6_brown.RData");
#The variable lnames contains the names of loaded variables.
beanlnames

#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=35, by=2))
# Call the network topology analysis function
beansft = pickSoftThreshold(datExprbean, powerVector = powers, 
                            verbose = 5, networkType = "signed")
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(beansft$fitIndices[,1], -sign(beansft$fitIndices[,3])*beansft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"),
     ylim = c(-1,1));
text(beansft$fitIndices[,1], -sign(beansft$fitIndices[,3])*beansft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(beansft$fitIndices[,1], beansft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(beansft$fitIndices[,1], beansft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# For the soak subset:
        # The optimum power for the combined, filtered subset is 24.
        # The optimum power for the combined, unfiltered dataset is 14.
        # The optimum power for the yellow and brown subsets with DESeq2 filters is 30.
        # The optimum power for the unfiltered brown subset is 22.
        # The optimum power for the unfiltered yellow subset is 18.

#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================

datExprbean[] <- lapply(datExprbean, as.numeric)

# Cut height 0.95, 6 modules
netbean2 = blockwiseModules(datExprbean, power = 14, 
                           TOMType = "signed", minModuleSize = 30,
                           networkType = "signed",
                           detectCutHeight = 0.95, #Module tree cutoff
                           reassignThreshold = 0, #Gene significance cutoff
                           mergeCutHeight = 0.25,
                           numericLabels = TRUE, pamRespectsDendro = FALSE,
                           saveTOMs = TRUE,
                           saveTOMFileBase = "beanTOM", 
                           verbose = 3)

 #=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
beanmergedColors = labels2colors(netbean2$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(netbean2$dendrograms[[1]], beanmergedColors[netbean2$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


# Cutoff 0.95
beanmoduleLabels = netbean2$colors
View(beanmoduleLabels)
beanmoduleColors = labels2colors(netbean2$colors)
View(beanmoduleColors)
beanMEs = netbean2$MEs;
beangeneTree = netbean2$dendrograms[[1]];
save(beanMEs, beanmoduleLabels, beanmoduleColors, beangeneTree, 
     file = "Bean-02.95-signed-networkConstruction-auto.RData")
     #file = "Bean-02.95-signed-networkConstruction-auto-yellow.RData")
     #file = "Bean-02.95-signed-networkConstruction-auto-brown.RData")
     

