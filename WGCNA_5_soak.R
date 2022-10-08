#Credit: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/

#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================

rm(list=ls())

# Cutoff 0.95

# Display the current working directory
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Load the expression and trait data saved in the first part
#beanlnames = load(file = "Bean_cooking_data_WGCNA_10.6.RData");
beanlnames = load(file = "Bean_cooking_data_WGCNA_10.6_yellow.RData");
#beanlnames = load(file = "Bean_cooking_data_WGCNA_10.6_brown.RData");
#The variable beanlnames contains the names of loaded variables.
beanlnames
# Load network data saved in the second part.
#beanlnames = load(file = "Bean-02.95-signed-networkConstruction-auto.RData");
beanlnames = load(file = "Bean-02.95-signed-networkConstruction-auto-yellow.RData");
#beanlnames = load(file = "Bean-02.95-signed-networkConstruction-auto-brown.RData");
beanlnames

#=====================================================================================
#
#  Code chunk 2
#
#
#=====================================================================================


# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
datExprbean[]<-lapply(datExprbean, as.numeric)
View(datExprbean)
beandissTOM = 1-TOMsimilarityFromExpr(datExprbean, networkType = "signed", power = 14);
beannSelect = 704
# For reproducibility, we set the random seed
set.seed(57);

library(WGCNA)
beannGenes = ncol(datExprbean);
beannSamples = nrow(datExprbean);
# Recalculate MEs with color labels
beanMEs0 = moduleEigengenes(datExprbean, beanmoduleColors)$eigengenes
beanMEs = orderMEs(beanMEs0)
beanmoduleTraitCor = cor(beanMEs, beandatTraits, use = "p");
beanmoduleTraitPvalue = corPvalueStudent(beanmoduleTraitCor, beannSamples);

select = sample(beannGenes, size = beannSelect);
selectTOM = beandissTOM[select, select];
# There's no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = beanmoduleColors[select];
# Open a graphical window
sizeGrWindow(9,9)
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing 
# the color palette; setting the diagonal to NA also improves the clarity of the plot
beanplotDiss = selectTOM^7;
diag(beanplotDiss) = NA;
TOMplot(beanplotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")

#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================

#CookT_Y
# Recalculate module eigengenes
beanMEs = moduleEigengenes(datExprbean, beanmoduleColors)$eigengenes
# Isolate weight from the clinical traits
CookT_Y = as.data.frame(beandatTraits$CookT_Y);
names(CookT_Y) = "CookT_Y"
# Add the weight to existing module eigengenes
beanMET = orderMEs(cbind(beanMEs, CookT_Y))
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5);
par(cex = 0.9)
plotEigengeneNetworks(beanMET, "", 
                      marDendro = c(1,3.5,1,3.5), 
                      marHeatmap = c(3,4,1,1), 
                      cex.lab = 0.6, 
                      xLabelsAngle = 10)

#CookT_B
# Recalculate module eigengenes
beanMEs = moduleEigengenes(datExprbean, beanmoduleColors)$eigengenes
# Isolate weight from the clinical traits
CookT_B = as.data.frame(beandatTraits$CookT_B);
names(CookT_B) = "CookT_B"
# Add the weight to existing module eigengenes
beanMET = orderMEs(cbind(beanMEs, CookT_B))
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5);
par(cex = 0.9)
plotEigengeneNetworks(beanMET, "", 
                      marDendro = c(1,3.5,1,3.5), 
                      marHeatmap = c(3,4,1,1), 
                      cex.lab = 0.6, 
                      xLabelsAngle = 10)

#SoakT_Y
# Recalculate module eigengenes
beanMEs = moduleEigengenes(datExprbean, beanmoduleColors)$eigengenes
# Isolate weight from the clinical traits
SoakT_Y = as.data.frame(beandatTraits$SoakT_Y);
names(SoakT_Y) = "SoakT_Y"
# Add the weight to existing module eigengenes
beanMET = orderMEs(cbind(beanMEs, SoakT_Y))
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(9,9);
par(cex = 0.9)
plotEigengeneNetworks(beanMET, "", 
                      marDendro = c(1,3.5,1,3.5), 
                      marHeatmap = c(3,4,1,1), 
                      cex.lab = 0.6, 
                      xLabelsAngle = 10)

#SoakT_B
# Recalculate module eigengenes
beanMEs = moduleEigengenes(datExprbean, beanmoduleColors)$eigengenes
# Isolate weight from the clinical traits
SoakT_B = as.data.frame(beandatTraits$SoakT_B);
names(SoakT_B) = "SoakT_B"
# Add the weight to existing module eigengenes
beanMET = orderMEs(cbind(beanMEs, SoakT_B))
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5);
par(cex = 0.9)
plotEigengeneNetworks(beanMET, "", 
                      marDendro = c(1,4,1,3.5), 
                      marHeatmap = c(3,4,1,1), 
                      cex.lab = 0.6, 
                      xLabelsAngle = 10)

#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================


# Plot the dendrogram
sizeGrWindow(6,6);
par(cex = 1.0)
plotEigengeneNetworks(beanMET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
plotEigengeneNetworks(beanMET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 360)

