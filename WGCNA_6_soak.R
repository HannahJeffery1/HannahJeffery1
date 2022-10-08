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
beanlnames = load(file = "Bean_cooking_data_WGCNA_10.6.RData");
#beanlnames = load(file = "Bean_cooking_data_WGCNA_10.6_yellow.RData");
#beanlnames = load(file = "Bean_cooking_data_WGCNA_10.6_brown.RData");
#The variable beanlnames contains the names of loaded variables.
beanlnames
# Load network data saved in the second part.
beanlnames = load(file = "Bean-02.95-signed-networkConstruction-auto.RData");
#beanlnames = load(file = "Bean-02.95-signed-networkConstruction-auto-yellow.RData");
#beanlnames = load(file = "Bean-02.95-signed-networkConstruction-auto-brown.RData");
beanlnames

#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================


# Recalculate topological overlap
datExprbean[]<-lapply(datExprbean, as.numeric)
library(WGCNA)
beanTOM = TOMsimilarityFromExpr(datExprbean, networkType = "signed", power = 14);

# Read in the annotation file
#beanannot = read.csv(file = "soakDEGs_2.csv");
beanannot = read.csv("WGCNA_norm10.6_soak.csv");
#beanannot = read.csv("WGCNA_norm10.6_soak_yellow.csv");
#beanannot = read.csv("WGCNA_norm10.6_soak_brown.csv");
# Select a module
#module = c("greenyellow");
#module = c("black");
#module = c("blue");
#module = c("brown");
#module = c("grey60");
#module = c("magenta");
module = c();"turquoise"
# Select module probes
probes = names(datExprbean)
inModule = (beanmoduleColors==module);
modProbes = probes[inModule];
# Select the corresponding Topological Overlap
modTOM = beanTOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into an edge list file VisANT can read
#View(beanannot)
vis = exportNetworkToVisANT(modTOM,
                            file = paste("VisANTInput-", module, ".txt", sep=""),
                            weighted = TRUE,
                            threshold = 0,
                            probeToGene = data.frame(beanannot$V3, 
                                                     beanannot$V2) )

#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================


nTop = 30;
IMConn = softConnectivity(datExprbean[, modProbes]);
top = (rank(-IMConn) <= nTop)
vis = exportNetworkToVisANT(modTOM[top, top],
                            file = paste("VisANTInput-", module, "-top30.txt", sep=""),
                            weighted = TRUE,
                            threshold = 0,
                            probeToGene = data.frame(beanannot$V3, 
                                                     beanannot$V2) )


#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================


# Select modules
#modules = c("greenyellow");
#modules = c("black");
#modules = c("blue");
#modules = c("brown");
#modules = c("grey60");
#modules = c("magenta");
modules = c("turquoise");
# Select module probes
probes = names(datExprbean)
inModule = is.finite(match(beanmoduleColors, modules));
modProbes = probes[inModule];
modGenes = beanannot$V2[match(modProbes, beanannot$V3)];
# Select the corresponding Topological Overlap
modTOM = beanTOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = beanmoduleColors[inModule]);

