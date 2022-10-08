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
#beanlnames = load(file = "Bean_cooking_data_WGCNA_10.6_yellow.RData");
beanlnames = load(file = "Bean_cooking_data_WGCNA_10.6_brown.RData");
#The variable beanlnames contains the names of loaded variables.
beanlnames
# Load network data saved in the second part.
#beanlnames = load(file = "Bean-02.95-signed-networkConstruction-auto.RData");
#beanlnames = load(file = "Bean-02.95-signed-networkConstruction-auto-yellow.RData");
beanlnames = load(file = "Bean-02.95-signed-networkConstruction-auto-brown.RData");
beanlnames
View(beangeneTree)

#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================

View(datExprbean)
# Define numbers of genes and samples
library(WGCNA)
beannGenes = ncol(datExprbean);
beannSamples = nrow(datExprbean);
# Recalculate MEs with color labels
beanMEs0 = moduleEigengenes(datExprbean, beanmoduleColors)$eigengenes
beanMEs = orderMEs(beanMEs0)
beanmoduleTraitCor = cor(beanMEs, beandatTraits, use = "p");
beanmoduleTraitPvalue = corPvalueStudent(beanmoduleTraitCor, beannSamples);

#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================

sizeGrWindow(10,6)
# Will display correlations and their p-values
beantextMatrix =  paste(signif(beanmoduleTraitCor, 2), "\n(",
                    signif(beanmoduleTraitPvalue, 1), ")", sep = "");
dim(beantextMatrix) = dim(beanmoduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = beanmoduleTraitCor,
               xLabels = names(beandatTraits),
               yLabels = names(beanMEs),
               ySymbols = names(beanMEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = beantextMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-Trait Relationships"))


#=====================================================================================
#
#  Code chunk 4
#
#=====================================================================================

# Define variable containing the trait using the datTraitbean data set
# Do each trait one-at-a-time
#CookT_Y
View(beandatTraits)
binary_data_1 = as.data.frame(beandatTraits$Soak_time_hr);
names(binary_data_1) = "binary_data_1"
# names (colors) of the modules
beanmodNames = substring(names(beanMEs), 3)

beangeneModuleMembership = as.data.frame(cor(datExprbean, beanMEs, use = "p"));
beanMMPvalue = as.data.frame(corPvalueStudent(as.matrix(beangeneModuleMembership), beannSamples));

names(beangeneModuleMembership) = paste("beanMM", beanmodNames, sep="");
names(beanMMPvalue) = paste("beanp.MM", beanmodNames, sep="");

beangeneTraitSignificance = as.data.frame(cor(datExprbean, binary_data_1, use = "p"));
beanGSPvalue = as.data.frame(corPvalueStudent(as.matrix(beangeneTraitSignificance), beannSamples));

names(beangeneTraitSignificance) = paste("beanGS.", names(binary_data_1), sep="");
names(beanGSPvalue) = paste("beanp.GS.", names(binary_data_1), sep="");

#############################################################

#CookT_B
View(beandatTraits)
binary_data_2 = as.data.frame(beandatTraits$Cook_time_min);
names(binary_data_2) = "binary_data_2"
# names (colors) of the modules
beanmodNames = substring(names(beanMEs), 3)

beangeneModuleMembership = as.data.frame(cor(datExprbean, beanMEs, use = "p"));
beanMMPvalue = as.data.frame(corPvalueStudent(as.matrix(beangeneModuleMembership), beannSamples));

names(beangeneModuleMembership) = paste("beanMM", beanmodNames, sep="");
names(beanMMPvalue) = paste("beanp.MM", beanmodNames, sep="");

beangeneTraitSignificance = as.data.frame(cor(datExprbean, binary_data_2, use = "p"));
beanGSPvalue = as.data.frame(corPvalueStudent(as.matrix(beangeneTraitSignificance), beannSamples));

names(beangeneTraitSignificance) = paste("beanGS.", names(binary_data_2), sep="");
names(beanGSPvalue) = paste("beanp.GS.", names(binary_data_2), sep="");

#############################################################

#SoakT_Y
View(beandatTraits)
binary_data_3 = as.data.frame(beandatTraits$CookT_B);
names(binary_data_3) = "binary_data_3"
# names (colors) of the modules
beanmodNames = substring(names(beanMEs), 3)

beangeneModuleMembership = as.data.frame(cor(datExprbean, beanMEs, use = "p"));
beanMMPvalue = as.data.frame(corPvalueStudent(as.matrix(beangeneModuleMembership), beannSamples));

names(beangeneModuleMembership) = paste("beanMM", beanmodNames, sep="");
names(beanMMPvalue) = paste("beanp.MM", beanmodNames, sep="");

beangeneTraitSignificance = as.data.frame(cor(datExprbean, binary_data_3, use = "p"));
beanGSPvalue = as.data.frame(corPvalueStudent(as.matrix(beangeneTraitSignificance), beannSamples));

names(beangeneTraitSignificance) = paste("beanGS.", names(binary_data_3), sep="");
names(beanGSPvalue) = paste("beanp.GS.", names(binary_data_3), sep="");

#############################################################

#SoakT_B
View(beandatTraits)
binary_data_4 = as.data.frame(beandatTraits$SoakT_B);
names(binary_data_4) = "binary_data_4"
# names (colors) of the modules
beanmodNames = substring(names(beanMEs), 3)

beangeneModuleMembership = as.data.frame(cor(datExprbean, beanMEs, use = "p"));
beanMMPvalue = as.data.frame(corPvalueStudent(as.matrix(beangeneModuleMembership), beannSamples));

names(beangeneModuleMembership) = paste("beanMM", beanmodNames, sep="");
names(beanMMPvalue) = paste("beanp.MM", beanmodNames, sep="");

beangeneTraitSignificance = as.data.frame(cor(datExprbean, binary_data_4, use = "p"));
beanGSPvalue = as.data.frame(corPvalueStudent(as.matrix(beangeneTraitSignificance), beannSamples));

names(beangeneTraitSignificance) = paste("beanGS.", names(binary_data_4), sep="");
names(beanGSPvalue) = paste("beanp.GS.", names(binary_data_4), sep="");

#=====================================================================================
#
#  Code chunk 5
#
#=====================================================================================

module = "pink"
column = match(module, beanmodNames);
beanmoduleGenes = beanmoduleColors==module;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(beangeneModuleMembership[beanmoduleGenes, column]),
                   abs(beangeneTraitSignificance[beanmoduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene Significance for Brown Bean Cook Time",
                   main = paste("Module Membership vs. Gene Significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "black")

#########################################

module = "pink"
column = match(module, beanmodNames);
beanmoduleGenes = beanmoduleColors==module;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(beangeneModuleMembership[beanmoduleGenes, column]),
                   abs(beangeneTraitSignificance[beanmoduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene Significance for Brown Bean Soak Time",
                   main = paste("Module Membership vs. Gene Significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "black")

#=====================================================================================
#
#  Code chunk 6
#
#=====================================================================================

names(datExprbean)

#=====================================================================================
#
#  Code chunk 7
#
#====================================================================================

#  These identifiers will depend on the module color and trait you are working with
names(datExprbean)[beanmoduleColors=="turquoise"]

#=====================================================================================
#
#  Code chunk 8
#
#=====================================================================================

# I made WGCNA_file by hand by averaging all of the gene expression values for the replicates.  
# Add V1, V2, V3, etc. to an extra column in WGCNA_file.  Rename it WGCNA_file_annot.
# Merge WGCNA_file_annot with "beanmoduleLabels".
#beanannot = read.csv("soakDEGs_2.csv");
beanannot = read.csv("WGCNA_norm10.6_soak_yellow.csv");
#beanannot = read.csv("WGCNA_norm10.6_soak_brown.csv");
beanprobes = as.data.frame(names(datExprbean))
#View(beanprobes$`names(datExprbean)`)
#View(beanannot)
beanprobes2annot = match(beanprobes$`names(datExprbean)`, beanannot$V3)
# The following is the number or probes without annotation:
sum(is.na(beanprobes2annot))
# Should return 0.
# R excludes duplicates.

#=====================================================================================
#
#  Code chunk 9
#
#=====================================================================================

# Create the starting data frame (one trait at-a-time)
beangeneInfo = data.frame(annotation = beanprobes,
                       Annotation  = beanannot$V2[beanprobes2annot],
                       beanmoduleColor = beanmoduleColors,
                       beangeneTraitSignificance,
                       beanGSPvalue)
#View(beangeneInfo)
#View(beanMEs)
# Order modules by their significance (one-at-a-time)
# CookT_Y
beanmodOrder = order(-abs(cor(beanMEs, beandatTraits$CookT_Y, use = "p")));
# CookT_B
beanmodOrder = order(-abs(cor(beanMEs, beandatTraits$CookT_B, use = "p")));
# SoakT_Y
beanmodOrder = order(-abs(cor(beanMEs, beandatTraits$SoakT_Y, use = "p")));
# SoakT_B
beanmodOrder = order(-abs(cor(beanMEs, beandatTraits$SoakT_B, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(beangeneModuleMembership))
{
  beanoldNames = names(beangeneInfo)
  beangeneInfo = data.frame(beangeneInfo, beangeneModuleMembership[, beanmodOrder[mod]], 
                         beanMMPvalue[, beanmodOrder[mod]]);
  names(beangeneInfo) = c(beanoldNames, paste("MM.", beanmodNames[beanmodOrder[mod]], sep=""),
                       paste("p.MM.", beanmodNames[beanmodOrder[mod]], sep=""))
}

# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
#CookT_Y
beangeneOrder = order(beangeneInfo$beanmoduleColor, 
                      -abs(beangeneInfo$beanGS.binary_data_1));
beangeneInfo = beangeneInfo[beangeneOrder, ];

#CookT_B
beangeneOrder = order(beangeneInfo$beanmoduleColor, 
                      -abs(beangeneInfo$beanGS.binary_data_2));
beangeneInfo = beangeneInfo[beangeneOrder, ];

#SoakT_Y
beangeneOrder = order(beangeneInfo$beanmoduleColor, 
                      -abs(beangeneInfo$beanGS.binary_data_3));
beangeneInfo = beangeneInfo[beangeneOrder, ];

#SoakT_B
beangeneOrder = order(beangeneInfo$beanmoduleColor, 
                      -abs(beangeneInfo$beanGS.binary_data_4));
beangeneInfo = beangeneInfo[beangeneOrder, ]
View(beangeneInfo)

#=====================================================================================
#
#  Code chunk 10
#
#=====================================================================================

# CookT_Y
write.csv(beangeneInfo, file = "beangeneInfo_CookT_Y.csv")
# CookT_B
write.csv(beangeneInfo, file = "beangeneInfo_CookT_B.csv")
# SoakT_Y
write.csv(beangeneInfo, file = "beangeneInfo_SoakT_Y.csv")
# SoakT_B
write.csv(beangeneInfo, file = "beangeneInfo_SoakT_B.csv")
