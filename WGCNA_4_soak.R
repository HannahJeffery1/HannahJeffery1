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
beanlnames
# Load network data saved in the second part.
beanlnames = load(file = "Bean-02.95-signed-networkConstruction-auto.RData");
#beanlnames = load(file = "Bean-02.95-signed-networkConstruction-auto-yellow.RData");
#beanlnames = load(file = "Bean-02.95-signed-networkConstruction-auto-brown.RData");
beanlnames

#=====================================================================================
#
#  Code chunk 2
# Output gene lists for use with online software and services.
#
#=====================================================================================

#beanannot = read.csv("soakDEGs_2.csv");
beanannot = read.csv("WGCNA_norm10.6_soak.csv");
#beanannot = read.csv("WGCNA_norm10.6_soak_yellow.csv");
#beanannot = read.csv("WGCNA_norm10.6_soak_brown.csv");
beanprobes = as.data.frame(names(datExprbean))
View(beanprobes$`names(datExprbean)`)
beanprobes2annot = match(beanprobes$`names(datExprbean)`, beanannot$V3)
# The following is the number or probes without annotation:
sum(is.na(beanprobes2annot))
# Should return 0.
# R excludes duplicates.
# Get the corresponding Locus Link IDs
beanallLLIDs = beanannot$V3[beanprobes2annot];

# $ Choose interesting modules
###############################################################################################
#turquoise
beanintModules = c("turquoise")

for (module in beanintModules)
{
  # Select module probes
  beanmodGenes = (beanmoduleColors==module)
  # Get their ID codes
  beanmodLLIDs = beanallLLIDs[beanmodGenes];
  # Write them into a file
  fileName = paste("transcript-", module, ".txt", sep="");
  turquoise<-as.data.frame(beanmodLLIDs);
  turquoise$color<-c("turquoise")
}

View(turquoise)

###############################################################################################

# As background in the enrichment analysis, we will use all probes in the analysis.

fileName = paste("transcript.txt", sep="");
write.table(as.data.frame(beanallLLIDs), file = "all-0.95-turquoise",
            row.names = FALSE, col.names = FALSE)

# Now, create specific module files

fileName = paste("transcript-turquoise.txt", sep="");
write.table(as.data.frame(turquoise), file = "turquoise-0.95-turquoise",
            row.names = FALSE, col.names = FALSE)


