#######################################
#######################################
#TopGO analysis
#######################################
#######################################
#Credit: https://bioconductor.org/packages/devel/bioc/vignettes/topGO/inst/doc/topGO.pdf

rm(list = ls())
#Set parameters
options(stringsAsFactors = FALSE)

#Install/update topGO base package
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

BiocManager::install("topGO")

#Data preparation
library(topGO)
library("biomaRt")


listMarts()
listMarts(host = "plants.ensembl.org")

#https://www.biostars.org/p/250927/
#Collect gene names from biomart
mart <- biomaRt::useMart(biomart = "plants_mart",
                         dataset = "pvulgaris_eg_gene",
                         host = 'plants.ensembl.org')

#Get ensembl gene ids and GO terms
GTOGO <- biomaRt::getBM(attributes = c( "ensembl_gene_id",
                                        "go_id"), mart = mart)

#Examine result
head(GTOGO)
##write.csv(GTOGO$ensembl_gene_id, file = "all_genes_orig")
#Upload manually edited GO file
#GOnames<-read.csv(file = "all_genes_orig")
#GTOGO$ensembl_gene_id<-GOnames

#Remove blank entries
GTOGO <- GTOGO[GTOGO$go_id != '',]

#Convert from table format to list format
geneID2GO <- by(GTOGO$go_id,
                GTOGO$ensembl_gene_id,
                function(x) as.character(x))

#Examine result
head(geneID2GO)

###############START HERE##################
###########################################

#Edit file (change Phvul. to PHAVU_ and add a "g" to the end of the name)
#int.genes<-read.csv(file = "TurquoiseB_pvalues.csv")

#convert<-read.csv(file = "beangeneInfo_SoakT_Y.csv", header = T)
#convert<-convert[-c(1)]
#convert<-convert[-c(3:55)]
#convert<-read.csv(file = "beangeneInfo_SoakT_Y.csv", header = T)
#convert<-convert[-c(1)]
#convert<-convert[-c(3:69)]

#int.genes<-merge(convert, int.genes)
#int.genes<-int.genes[-c(1)]
#write.csv(int.genes, file = "int.genes.csv")

###Finished files
int.genes<-read.csv(file = "LightgreenB_pvalues.csv")
#int.genes<-read.csv(file = "Individuals_DESeq_alltimes.csv")

#Remove duplicate genes (remnants of the GO file)
names(int.genes)<-c("x", "pvalue")
int.genes<-int.genes[!duplicated(int.genes$x), ]

View(int.genes)

#Create a list of all the genes
all_genes<-GTOGO$ensembl_gene_id
all_genes<-as.data.frame(all_genes)
all_genes<-all_genes[!duplicated(all_genes$all_genes), ]

#Create the parent dataset
int.genes_1<-merge(all_genes, int.genes)

#Create a list of gene IDs
col<-int.genes_1$x

#int.genes_2<-sort(unique(as.character(int.genes_1$pvalue)))
int.genes_2<-int.genes_1$pvalue
int.genes_2<-as.numeric(int.genes_2)
int.genes_2<-as.vector(int.genes_2)
names(int.genes_2)<-c("pvalue")
#https://stackoverflow.com/questions/12923809/retain-rownames-as-names-of-vector-in-r
names(int.genes_2)<-col
View(int.genes_2)

#https://rdrr.io/bioc/topGO/src/inst/scripts/GO_ALL.R
topDiffGenes <- function(allScore) {
  return(allScore < 0.01)
}

####Change the ontology (BP --> CC --> MF)
#https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/friday/enrichment.html
sampleGOdata <- new("topGOdata",
                    #ontology = "BP",
                    #ontology = "CC",
                    ontology = "MF",
                    allGenes = int.genes_2,
                    geneSel = topDiffGenes,
                    nodeSize = 10,
                    annot = annFUN.gene2GO,
                    gene2GO = geneID2GO)


#Fisher's exact test
resultFisher <- runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")

resultFisher

resultKS <- runTest(sampleGOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(sampleGOdata, algorithm = "elim", statistic = "ks")

#Analysis of Fisher's exact test results
allRes<-GenTable(sampleGOdata, classicFisher = resultFisher,
                 classicKS = resultKS, elimKS = resultKS.elim,
                 orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 10)

pValue.classic<-score(resultKS)
pValue.elim<-score(resultKS.elim)[names(pValue.classic)]
gstat<-termStat(sampleGOdata, names(pValue.classic))
gSize<-gstat$Annotated/max(gstat$Annotated)*4
#Create a map of the results
#https://rpubs.com/aemoore62/TopGo_colMap_Func_Troubleshoot
colMap <- function(x) {
  .col <- rep(rev(heat.colors(length(unique(x)))), time = table(x))
  return(.col[match(1:length(x), order(x))])
}
gCol<-colMap(gstat$Significant)
plot(pValue.classic, pValue.elim, xlab = "p-value classic",  ylab = "p-value elim", 
     main = "Stringent (elim) against Classic p-value Selection", pch = 19, cex = gSize, col = gCol)

#Find the number of GO terms with non-monotonic behavior (these are not likely to be significant at all)
sel.go<-names(pValue.classic)[pValue.elim<pValue.classic]
cbind(termStat(sampleGOdata, sel.go),
      elim = pValue.elim[sel.go],
      classic = pValue.classic[sel.go])

#Create graph of interconnected GO terms (save as a PDF)
showSigOfNodes(sampleGOdata, score(resultKS.elim), firstSigNodes = 5, useInfo = 'all')

#Create a table of the results
allRes<-GenTable(sampleGOdata, classic = resultFisher,
                 KS = resultKS, 
                 elimKS = resultKS.elim,
                 orderBy = "elimKS",
                 topNodes = 30)
View(allRes)

#write.csv(allRes, file = "Lightgreen_GO_BP.csv")
#write.csv(allRes, file = "Lightgreen_GO_CC.csv")
write.csv(allRes, file = "Lightgreen_GO_MF.csv")
