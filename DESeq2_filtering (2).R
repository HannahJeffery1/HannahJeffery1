############################
#DESeq2 filtering
############################

#Perform this step once:

TurB<-read.csv("CT_Turquoise_brown_genes_pvalue.csv")
TurB<-TurB[-c(1)]
TurB<-TurB[-c(2)]
names(TurB)<-c("V2")

BroB<-read.csv("CT_Brown_brown_genes_pvalue.csv")
BroB<-BroB[-c(1)]
BroB<-BroB[-c(2)]
names(BroB)<-c("V2")

RedB<-read.csv("CT_Red_brown_genes_pvalue.csv")
RedB<-RedB[-c(1)]
RedB<-RedB[-c(2)]
names(RedB)<-c("V2")

PinkB<-read.csv("CT_Pink_brown_genes_pvalue.csv")
PinkB<-PinkB[-c(1)]
PinkB<-PinkB[-c(2)]
names(PinkB)<-c("V2")

MbluB<-read.csv("CT_Midnightblue_brown_genes_pvalue.csv")
MbluB<-MbluB[-c(1)]
MbluB<-MbluB[-c(2)]
names(MbluB)<-c("V2")

########################

TurY<-read.csv("CT_Turquoise_yellow_genes_pvalue.csv")
TurY<-TurY[-c(1)]
TurY<-TurY[-c(2)]
names(TurY)<-c("V2")

BroY<-read.csv("CT_Brown_yellow_genes_pvalue.csv")
BroY<-BroY[-c(1)]
BroY<-BroY[-c(2)]
names(BroY)<-c("V2")

RedY<-read.csv("CT_Red_yellow_genes_pvalue.csv")
RedY<-RedY[-c(1)]
RedY<-RedY[-c(2)]
names(RedY)<-c("V2")

BluY<-read.csv("CT_Blue_yellow_genes_pvalue.csv")
BluY<-BluY[-c(1)]
BluY<-BluY[-c(2)]
names(BluY)<-c("V2")

WhiY<-read.csv("CT_White_yellow_genes_pvalue.csv")
WhiY<-WhiY[-c(1)]
WhiY<-WhiY[-c(2)]
names(WhiY)<-c("V2")

########################

#Load data:
Filter<-read.table(file = "sigGenes_1.txt", header = TRUE)
Filter<-Filter[-c(2:6)]

Filtered<-merge(TurB, Filter)
names(Filtered)<-c("Gene", "pvalue")
write.csv(Filtered, file = "Filtered_TurB_pvalues.csv")

Filtered<-merge(BroB, Filter)
names(Filtered)<-c("Gene", "pvalue")
write.csv(Filtered, file = "Filtered_BroB_pvalues.csv")

Filtered<-merge(RedB, Filter)
names(Filtered)<-c("Gene", "pvalue")
write.csv(Filtered, file = "Filtered_RedB_pvalues.csv")

Filtered<-merge(PinkB, Filter)
names(Filtered)<-c("Gene", "pvalue")
write.csv(Filtered, file = "Filtered_PinkB_pvalues.csv")

Filtered<-merge(MbluB, Filter)
names(Filtered)<-c("Gene", "pvalue")
write.csv(Filtered, file = "Filtered_MbluB_pvalues.csv")

Filtered<-merge(TurY, Filter)
names(Filtered)<-c("Gene", "pvalue")
write.csv(Filtered, file = "Filtered_TurY_pvalues.csv")

Filtered<-merge(BroY, Filter)
names(Filtered)<-c("Gene", "pvalue")
write.csv(Filtered, file = "Filtered_BroY_pvalues.csv")

Filtered<-merge(RedY, Filter)
names(Filtered)<-c("Gene", "pvalue")
write.csv(Filtered, file = "Filtered_RedY_pvalues.csv")

Filtered<-merge(BluY, Filter)
names(Filtered)<-c("Gene", "pvalue")
write.csv(Filtered, file = "Filtered_BluY_pvalues.csv")

Filtered<-merge(WhiY, Filter)
names(Filtered)<-c("Gene", "pvalue")
write.csv(Filtered, file = "Filtered_WhiY_pvalues.csv")








