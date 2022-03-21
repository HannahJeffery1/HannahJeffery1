#************************************************************************;
#*Make QTL files;
#************************************************************************;
#remove(annotated_sigGenes_dev_2)
#remove(annotated_sigGenes_soak_4)

# DESeq2 files are saved directly in the "subset" folders.
#library(readxl)
#annotated_sigGenes_soak_4<-read_xlsx("annotated_sigGenes_dev_2.xlsx")
#annotated_sigGenes_soak_4<-read_xlsx("annotated_sigGenes_soak_4.xlsx")
annotated_sigGenes_soak_4<-read.csv(file = "WGCNA_norm10.6_soak.csv") #Found in Work_dir

Annotation<-read.table(file = "Pvulgaris_442_v2.1.annotation_info_dedup.txt") #In DESeq2 folder
Annotation<-Annotation[-c(1)]
names(Annotation)<-c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12")
Annotation<-Annotation[-c(2:3)]
Annot<-merge(Annotation, annotated_sigGenes_soak_4)

AnnotMod<-Annot


#Upload brown bean module data
Btur<-read.csv(file = "beangeneInfo_CookT_B.csv")
keep<-c("turquoise")
Btur<-Btur[Btur$beanmoduleColor %in% keep, ]
Btur<-Btur[-c(1:2)]
Btur<-Btur[-c(2:54)]
names(Btur)<-c("V2")
turB<-merge(Annot, Btur)

Bbro<-read.csv(file = "beangeneInfo_CookT_B.csv")
keep<-c("brown")
Bbro<-Bbro[Bbro$beanmoduleColor %in% keep, ]
Bbro<-Bbro[-c(1:2)]
Bbro<-Bbro[-c(2:54)]
names(Bbro)<-c("V2")
broB<-merge(Annot, Bbro)

Bmidblu<-read.csv(file = "beangeneInfo_CookT_B.csv")
keep<-c("midnightblue")
Bmidblu<-Bmidblu[Bmidblu$beanmoduleColor %in% keep, ]
Bmidblu<-Bmidblu[-c(1:2)]
Bmidblu<-Bmidblu[-c(2:54)]
names(Bmidblu)<-c("V2")
midbluB<-merge(Annot, Bmidblu)

Bred<-read.csv(file = "beangeneInfo_CookT_B.csv")
keep<-c("red")
Bred<-Bred[Bred$beanmoduleColor %in% keep, ]
Bred<-Bred[-c(1:2)]
Bred<-Bred[-c(2:54)]
names(Bred)<-c("V2")
redB<-merge(Annot, Bred)

Bpink<-read.csv(file = "beangeneInfo_CookT_B.csv")
keep<-c("pink")
Bpink<-Bpink[Bpink$beanmoduleColor %in% keep, ]
Bpink<-Bpink[-c(1:2)]
Bpink<-Bpink[-c(2:54)]
names(Bpink)<-c("V2")
pinkB<-merge(Annot, Bpink)


#Upload yellow bean module data
Ytur<-read.csv(file = "beangeneInfo_CookT_Y.csv")
keep<-c("turquoise")
Ytur<-Ytur[Ytur$beanmoduleColor %in% keep, ]
Ytur<-Ytur[-c(1:2)]
Ytur<-Ytur[-c(2:68)]
names(Ytur)<-c("V2")
turY<-merge(Annot, Ytur)

Ybro<-read.csv(file = "beangeneInfo_CookT_Y.csv")
keep<-c("brown")
Ybro<-Ybro[Ybro$beanmoduleColor %in% keep, ]
Ybro<-Ybro[-c(1:2)]
Ybro<-Ybro[-c(2:68)]
names(Ybro)<-c("V2")
broY<-merge(Annot, Ybro)

Yblu<-read.csv(file = "beangeneInfo_CookT_Y.csv")
keep<-c("blue")
Yblu<-Yblu[Yblu$beanmoduleColor %in% keep, ]
Yblu<-Yblu[-c(1:2)]
Yblu<-Yblu[-c(2:68)]
names(Yblu)<-c("V2")
bluY<-merge(Annot, Yblu)

Yred<-read.csv(file = "beangeneInfo_CookT_Y.csv")
keep<-c("red")
Yred<-Yred[Yred$beanmoduleColor %in% keep, ]
Yred<-Yred[-c(1:2)]
Yred<-Yred[-c(2:68)]
names(Yred)<-c("V2")
redY<-merge(Annot, Yred)

Ywhi<-read.csv(file = "beangeneInfo_CookT_Y.csv")
keep<-c("white")
Ywhi<-Ywhi[Ywhi$beanmoduleColor %in% keep, ]
Ywhi<-Ywhi[-c(1:2)]
Ywhi<-Ywhi[-c(2:68)]
names(Ywhi)<-c("V2")
whiY<-merge(Annot, Ywhi)

write.csv(turB, file="CT-Turquoise-brown-genes.csv")
write.csv(broB, file="CT-Brown-brown-genes.csv")
write.csv(redB, file="CT-Red-brown-genes.csv")
write.csv(pinkB, file="CT-Pink-brown-genes.csv")
write.csv(midbluB, file="CT-Midnightblue-brown-genes.csv")
write.csv(turY, file="CT-Turquoise-yellow-genes.csv")
write.csv(broY, file="CT-Brown-yellow-genes.csv")
write.csv(redY, file="CT-Red-yellow-genes.csv")
write.csv(bluY, file="CT-Blue-yellow-genes.csv")
write.csv(whiY, file="CT-White-yellow-genes.csv")



#Upload QTL files from YYRILCooktime
CG1.genes<-read.table(file = "CG2-1-1.txt", header = TRUE)
CG3.genes<-read.table(file = "CG2-2-1.txt", header = TRUE)
CG4.genes<-read.table(file = "CG2-2-2.txt", header = TRUE)
CG5.genes<-read.table(file = "CG5-1-1.txt", header = TRUE)
CG7.genes<-read.table(file = "CG8-1-1.txt", header = TRUE)
CG8.genes<-read.table(file = "CG8-1-2.txt", header = TRUE)
CG9.genes<-read.table(file = "CG8-2-1.txt", header = TRUE)
CG10.genes<-read.table(file = "CG8-2-2.txt", header = TRUE)
CG11.genes<-read.table(file = "CG8-2-3.txt", header = TRUE)
CG12.genes<-read.table(file = "CG10-1-1.txt", header = TRUE)
CG13.genes<-read.table(file = "CG10-1-2.txt", header = TRUE)

#Upload QTL files from YYRILwuptake
WU1<-read.table(file = "WU4-1.txt", header = TRUE)
WU2<-read.table(file = "WU4-2.txt", header = TRUE)
WU3<-read.table(file = "WU9-1-1.txt", header = TRUE)
WU4<-read.table(file = "WU9-1-2.txt", header = TRUE)
WU5<-read.table(file = "WU9-1-3.txt", header = TRUE)
WU6<-read.table(file = "WU10-1-1.txt", header = TRUE)
WU8<-read.table(file = "WU10-1-3.txt", header = TRUE)



##############
#Read in list of WGCNA/WGCNA+DESeq2 genes
turY<-read.csv("CT_Turquoise_yellow_genes_pvalue.csv")
broY<-read.csv("CT_Blue_yellow_genes_pvalue.csv")
redY<-read.csv("CT_Brown_yellow_genes_pvalue.csv")
bluY<-read.csv("CT_Red_yellow_genes_pvalue.csv")
whiY<-read.csv("CT_White_yellow_genes_pvalue.csv")

turY<-read.csv("Filtered_TurY_pvalues.csv")
broY<-read.csv("Filtered_BroY_pvalues.csv")
redY<-read.csv("Filtered_RedY_pvalues.csv")
bluY<-read.csv("Filtered_BluY_pvalues.csv")
whiY<-read.csv("Filtered_WhiY_pvalues.csv")
## ## ## ## ##
#How many QTL genes are in each module?
#Remember to upload new "filtered" files in EACH folder
#AnnotMod<-turB
#AnnotMod<-broB
#AnnotMod<-redB
#AnnotMod<-pinkB
#AnnotMod<-mbluB

AnnotMod<-turY
AnnotMod<-broY
AnnotMod<-redY
AnnotMod<-bluY
AnnotMod<-whiY


AnnotMod<-AnnotMod[-c(1)]
AnnotMod<-AnnotMod[-c(2)]
names(AnnotMod)<-c("V1")


# Cooking time QTLs
names(CG1.genes)<-c("V1")
CG1.genes_new<-merge(CG1.genes, AnnotMod)
##View(CG1.genes)
write.csv(CG1.genes_new, "CG1.genes.csv")

names(CG3.genes)<-c("V1")
CG3.genes_new<-merge(CG3.genes, AnnotMod)
#View(CG3.genes)
write.csv(CG3.genes_new, "CG3.genes.csv")

names(CG4.genes)<-c("V1")
CG4.genes_new<-merge(CG4.genes, AnnotMod)
#View(CG4.genes)
write.csv(CG4.genes_new, "CG4.genes.csv")

names(CG5.genes)<-c("V1")
CG5.genes_new<-merge(CG5.genes, AnnotMod)
#View(CG5.genes)
write.csv(CG5.genes_new, "CG5.genes.csv")

names(CG7.genes)<-c("V1")
CG7.genes_new<-merge(CG7.genes, AnnotMod)
#View(CG7.genes)
write.csv(CG7.genes_new, "CG7.genes.csv")

names(CG8.genes)<-c("V1")
CG8.genes_new<-merge(CG8.genes, AnnotMod)
#View(CG8.genes)
write.csv(CG8.genes_new, "CG8.genes.csv")

names(CG9.genes)<-c("V1")
CG9.genes_new<-merge(CG9.genes, AnnotMod)
#View(CG9.genes)
write.csv(CG9.genes_new, "CG9.genes.csv")

names(CG10.genes)<-c("V1")
CG10.genes_new<-merge(CG10.genes, AnnotMod)
#View(CG10.genes)
write.csv(CG10.genes_new, "CG10.genes.csv")

names(CG11.genes)<-c("V1")
CG11.genes_new<-merge(CG11.genes, AnnotMod)
#View(CG11.genes)
write.csv(CG11.genes_new, "CG11.genes.csv")

names(CG12.genes)<-c("V1")
CG12.genes_new<-merge(CG12.genes, AnnotMod)
#View(CG12.genes)
write.csv(CG12.genes_new, "CG12.genes.csv")

names(CG13.genes)<-c("V1")
CG13.genes_new<-merge(CG13.genes, AnnotMod)
#View(CG13.genes)
write.csv(CG13.genes_new, "CG13.genes.csv")


# Water uptake QTLs
names(WU1)<-c("V1")
WU1_new<-merge(WU1, AnnotMod)
#View(WU1)
write.csv(WU1_new, "WU1.csv")

names(WU2)<-c("V1")
WU2_new<-merge(WU2, AnnotMod)
#View(WU1)
write.csv(WU2_new, "WU2.csv")

names(WU3)<-c("V1")
WU3_new<-merge(WU3, AnnotMod)
#View(WU13)
write.csv(WU3_new, "WU3.csv")

names(WU4)<-c("V1")
WU4_new<-merge(WU4, AnnotMod)
#View(WU4)
write.csv(WU4_new, "WU4.csv")

names(WU5)<-c("V1")
WU5_new<-merge(WU5, AnnotMod)
#View(WU5)
write.csv(WU5_new, "WU5.csv")

names(WU6)<-c("V1")
WU6_new<-merge(WU6, AnnotMod)
#View(WU6)
write.csv(WU6_new, "WU6.csv")

names(WU8)<-c("V1")
WU8_new<-merge(WU8, AnnotMod)
#View(WU8)
write.csv(WU8_new, "WU8.csv")


#Combine large QTL datasets (example)
#Upload QTL files from YYRILCooktime
CG1.genes<-read.csv(file = "CG1.genes.csv", header = TRUE)
CG3.genes<-read.csv(file = "CG3.genes.csv", header = TRUE)
CG4.genes<-read.csv(file = "CG4.genes.csv", header = TRUE)
CG5.genes<-read.csv(file = "CG5.genes.csv", header = TRUE)
CG7.genes<-read.csv(file = "CG7.genes.csv", header = TRUE)
CG8.genes<-read.csv(file = "CG8.genes.csv", header = TRUE)
CG9.genes<-read.csv(file = "CG9.genes.csv", header = TRUE)
CG10.genes<-read.csv(file = "CG10.genes.csv", header = TRUE)
CG11.genes<-read.csv(file = "CG11.genes.csv", header = TRUE)
CG12.genes<-read.csv(file = "CG12.genes.csv", header = TRUE)
CG13.genes<-read.csv(file = "CG13.genes.csv", header = TRUE)

CG1.genes<-CG1.genes[-c(1)]
CG3.genes<-CG3.genes[-c(1)]
CG4.genes<-CG4.genes[-c(1)]
CG5.genes<-CG5.genes[-c(1)]
CG7.genes<-CG7.genes[-c(1)]
CG8.genes<-CG8.genes[-c(1)]
CG9.genes<-CG9.genes[-c(1)]
CG10.genes<-CG10.genes[-c(1)]
CG11.genes<-CG11.genes[-c(1)]
CG12.genes<-CG12.genes[-c(1)]
CG13.genes<-CG13.genes[-c(1)]

CT_QTL_genes<-rbind(CG1.genes, CG3.genes, 
                    CG4.genes, CG5.genes, CG7.genes, 
                    CG8.genes, CG9.genes, CG11.genes,
                    CG10.genes, CG12.genes, CG13.genes)
write.csv(CT_QTL_genes, file = "CT_QTL_genes.csv")

#Upload QTL files from YYRILwuptake
WU1<-read.csv(file = "WU1.csv", header = TRUE)
WU2<-read.csv(file = "WU2.csv", header = TRUE)
WU3<-read.csv(file = "WU3.csv", header = TRUE)
WU4<-read.csv(file = "WU4.csv", header = TRUE)
WU5<-read.csv(file = "WU5.csv", header = TRUE)
WU6<-read.csv(file = "WU6.csv", header = TRUE)
WU8<-read.csv(file = "WU8.csv", header = TRUE)

WU1<-WU1[-c(1)]
WU2<-WU2[-c(1)]
WU3<-WU3[-c(1)]
WU4<-WU4[-c(1)]
WU5<-WU5[-c(1)]
WU6<-WU6[-c(1)]
WU8<-WU8[-c(1)]

WU_QTL_genes<-rbind(WU1, WU2, WU3, WU4, WU5, WU6,
WU8)
write.csv(WU_QTL_genes, file = "WU_QTL_genes.csv")

