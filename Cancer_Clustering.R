#11/13/19, Shrikant Pawar
#Clustering analysis for different cancer type.
library(GEOquery)
library(Biobase)
library(preprocessCore)
library(multiClust)
suppressPackageStartupMessages(library(ctc))
suppressPackageStartupMessages(library(gplots))
suppressPackageStartupMessages(library(dendextend))
library(ctc)
library(gplots)
library(dendextend)
library(graphics)
library(grDevices)
library(amap)
setwd("C:/Users/Bio-user/Desktop/IWBBIO 2020/CSV") 
Breast <- read.csv("Beast.csv", header=TRUE, sep=",")
Colon1 <- read.csv("Colon1.csv", header=TRUE, sep=",")
#Colon2 <- read.csv("Colon2.csv", header=TRUE, sep=",")
Lung1 <- read.csv("Lung1.csv", header=TRUE, sep=",")
#Lung2 <- read.csv("Lung2.csv", header=TRUE, sep=",")
MM <- read.csv("MM.csv", header=TRUE, sep=",")
Oesophageal <- read.csv("Oesophageal.csv", header=TRUE, sep=",")
Ovarian <- read.csv("Ovarian.csv", header=TRUE, sep=",")

names(Breast)[2:287]<-"B"
names(Colon1)[2:586]<-"C"
#names(Colon2)[2:291]<-"C"
names(Lung1)[2:308]<-"L"
#names(Lung2)[2:463]<-"L"
names(MM)[2:489]<-"M"
names(Oesophageal)[2:9]<-"O"
names(Ovarian)[2:9]<-"OV"

names(Colon1)[1] <- "ID_REF"
#names(Colon2)[1] <- "ID_REF"
names(Lung1)[1] <- "ID_REF"
#names(Lung2)[1] <- "ID_REF"
names(MM)[1] <- "ID_REF"
names(Oesophageal)[1] <- "ID_REF"
names(Ovarian)[1] <- "ID_REF"

###MM Skipped because the probe names are different###

#Merge datasets to form one object with all the chips

one <- merge(Colon1, Breast, by="ID_REF", all=TRUE)
two <- merge(one, Lung1, by="ID_REF", all=TRUE)
three <- merge(two, Oesophageal, by="ID_REF", all=TRUE)
four <- merge(three, Ovarian, by="ID_REF", all=TRUE)
#five <- merge(four, Oesophageal, by="ID_REF", all=TRUE)
#six <- merge(five, Ovarian, by="ID_REF", all=TRUE)

dataclusterno2 <- as.data.frame(lapply(four, as.numeric))
dataclusterno2[is.na(dataclusterno2)] <- 0

# Call the number_clusters function
## # data.exp is the original expression matrix object ouputted from
## # the input_file function
## # User chooses the gap_statistic option by making gap_statistic equal TRUE
## # The Fixed argument is also set to NULL
# OBTAIN CLUSTER NUMBER FROM BELOW

cluster_num <- number_clusters(data.exp=dataclusterno2,Fixed=NULL,
     gap_statistic=TRUE)

Clustering k = 1,2,..., K.max (= 8): .. done
Bootstrapping, b = 1,2,..., B (= 100)  [one "." per sample]:
.................................................. 50 
.................................................. 100 
[1] "The gap statistic cluster number is: 5"

# K means run

breast <- Breast[,2:287]
colon <- Colon1[,2:287]
lung <- Lung1[,2:308]
mm <- MM[,2:287]
oesoph <- Oesophageal[,2:9]
ovary <- Ovarian[,2:296]

kmeans_analysis <- cluster_analysis(sel.exp=na.omit(l),
    cluster_type="Kmeans", seed=5,
    distance=NULL, linkage_type=NULL, 
    gene_distance=NULL, num_clusters=100,
    data_name="LungClusters", probe_rank="SD_Rank",
    probe_num_selection="Fixed_Probe_Num",
    cluster_num_selection="Fixed_Clust_Num")


hclust_analysis <- cluster_analysis(sel.exp=breast[,1:20],
    cluster_type="HClust",
    distance="euclidean", linkage_type="ward.D2", 
    gene_distance="correlation",
    num_clusters=5, data_name="GSE2034 Breast", 
    probe_rank="SD_Rank", probe_num_selection="Fixed_Probe_Num",
    cluster_num_selection="Fixed_Clust_Num")
---------------------------------------------------------

clusters <- hclust(dist(t(ovary)))
plot(clusters, cex = 0.6, hang = -1)
sub_grp <- cutree(clusters, k = 5)
plot(clusters, cex = 0.6)
rect.hclust(clusters, k = 5, border = 2:5)

write.table(sub_grp, "ovaryid.txt", sep="\t")


----------------------Heatmap on genes accross cancers-----------------

All <- read.csv("Heatmap.csv", header=TRUE, sep=",")
rownames(All) <- All[,1]
All <- All[,2:7]
Alll <- data.matrix(All)
Alll[,6] <- as.numeric(as.character(All[,6]))
Alll[is.na(Alll)] <- 0

newdataup <- Alll[which(Alll[,1]>1.5
& Alll[,2]>1.5
& Alll[,3]>1.5
& Alll[,4]>1.5
& Alll[,5]>1.5
& Alll[,6]>1.5),]

newdatadown <- Alll[which(Alll[,1] < -1.5
& Alll[,2] < -1.5
& Alll[,3] < -1.5
& Alll[,4] < -1.5
& Alll[,5] < -1.5
& Alll[,6] < -1.5),]

total <- rbind(newdataup, newdatadown)

library("gplots")
heatmap.2(total, scale = "column", col = bluered(100), 
          trace = "none", density.info = "none")





