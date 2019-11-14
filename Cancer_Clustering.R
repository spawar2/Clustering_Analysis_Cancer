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
Colon2 <- read.csv("Colon2.csv", header=TRUE, sep=",")
Lung1 <- read.csv("Lung1.csv", header=TRUE, sep=",")
Lung2 <- read.csv("Lung2.csv", header=TRUE, sep=",")
MM <- read.csv("MM.csv", header=TRUE, sep=",")
Oesophageal <- read.csv("Oesophageal.csv", header=TRUE, sep=",")

names(Breast)[2:287]<-"B"
names(Colon1)[2:586]<-"C"
names(Colon2)[2:291]<-"C"
names(Lung1)[2:308]<-"L"
names(Lung2)[2:463]<-"L"
names(MM)[2:489]<-"M"
names(Oesophageal)[2:9]<-"O"

names(Colon1)[1] <- "ID_REF"
names(Colon2)[1] <- "ID_REF"
names(Lung1)[1] <- "ID_REF"
names(Lung2)[1] <- "ID_REF"
names(MM)[1] <- "ID_REF"
names(Oesophageal)[1] <- "ID_REF"

###MM Skipped because the probe names are different###

#Merge datasets to form one object with all the chips

one <- merge(Colon1, Breast, by="ID_REF", all=TRUE)
two <- merge(one, Colon2, by="ID_REF", all=TRUE)
three <- merge(two, Lung1, by="ID_REF", all=TRUE)
four <- merge(three, Lung2, by="ID_REF", all=TRUE)
five <- merge(four, Oesophageal, by="ID_REF", all=TRUE)

dataclusterno2 <- as.data.frame(lapply(five, as.numeric))
dataclusterno2[is.na(dataclusterno2)] <- 0

# Call the number_clusters function
## # data.exp is the original expression matrix object ouputted from
## # the input_file function
## # User chooses the gap_statistic option by making gap_statistic equal TRUE
## # The Fixed argument is also set to NULL
# OBTAIN CLUSTER NUMBER FROM BELOW

cluster_num <- number_clusters(data.exp=dataclusterno2, Fixed=NULL,
     gap_statistic=TRUE)









