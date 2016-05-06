######################################################################

                        ##R code for task 13

######################################################################

## Combine tasks 1 to 12 into this document

######################################################################
## pre-analysis stuff
setwd('~/Documents/masters/data/task13/')

#source the file to get all the functions from task10.R
source('~/Documents/codes/bmi-cancer-code/task13/functions.R')

# install packages
source('http://bioconductor.org/biocLite.R')
install.packages("data.table")
install.packages("gplots")
install.packages("WGCNA")
install.packages("hgu133a.db")
install.packages("tidyr")
install.packages("dplyr")
install.packages("limma")
install.packages("edgeR")
install.packages("DESeq")
install.packages("affy")
install.packages("org.Hs.eg.db")
install.packages("KEGG.db")
install.packages("GO.db")
install.packages("reactome.db")

# load libraries:

library(data.table)
library(gplots)

#libraries for data wrangling
library(WGCNA)
library(hgu133a.db)
library(tidyr)
library(dplyr)

#libraries for DEG analysis
library(limma)
library(edgeR)
library(DESeq)
library(affy)

#libraries for pathway analysis
library(org.Hs.eg.db)
library(KEGG.db)
library(GO.db)
library(reactome.db)

##for venn diagram
source('http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/overLapper.R')

######################################################################

## First load all of the necessary data into R:

## Creighton et al data:





## Fuentes-Mattei et al data:





## ICGC data:





## Pathway data base:






