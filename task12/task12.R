######################################################################

                        ##R code for task 12

######################################################################

## Go back to Creighton et al data and double check if any other
## clinical variables are affecting the BMI marker.

######################################################################
## pre-analysis stuff
setwd('~/Documents/masters/data/task11/')

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

#pathway analysis packages from bioconductor
source('http://bioconductor.org/biocLite.R')
biocLite(c('topGO','GOstats','ReactomePA','KEGGprofile','PathNet','safe'))
library(topGO)
library(GOstats)
library(ReactomePA)
library(KEGGprofile)
library(PathNet)
library(safe)

