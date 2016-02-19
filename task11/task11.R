######################################################################

                        ##R code for task 11

######################################################################

## Do similar stuff as task 10, but with pathway enrichment analysis.

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
######################################################################

## R code from Tom:

##Need to search for the pathways using the gene symbols
## obtained from DEG (or somewhere else).

## Search/match the gene symbols with the entrez ID by pulling out the
## ID using the ****.list variable.

# Once you know the entrez ID, you can search the pathways that match
# the entrez ID, using the path* variables.

######################################################################

#Import Human Gene Symbols
SYMBOL.list<-as.list(org.Hs.egSYMBOL)

##Import KEGG pathways:
KEGG.list<-as.list(org.Hs.egPATH)
#name KEGG pathway lists with corresponding gene symbols names(KEGG.list)<-unlist(SYMBOL.list)
#mapping KEGG path IDs to human read pathway name
path<-as.list(KEGGPATHID2NAME)

##Import GO pathways:
GO.list<-as.list(org.Hs.egGO)
#... and reformat so matching KEGG.list
GO.list2<-list()
for(i in 1:length(GO.list))  GO.list2[[i]]<-names(GO.list[[i]])
rm(GO.list)
#name GO pathway lists with corresponding gene symbols
names(GO.list2)<-SYMBOL.list[names(GO.list2)]
#mapping GO IDs to human read pathway name
path2<-as.list(GOTERM)
#... and reformat so matching KEGG.list
pathx<-list()
for(i in 1:length(path2)) pathx[[i]]<-path2[[i]]@Term
names(pathx)<-names(path2)
path2<-pathx
rm(pathx)

#Import Human Reactome pathways
reactome.list<-as.list(reactomeEXTID2PATHID)
#name reactome pathway lists with corresponding gene symbols
reactome.list<-reactome.list[match(names(SYMBOL.list), names(reactome.list))]
names(reactome.list)<-unlist(SYMBOL.list)
#mapping readtome path IDs to human read pathway name
path3<-as.list(reactomePATHID2NAME)
#extract human pathways
path3<-lapply(path3[grep("Homo sapiens", path4)], function(x) strsplit(x, split=": ")[[1]][2:length(strsplit(x, split=": ")[[1]])])
path3<-lapply(path3, function(x) if(length(x)>1) paste0(x[1], ": ", x[2]) else x)













































