###############################################################################

                        ##R code for task 13

###############################################################################

## Combine tasks 1 to 12 into this document

###############################################################################
## pre-analysis stuff
setwd('~/Documents/masters/data/task13/')

#source the file to get all the functions from task10.R
source('~/Documents/codes/bmi-cancer-code/task13/functions.R')

# install packages
source('http://bioconductor.org/biocLite.R')
install.packages("data.table")
install.packages("gplots")
biocLite()
biocLite("WGCNA")
biocLite("hgu133a.db")
biocLite("tidyr")
biocLite("dplyr")
biocLite("limma")
biocLite("edgeR")
biocLite("DESeq")
biocLite("affy")
biocLite("org.Hs.eg.db")
biocLite("KEGG.db")
biocLite("GO.db")
biocLite("reactome.db")

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

###############################################################################

## First load all of the necessary data into R:

###############################################################################
## Creighton et al data:

files = readLines('*****/files.txt')
files = paste('****/', files, sep='')
cr_raw = ReadAffy(filenames = files)
cr_raw = rma(cr_raw) ## RMA normalise the data
cr_raw = exprs(cr_raw) ## change the format into matrix

cr_obsgene = read.csv('*****')

###############################################################################
## Fuentes-Mattei et al data:

files = readLines('*****/files.txt')
files = paste('****/', files, sep='')
fm_raw = ReadAffy(filenames = files)
fm_raw = rma(fm_raw) ## RMA normalise the data
fm_raw = exprs(fm_raw) ## change the format into matrix

fm_obsgene = read.csv('*****')

###############################################################################
## ICGC data:

## TODO: make sure you do this for all cancer types
## TODO: process the data using whatever method I have used
files = readLines('*****/files.txt')
for (i in 1:length(files)) {
	seq = read.table(files[i], sep = '\t', header = T)
}

###############################################################################
## Pathway data base:

#Import Human Gene Symbols
SYMBOL.list<-as.list(org.Hs.egSYMBOL)

KEGG.list<-as.list(org.Hs.egPATH) ##Import KEGG pathways:
names(KEGG.list)<-unlist(SYMBOL.list) #name KEGG pathway lists with corresponding gene symbols
keggpath<-as.list(KEGGPATHID2NAME) #mapping KEGG path IDs to human read pathway name

GO.list<-as.list(org.Hs.egGO) ##Import GO pathways:
tmp = GO.list #... and reformat so matching KEGG.list

# pull out the GOID
for (i in 1:length(tmp)) {
    ind = names(tmp)[i]
    if(class(tmp[[ind]]) == 'list') {
        tmp[[ind]] = names(tmp[[ind]])
    }
}
tmp = tmp[which(!is.na(tmp))]##filter out the NA values
GO.list = tmp
names(GO.list)<-SYMBOL.list[names(GO.list)] #name GO pathway lists with corresponding gene symbols
gopath<-as.list(GOTERM) #mapping GO IDs to human read pathway name
#... and reformat so matching KEGG.list
tmp = list()
for(i in 1:length(gopath)) tmp[[i]]<-gopath[[i]]@Term
names(tmp) = names(gopath)
gopath = tmp

#Import Human Reactome pathways
reactome.list <-as.list(reactomeEXTID2PATHID)
reactome.list = reactome.list[order(as.numeric(names(reactome.list)))] #Sort the names in the list:
reactome.list = reactome.list[which(names(reactome.list) %in% names(SYMBOL.list))] #get paths that have gene symbols:
names(reactome.list) = unlist(SYMBOL.list[names(reactome.list)]) #rename the entrez gene ID into gene symbol
reactomepath <- as.list(reactomePATHID2NAME) #mapping readtome path IDs to human read pathway name
reactomepath = reactomepath[grep('Homo sapiens',reactomepath)] #pull out all human-related pathways
reactomepath = lapply(reactomepath, function(x) gsub('Homo sapiens: ', '', x)) #cut out the 'Homo sapiens: ' bit so it's only the pathway names.

###############################################################################
## Creighton metagene stuff

## make matrix with just obesity-associated genes:
cr_obsmat = cr_raw[cr_obsgene,]

cr_genesymbol = mapIds(hgu133a.db, keys = cr_obsgene, column = 'SYMBOL', keytype = "PROBEID", multiVals = 'first')
rownames(cr_obsgene) = as.vector(cr_genesymbol)

##remove NAs and duplicate genes
cr_obsmat = cr_obsmat[-(which(is.na(rownames(cr_obsmat)))),]
cr_obsmat = collapseRows(cr_obsmat, rowGroup = unique(rownames(cr_obsmat)), rowID = unique(rownames(cr_obsmat)))
cr_obsmat = cr_obsmat$datETcollapsed

## TODO: Do the heatmap stuff and add p-values to the plots
cr_rawmeta = svd(cr_obsmat)
cr_rawmeta = rank(cr_rawmeta$v[,1])/ncol(cr_obsmat)
cr_rawmeta = 1-cr_rawmeta

cr_obsmat_adj = t(apply(cr_obsmat, 1, function(x) (x-mean(x))/sd(x)))
cr_obsmat_adj[cr_obsmat_adj > 3] = 3
cr_obsmat_adj[cr_obsmat_adj < -3] = -3

cr_adjmeta = svd(cr_obsmat)
cr_adjmeta = rank(cr_adjmeta$v[,1])/ncol(cr_obsmat)
cr_adjmeta = 1-cr_adjmeta

## get genes that are common from the creighton gene sets and the ICGC data set









###############################################################################
## Fuentes-Mattei metagene stuff










###############################################################################
## Creighton gene expression analysis stuff

creighton_raw



###############################################################################
## Common genes across ICGC data stuff

## use codes from task11?



###############################################################################
## Pathway enrichment stuff















