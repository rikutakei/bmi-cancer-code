######################################################################

                        ##R code for task5

######################################################################

## New dataset (Fuentes-Mattei et al, 2014) to see whether metagene correlates with BMI

######################################################################
## pre-analysis stuff
setwd('~/Documents/masters/data/task05/')

library(limma)
library(edgeR)
library(DESeq)
library(dplyr)
library(gplots)
library(affy)
library(WGCNA)
library(hgu133a.db)
######################################################################

## Read in Fuentes-Mattei data:
files = readLines('fmraw/files.txt')
files = paste('fmraw/',files,sep='')
raw = ReadAffy(filenames = files)
raw = rma(raw) ## RMA normalise the data
raw = exprs(raw) ## change the format into matrix

obsgene = read.csv('supTab3-fixed.csv') ## list of probe IDs associated with obesity

probeid = obsgene$Probe ## get Probe IDs relevant to obesity
obsgenemat = raw[probeid,] ## get the obesity related probes

## There are 130 gene probes, but only 111 unique genes
unique(obsgene$Gene.Symbol) %>% length ## 111 genes
length(probeid) ## 130 probes

## Collapse the repeated genes using collapseRows from WGCNA library
collapse = collapseRows(obsgenemat, rowGroup=as.vector(obsgene$Gene.Symbol),rowID=rownames(obsgenemat))
obsgenemat = as.matrix(collapse$datETcollapsed)

obsgenemat = t(apply(obsgenemat, 1, function(x) (x-mean(x))/sd(x)))
obsgenemat[obsgenemat > 3] = 3
obsgenemat[obsgenemat < -3] = -3

heatmap.2(obsgenemat,trace='none',scale='none',col='bluered')

## Apply SVD to obsgenemat:
svd = svd(obsgenemat)
meta = rank(svd$v[,1])/ncol(obsgenemat)
ord = order(meta)

heatmap.2(obsgenemat[,ord],trace='none',scale='none',col='bluered', ColSideColors = bluered(length(meta))[rank(meta)][ord],Colv=F,Rowv=T)

## Make transformation matrix:
transmat = diag(1/svd$d) %*% t(svd$u)

######################################################################

## Transform Creighton data with Fuentes-Mattei metagene

######################################################################

## Read in Cregihton data:
files = readLines('creightonraw/files.txt')
files = paste('creightonraw/',files,sep='')
crraw = ReadAffy(filenames = files)
crraw = rma(crraw) ## RMA normalise the data
crraw = exprs(crraw) ## change the format into matrix

## make matrix using the same obesity probe IDs as in Fuentes-Mattei data
crmat = crraw[probeid,]

## Collapse the repeated genes using collapseRows from WGCNA library
collapse = collapseRows(crmat, rowGroup=as.vector(obsgene$Gene.Symbol),rowID=rownames(crmat))
crmat = as.matrix(collapse$datETcollapsed)

crmat = t(apply(crmat, 1, function(x) (x-mean(x))/sd(x)))
crmat[crmat > 3] = 3
crmat[crmat < -3] = -3

heatmap.2(crmat,trace='none',scale='none',col='bluered')

crmeta = t(transmat %*% crmat) ## apply transformation matrix
crmeta = crmeta[,1]
crmeta = rank(crmeta)/ncol(crmat)
ord = order(crmeta)

heatmap.2(crmat,trace='none',scale='none',col='bluered',ColSideColors=bluered(length(crmeta))[rank(crmeta)])

heatmap.2(crmat[,ord],trace='none',scale='none',col='bluered',ColSideColors=bluered(length(crmeta))[rank(crmeta)][ord],Colv=F,Rowv=T)

## load clinical data for Creighton data
crclin = read.csv('GSE24185_clinical.csv')







