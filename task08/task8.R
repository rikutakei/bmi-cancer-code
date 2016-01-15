######################################################################

                        ##R code for task6

######################################################################

## Start analysis of other cancer data produced from tidyr

######################################################################
## pre-analysis stuff
setwd('~/Documents/masters/data/task08/')

library(data.table)
library(limma)
library(edgeR)
library(DESeq)
library(dplyr)
library(gplots)
library(affy)
library(WGCNA)
library(hgu133a.db)
library(tidyr)
######################################################################

#read in the cancer data created from tidyr
files = readLines('files.txt')

#create variables for cancer data matrix
for (i in 1:length(files)) {
    txt = gsub('.txt','',files[i])
    assign(txt, dget(files[i]))
}

#need bmi data for all the cancer data (use data from task6)
bmifiles = gsub('mat','bmi',files)
bmifiles = paste('../task06/',bmifiles,sep='')

#read in all the bmi data and give it a variable name
for (i in 1:length(bmifiles)) {
    txt = gsub('.txt','',bmifiles[i])
    txt = gsub('../task06/','',txt)
    assign(txt, dget(bmifiles[i]))
}

## function to strip ICGC IDs into TCGA ID:
## (function copied from task4)
icgc_to_tcga = function(x) {
    tcgaid = gsub(pattern = '-[[:alnum:]]{3}-[[:alnum:]]{3}-[[:alnum:]]{4}-[[:alnum:]]{2}',
                                   replacement = '', rownames(x))
    rownames(x) = tcgaid
    return(x)
}

#use icgc_to_tcga function on all the cancer data
files = gsub('.txt','',files )

for (i in 1:length(files )) {
    t = get(files [i])
    t = t(t) ##transpose the matrix so samples ids are the rows, not columns
    t = icgc_to_tcga(t)
    assign(files[i],t)
}

#identify which samples have bmi data for every cancer type
bmifiles = gsub('mat','bmi',files) #change mat to bmi

for (i in 1:length(files)) {
    bmi = get(bmifiles[i])
    cancer = get(files[i])
    v = which(rownames(cancer) %in% rownames(bmi)) #find which sample ids overlap in both cancer and bmi data
    cancer = cancer[v,]
    assign(files[i], cancer)
}

#find/merge duplicate samples:
for (i in 1:length(files)) {
    t = get(files[i])

    if (length(rownames(t)) > length(unique(rownames(t)))) {
        t = collapseRows(t, unique(rownames(t)), unique(rownames(t)))
        t = t$datETcollapsed
        print(paste(files[i],'had duplicates'))
    }

    assign(files[i],t)
}

#need to know the obesity genes to make metagene (import Creighton et al obesity gene)
crobsgene = readLines('../task05/obsGenes.txt') #import Creighton obesity genes

#convert it into gene symbols:
crobsgene = mapIds(hgu133a.db, keys = crobsgene, column = 'SYMBOL', keytype = "PROBEID", multiVals = 'first')

crobsgene = crobsgene[-(which(is.na(crobsgene)))] #remove any NAs
crobsgene = unique(crobsgene) #include only unique genes.

#some genes aren't int he cancer data, so remove these genes from the obesity-related genes:
crobsgene = crobsgene[which(crobsgene %in% colnames(BLCAmat))]
crobsgene = crobsgene[order(crobsgene)] #order the obesity related genes.

#pull out the obesity related genes from each cancer type
for (i in 1:length(files)) {
    t = get(files[i])
    t = t[,crobsgene]
    txt = paste(files[i],'cr',sep='')
    assign(txt,t)
}

#start SVD analysis
files = paste(files,'cr',sep='')

#Make a function to log and standardise cancer data
# make sure that the cancer data have samples as their columns, not rows
# if the matrix to be standardised doesn't have to be logged, then set the log variable to F
standardise_data = function(x, log = T) {
    if (log) {
        x = log10(x + 1)
    }
    x = t(apply(x, 2, function(x) (x-mean(x))/sd(x)))
    x[x < -3] = -3
    x[x > 3] = 3
    return(x)
}

#need to get metagene from Creighton et al data first:
crfiles = readLines('../task05/creightonraw/files.txt')
crfiles = paste('../task05/creightonraw/', crfiles, sep='')
raw = ReadAffy(filenames = crfiles)
raw = rma(raw)
raw = exprs(raw)
colnames(raw) = gsub('.CEL','', colnames(raw))
crobsmat = raw

#need to pull out the obesity related genes in Creighton data
crgenes = rownames(crobsmat)
crgenes = mapIds(hgu133a.db, keys = crgenes, column = 'SYMBOL', keytype = 'PROBEID', multiVals = 'first')
rownames(crobsmat) = as.vector(crgenes)
crobsmat = crobsmat[crobsgene,] #pull out the onesity related genes.











