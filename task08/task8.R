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








