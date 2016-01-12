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

for (i in 1:length(bmifiles)) {
    txt = gsub('.txt','',bmifiles[i])
    txt = gsub('../task06/','',txt)
    assign(txt, dget(bmifiles[i]))
}




