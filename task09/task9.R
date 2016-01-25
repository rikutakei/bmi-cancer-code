######################################################################

                        ##R code for task9

######################################################################

## Start doing differential gene expression analysis with all the
## cancer types and find commonly differentially expressed genes.

######################################################################
## pre-analysis stuff
setwd('~/Documents/masters/data/task09/')

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

#read in Creighton et al data:
crfiles = readLines('../task05/creightonraw/files.txt')
crfiles = paste('../task05/creightonraw/', crfiles, sep='')
raw = ReadAffy(filenames = crfiles)
raw = rma(raw)
raw = exprs(raw)
colnames(raw) = gsub('.CEL','', colnames(raw))
#crobsmat = raw

#load clinical data for Creighton et al data.
crclin = read.csv('../task05/GSE24185_clinical.csv')

#start DEG analysis with Creighton et al data
#Need to make a model matrix first:
group = crclin$bmiStatus

#make >overweight group and >obese group
groupcrov = ifelse(group == 'normal', 'normal', 'overweight/obese')
groupcrob = ifelse(group == 'obese',  'obese','normal/overweight')

#make model matrix
modelcrov = model.matrix(~groupcrov)
colnames(modelcrov) = c('normal','overweight/obese')
modelcrob = model.matrix(~groupcrob)
colnames(modelcrob) = c('normal/overweight','obese')

#fit linear model
crovfit = lmFit(raw,modelcrov)
crovfit = eBayes(crovfit)
crovtt = topTable(crovfit, coef = 'overweight/obese', adjust = 'BH', n = nrow(raw))

crobfit = lmFit(raw,modelcrob)
crobfit = eBayes(crobfit)
crobtt = topTable(crobfit, coef = 'obese', adjust = 'BH', n = nrow(raw))

#see how many genes are significant:
sum(crovtt$P.Value < 0.05) #759 genes
sum(crovtt$P.Value < 0.01) #104 genes
sum(crobtt$P.Value < 0.05) #5278 genes
sum(crobtt$P.Value < 0.01) #1781 genes

sum(crovtt$adj.P.Val < 0.05) #0 genes
sum(crovtt$adj.P.Val < 0.01) #0 genes
sum(crobtt$adj.P.Val < 0.05) #9 genes
sum(crobtt$adj.P.Val < 0.01) #0 genes

#use this to convert gene probes into gene symbols
#crgenes = mapIds(hgu133a.db, keys = crgenes, column = 'SYMBOL', keytype = 'PROBEID', multiVals = 'first')





#read in cancer files
files = readLines('../task08/files.txt')

#create variables for cancer data matrix
for (i in 1:length(files)) {
    txt = gsub('.txt','',files[i])
    files[i] = txt #rename files to its variable names
    assign(txt, dget(files[i]))
}

icgc_to_tcga = function(x) {
    tcgaid = gsub(pattern = '-[[:alnum:]]{3}-[[:alnum:]]{3}-[[:alnum:]]{4}-[[:alnum:]]{2}',
    replacement = '', rownames(x))
    rownames(x) = tcgaid
    return(x)
}

















