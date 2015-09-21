######################################################################

                        ##R code for task6

######################################################################

## Now looking at other cancer data.

######################################################################
## pre-analysis stuff
setwd('~/Documents/masters/data/task06/')

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

##load RNA seq data:
file = readLines('raw/names.txt')
file = paste('raw/',file,sep='')

for (i in 1:length(file)) {
    txt = gsub('raw/exp_seq.','',file[i])
    txt = gsub('-US.tsv','',txt)
    txt = paste(txt, 'raw',sep='')

    assign(txt, read.table(file[i], sep = '\t', header = T))
}

seq = tbl_df(seq) ## change data frame into tbl class

## extract all the raw count data and put it in a matrix:
icgcToMatrix = function(icgc_seq) {
    seq = tbl_df(icgc_seq)
    seq = select(seq, submitted_sample_id, gene_id, raw_read_count)
    initsamplenames = as.vector(unique(seq$submitted_sample_id))
    initgenenames = as.vector(unique(seq$gene_id))
    if ('?' %in% initgenenames) {
        initgenenames = initgenenames[-which(initgenenames == '?')]
        seq = filter(seq, gene_id %in% initgenenames)
    }
    raw = matrix(nrow = length(initsamplenames), ncol = length(initgenenames))
    rownames(raw) = initsamplenames
    colnames(raw) = initgenenames

    for (i in 1:length(initsamplenames)) {
        sample = initsamplenames[i]
        tmp = filter(seq, submitted_sample_id == sample)[,2:3]

        for (j in 1:length(initgenenames)) {
            gene = initgenenames[j]
            tmp2 = filter(tmp, gene_id == gene)[,2]
            #print(tmp2)
            #print(j)
            if (dim(tmp2)[1] > 1) {
                raw[sample, gene] = as.numeric(max(tmp2))
            } else {
                raw[sample, gene] = as.numeric(tmp2)
            }
            #raw[sample, gene] = as.numeric(tmp2)
        }
        print(paste('Finished processing ', sample))
    }

    return(raw)
}

## Need to do icgcToMatrix to all the datasets:
data = c('BLCAraw','CESCraw','COADraw','KIRPraw','LIHCraw','READraw','SKCMraw')

for (i in 1:length(data)) {
    txt = gsub('raw', 'mat', data[i])
    tmp = get(data[i])
    tmp = icgcToMatrix(tmp)
    assign(txt, tmp)
    ## save the matrix:
    dput(get(txt), file = paste(data[i],'readcount.txt',sep =''))
}


## save the matrix:
## dput(raw, file = 'raw_read_count.txt')

##run icgc_to_tcga function from task 4

## load BMI_data.txt from task1
all_bmi = dget(file = 'BMI_data.txt')

names(all_bmi) ## all the names of the cancer in this list

## get BMI info from this list
tmpbmi = all_bmi$blca ##substitute the cancer type for a different one
BLCAbmi = matrix(nrow = nrow(tmpbmi), ncol = 3)
rownames(BLCAbmi) = tmpbmi$bcr_patient_barcode
BLCAbmi[,1] = as.numeric(as.character(tmpbmi$weight_kg_at_diagnosis))
BLCAbmi[,2] = as.numeric(as.character(tmpbmi$height_cm_at_diagnosis))/100
BLCAbmi[,3] = BLCAbmi[,1]/(BLCAbmi[,2]^2)
colnames(BLCAbmi) = c('weight', 'height', 'bmi')
dput(BLCAbmi, file = 'BLCAbmi.txt')

##identify the samples that have BMI data and pull it out
##(substitute the cancer code for a different cancer type)
v = which(rownames(BLCAmat) %in% rownames(BLCAbmi))
BLCAmat = BLCAmat[v,]








