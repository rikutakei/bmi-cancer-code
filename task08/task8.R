######################################################################

                        ##R code for task8

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

#some genes aren't in the cancer data, so remove these genes from the obesity-related genes:
crobsgene = crobsgene[which(crobsgene %in% colnames(BLCAmat))]
crobsgene = crobsgene[order(crobsgene)] #order the obesity related genes.
crobsgene = crobsgene[-(which(crobsgene =='SLC17A6'))] #This gene produced NaN in READ cancer data during standardisation - remove it from the gene list

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
# make sure that the cancer data have genes as their columns, not rows
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

#load clinical data for Creighton et al data.
crclin = read.csv('../task05/GSE24185_clinical.csv')

#standardise Creighton et al data without logging (already taken care by RMA)
crobsmat = t(crobsmat) #make sure that the columns are the genes
crobsmat = standardise_data(crobsmat, log = F)

#function that returns the metagene from a matrix
#make sure the samples are columns
make_meta_gene = function(x) {
    s = svd(x)
    meta = rank(s$v[,1])/ncol(x)
    meta = 1-meta
    return(meta)
}

crmeta = make_meta_gene(crobsmat)
ord = order(crmeta) #order the metagene so it looks better in heatmaps.

#Make a heatmap using the data (x), metagene (meta)
make_heatmap = function(x, meta, title = '') {
    ord = order(meta)
    heatmap.2(x[,ord], trace = 'none', scale = 'none', col = 'bluered', ColSideColors = bluered(length(meta))[rank(meta)][ord], Colv=F, Rowv=T, main=title)
}

#Make box plot and dot plot with the clinical data
# if the clinical variable name is not bmi_status (e.g. Creighton has bmiStatus), change it to fit the function
#name variable is used for the title of the plots
make_plots = function(meta, clin, name = '', variable= 'bmi_status') {
    boxplot(meta~clin[[variable]], main = paste(name, 'metagene vs. BMI status'))
    plot(clin$bmi, meta, pch = 20, main = paste(name, 'metagene vs. BMI value'))
}

make_heatmap(crobsmat, crmeta)
make_plots(1-crmeta, crclin, name='Creighton', variable = 'bmiStatus')

#make transformation matrix
s = svd(crobsmat)
crtransmat = diag(1/s$d) %*% t(s$u)

#identify which samples are obese/overweight/normal weight
for (i in 1:length(bmifiles)) {
    tmp = get(bmifiles[i])
    tmp = as.data.frame(tmp)
    x = vector()
    for(j in 1:nrow(tmp)) {
        if(tmp[j,3] >= 30) {
            x = append(x,'obese')
        } else if(tmp[j,3] <= 25) {
            x = append(x,'normal')
        } else {
            x = append(x,'overweight')
        }
    }
    tmp = cbind(tmp,x)
    colnames(tmp)[4] = 'bmi_status'
    assign(bmifiles[i],tmp)
}

#pull out samples that are in both cancer and bmi data
files = gsub('cr','',files)

for (i in 1:length(files)) {
    x = rownames(get(paste(files[i],'cr',sep=''))) ##rownames of obsmat
    y = rownames(get(gsub('mat','bmi',files[i]))) ##rownames of bmi
    if(length(x) > length(y)) {
        ind = which(x %in% y)
        assign(paste(files[i],'obs',sep=''), get(paste(files[i],'cr',sep=''))[ind,])
    } else {
        ind = which(y %in% x)
        assign(gsub('mat','bmi',files[i]), get(gsub('mat','bmi',files[i]))[ind,])
    }
}

#standardise cancer data
files = paste(files,'cr',sep='')

for (i in 1:length(files)) {
    t = get(files[i])
    t = standardise_data(t)
    txt = paste(files[i],'std',sep='')
    assign(txt,t)
}

#get metagenes for all the cancer types.
files = paste(files, 'std',sep='')
metafiles = gsub('matcrstd','meta',files)

for (i in 1:length(files)) {
    txt = metafiles[i]
    t = t(crtransmat %*% get(files[i])) #transform the matrix
    t = t[,1] #get metagene
    assign(txt,t) #save it
}

#make heatmaps and plots for all the cancer types:
pdf('crcancer.pdf')
for (i in 1:length(files)) {
    dat = get(files[i])
    meta = get(metafiles[i])
    bmi = get(bmifiles[i])
    txt = gsub('meta', ' (Creighton transformed)', metafiles[i])
    make_heatmap(dat, meta, title = txt)
    make_plots(meta, bmi, name = gsub('meta','', metafiles[i]))
}
dev.off()

######################################################################
# Analysis using Fuentes-Mattei et al data
######################################################################

## Read in Fuentes-Mattei data:
fmfiles = readLines('../task05/fmraw/files.txt')
fmfiles = paste('../task05/fmraw/',fmfiles,sep='')
fmraw = ReadAffy(filenames = fmfiles)
fmraw = rma(fmraw) ## RMA normalise the data
fmraw = exprs(fmraw) ## change the format into matrix

fmobsgene = read.csv('../task05/supTab3-fixed.csv') ## list of probe IDs associated with obesity
fmprobe = fmobsgene$Probe
fmobsmat = fmraw[fmprobe,] #Pull out the obesity related gene probes

#there are duplicates, so collapse it:
collapse = collapseRows(fmobsmat, rowGroup = as.vector(fmobsgene$Gene.Symbol), rowID = rownames(fmobsmat))
fmobsmat = as.matrix(collapse$datETcollapsed) #111 genes with 278 samples
fmobsgene = rownames(fmobsmat)
fmobsgene = fmobsgene[which(fmobsgene %in% colnames(BLCAmat))] #107 genes are in both FM data and the cancer data
fmobsmat = fmobsmat[fmobsgene,] #there's a " " gene - remove it

#standardise data:
fmobsmat = t(fmobsmat) #so the genes are columns
fmobsmat = standardise_data(fmobsmat, log = F)
fmmeta = make_meta_gene(fmobsmat)
make_heatmap(fmobsmat,fmmeta,title = 'Fuentes-Mattei') #make heatmap

fmtransmat = diag(1/svd(fmobsmat)$d) %*% t(svd(fmobsmat)$u)

#pre-process cancer data
files = gsub('crstd', '', files)
for (i in 1:length(files)) {
    txt = paste(files[i], 'fm', sep='')
    t = get(files[i])
    t = t[,fmobsgene]
    assign(txt, t)
}

#standardise and make metagene from the cancer data
files = paste(files, 'fm',sep = '')
for (i in 1:length(files)) {
    #standardise
    txt = paste(files[i],'std',sep='')
    t = get(files[i])
    t = standardise_data(t)
    assign(txt, t)

    #get metagene
    txt = gsub('matfmstd','fmmeta',txt)
    t = t(fmtransmat %*% t)
    t = t[,1]
    assign(txt,t)
}

#make heatmaps and plots
files = paste(files,'std',sep='')
metafiles = gsub('meta','fmmeta',metafiles)
pdf('fmcancer.pdf')
for (i in 1:length(files)) {
    dat = get(files[i])
    bmi = get(bmifiles[i])
    meta = get(metafiles[i])
    txt = gsub('fmmeta', ' (Fuentes-Mattei transformed)', metafiles[i])
    make_heatmap(dat, meta, title= txt)
    make_plots(meta,bmi,name = gsub('fmmeta','',metafiles[i]))
}
dev.off()






































