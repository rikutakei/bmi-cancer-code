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
files = paste('../task08/',files,sep='')

#create variables for cancer data matrix
for (i in 1:length(files)) {
    txt = gsub('.txt','',files[i])
    txt = gsub('../task08/','',txt)
    assign(txt, dget(files[i]))
    files[i] = txt #rename files to its variable names
}

icgc_to_tcga = function(x) {
    tcgaid = gsub(pattern = '-[[:alnum:]]{3}-[[:alnum:]]{3}-[[:alnum:]]{4}-[[:alnum:]]{2}',
    replacement = '', rownames(x))
    rownames(x) = tcgaid
    return(x)
}

#get BMI data
bmifiles = gsub('mat','bmi',files) #change mat to bmi
bmifiles = paste('../task06/',bmifiles,sep='')
bmifiles = paste(bmifiles,'.txt',sep='')

#read in all the bmi data and give it a variable name
for (i in 1:length(bmifiles)) {
    txt = gsub('.txt','',bmifiles[i])
    txt = gsub('../task06/','',txt)
    assign(txt, dget(bmifiles[i]))
    bmifiles[i] = txt
}

#rename the sample names
for (i in 1:length(files )) {
    t = get(files [i])
    t = t(t) ##transpose the matrix so samples ids are the rows, not columns
    t = icgc_to_tcga(t)
    assign(files[i],t)
}

##identify which samples have bmi data for every cancer type
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

#need to pull out only the samples present in the cancer data to make a design matrix
for (i in 1:length(bmifiles)) {
    t = rownames(get(files[i]))
    bmi = get(bmifiles[i])
    assign(bmifiles[i], bmi[t,])
}

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

#Need to make a model matrix for all the cancer types:

for (i in 1:length(bmifiles)) {
    group = get(bmifiles[i])
    group = group[,4]
    txt = gsub('bmi','',bmifiles[i])
    txt = paste('group',txt,sep='')
    txt = paste(txt,'ov',sep='')
    assign(txt, ifelse(group == 'normal', 'normal', 'overweight/obese'))

    model = get(txt)
    model = model.matrix(~model)
    colnames(model) = c('normal','overweight/obese')
    assign(gsub('group','model',txt), model)

    txt = gsub('ov','ob',txt)
    assign(txt, ifelse(group == 'obese', 'obese', 'normal/overweight'))
    model = get(txt)
    model = model.matrix(~model)
    colnames(model) = c('normal/overweight','obese')
    assign(gsub('group','model',txt), model)
}

#Need to normalise the data using voom:
#The RNA-seq data has to be genes by samples
normVoom = function(x, design) {
    a = x[rowSums(cpm(x)) > 9,]
    b = calcNormFactors(a)
    a = voom(a, design, lib.size = colSums(a)*b)
    return(a)
}

#function to make toptable, given an expression  matrix (gene by samples) and a corresponding model matrix.
make_tt = function(x, model) {
    fit = lmFit(x, model)
    fit = eBayes(fit)
    tt = topTable(fit, coef = 2, adjust = 'BH', n = nrow(x))
    return(tt)
}

#voom normalise the data and make top table
modelfiles = paste('model',files, sep='')
modelfiles = gsub('mat','',modelfiles)

for (i in 1:length(files)) {
    txt = gsub('mat','ovtt',files[i])
    mat = t(get(files[i]))
    model = get(paste(modelfiles[i],'ov',sep=''))
    mat = normVoom(mat,model)
    t = make_tt(mat, model)
    assign(txt, t)

    txt = gsub('mat','obtt',files[i])
    mat = t(get(files[i]))
    model = get(paste(modelfiles[i],'ob',sep=''))
    mat = normVoom(mat,model)
    t = make_tt(mat, model)
    assign(txt, t)
}

#make a function to pull out all the significant DEGs from the top table
# x is the toptable, adj is whether to use adjusted or unadjusted p value, and y is the p value to be used.
pull_deg = function(x, adj = F, y = 0.05) {
    if(adj) {
        x = x[x$adj.P.Val <= y,]
    } else {
        x = x[x$P.Value <= y,]
    }
    return(x)
}

#Put all the DEGs in a single list
oblist = as.list(gsub('mat','',files))
names(oblist) = gsub('mat','',files)

ovlist = as.list(gsub('mat','',files))
names(ovlist) = gsub('mat','',files)

#for p < 0.05
for (i in 1:length(files)) {
    t = get(gsub('mat','obtt',files[i]))
    oblist[[i]] = rownames(pull_deg(t))

    t = get(gsub('mat','ovtt',files[i]))
    ovlist[[i]] = rownames(pull_deg(t))
}

#for p < 0.01
for (i in 1:length(files)) {
    t = get(gsub('mat','obtt',files[i]))
    oblist[[i]] = rownames(pull_deg(t,y = 0.01))

    t = get(gsub('mat','ovtt',files[i]))
    ovlist[[i]] = rownames(pull_deg(t,y = 0.01))
}

#Make venn diagram
#source('http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/overLapper.R')
#ol = overLapper(setlist = oblist,sep='_',type='vennsets')
#counts = sapply(ol$Venn_List, length)
#vennPlot(counts = counts)

#Find all the genes that are commonly differentially expressed in all cancer types
obgenes = unlist(oblist)
names(t) = NULL
obgenes = names(which(table(obgenes) >= 3))

ovgenes = unlist(ovlist)
names(t) = NULL
ovgenes = names(which(table(ovgenes) >= 3))























