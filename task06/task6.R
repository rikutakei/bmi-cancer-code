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
data = c('BLCAmat','CESCmat','COADmat','KIRPmat','LIHCmat','READmat','SKCMmat','UCECmat')
for(i in 1:length(data)) {
    assign(data[i],icgc_to_tcga(get(data[i])))
}

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

##Need to remove duplicate samples in the data: (use WGCNA)
tmp = BLCAmat
tmp = collapseRows(tmp, unique(rownames(tmp)), unique(rownames(tmp)))
BLCAmat = tmp$datETcollapsed

##load Creighton data:
files = readLines('../task05/creightonraw/files.txt')
files = paste('../task05/creightonraw/', files, sep='')
raw = ReadAffy(filenames = files)
raw = rma(raw)
raw = exprs(raw)
colnames(raw) = gsub('.CEL','', colnames(raw))

obsgene = readLines('../task05/obsGenes.txt')
obsgenemat = raw[obsgene,]

genesymbol = mapIds(hgu133a.db, keys = obsgene, column = 'SYMBOL', keytype = "PROBEID", multiVals = 'first')
rownames(obsgenemat) = as.vector(genesymbol)

##remove NAs and duplicate genes
obsgenemat = obsgenemat[-(which(is.na(rownames(obsgenemat)))),]
obsgenemat = collapseRows(obsgenemat, rowGroup = unique(rownames(obsgenemat)), rowID = unique(rownames(obsgenemat)))
obsgenemat = obsgenemat$datETcollapsed

##I need to find which of the 695 gene symbols are present in each of the cancer type.
##the number of genes and the symbols seems to match between different cancer types:
dim(BLCAmat) ## 20501 genes in both BLCA and CESC
dim(CESCmat)
which((colnames(BLCAmat) %in% colnames(CESCmat)) == F) ##integer(0), so all 20501 genes match

##find which genes are in both data sets
tmp = as.vector(colnames(BLCAmat))
genelist = which(unique(genesymbol) %in% tmp)
genelist = unique(genesymbol)[genelist]

length(genelist) ##647 unique genes

genelist = genelist[-which(genelist == 'SLC17A6')]

obsgenemat = obsgenemat[genelist,]
crobsraw = obsgenemat

##Don't have to log the data as it has been rma normalised

##standardise the count data:
obsgenemat = t(apply(obsgenemat, 1, function(x) (x-mean(x))/sd(x)))
obsgenemat[obsgenemat > 3] = 3
obsgenemat[obsgenemat < -3] = -3

heatmap.2(obsgenemat, trace='none',scale='none', col='bluered')

##make metagene and check that Creighton data correlates with BMI
s = svd(obsgenemat)
crmeta = rank(s$v[,1])/ncol(obsgenemat)
ord = order(crmeta)

sraw = svd(crobsraw)
crmeta2 = -sraw$v[,1]

heatmap.2(obsgenemat, trace='none',scale='none', col='bluered', ColSideColors = bluered(length(crmeta))[rank(crmeta)])
heatmap.2(obsgenemat[,ord], trace='none',scale='none', col='bluered', ColSideColors = bluered(length(crmeta))[rank(crmeta)][ord], Colv=F, Rowv=T)

crclin = read.csv('GSE24185_clinical.csv')
par(mfrow = c(2,2))

boxplot(crmeta~crclin$bmiStatus)
plot(crclin$bmi, crmeta, pch = 20)

boxplot(crmeta2~crclin$bmiStatus)
plot(crclin$bmi, crmeta2, pch = 20)

par(mfrow = c(1,1))

##pull out the 647 genes from each cancer type
data = c('BLCAmat','CESCmat','COADmat','KIRPmat','LIHCmat','READmat','SKCMmat','UCECmat')

for (i in 1:length(data)) {
    x = get(data[i])
    x = x[,genelist]
    txt = paste('obs', data[i], sep='')
    assign(txt, x)
}

##need to apply transformation matrix to other cancer types
## make transformation matrix:
crtransmat = diag(1/s$d) %*% t(s$u) ## svd from normalised data
crrawtransmat = diag(1/sraw$d) %*% t(sraw$u)

##identify who's obese/overweight/normal weight in each cancer type (and put in data frame)
data = c('BLCAmat','CESCmat','COADmat','KIRPmat','LIHCmat','READmat','SKCMmat','UCECmat')
data = gsub('mat','bmi',data)

for (i in 1:length(data)) {
    tmp = get(data[i])
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
    assign(data[i],tmp)
}

##need to remove the samples that are not in the matrix
data = c('BLCAmat','CESCmat','COADmat','KIRPmat','LIHCmat','READmat','SKCMmat','UCECmat')

for (i in 1:length(data)) {
    x = rownames(get(paste('obs',data[i],sep=''))) ##rownames of obsmat
    y = rownames(get(gsub('mat','bmi',data[i]))) ##rownames of bmi
    if(length(x) > length(y)) {
        ind = which(x %in% y)
        assign(paste('obs',data[i],sep=''), get(paste('obs',data[i],sep=''))[ind,])
    } else {
        ind = which(y %in% x)
        assign(gsub('mat','bmi',data[i]), get(gsub('mat','bmi',data[i]))[ind,])
    }
}

##apply transformation matrix and get metagene for each cancer type

##standardise cancer data
data = paste('obs',data,sep='')
for(i in 1:length(data)) {
    x = get(data[i])
    x = log10(x + 1)
    x = t(apply(x, 2, function(x) (x-mean(x))/sd(x)))
    x[x < -3] = -3
    x[x > 3] = 3
    heatmap.2(x, trace='non',scale='none',col='bluered',main=data[i])
    txt = gsub('obs','std',data[i])
    assign(txt,x)
}

##obtain metagene:
data = c('BLCAmat','CESCmat','COADmat','KIRPmat','LIHCmat','READmat','SKCMmat','UCECmat')
for(i in 1:length(data)) {
    txt = gsub('mat','meta',data[i])
    x = get(paste('std',data[i],sep = ''))
    tmpmeta = t(crtransmat %*% x)
    tmpmeta = tmpmeta[,1]
    assign(txt,tmpmeta)
}

##make heatmaps with the metagene:
for(i in 1:length(data)) {
    tmp = get(paste('std',data[i],sep=''))
    m = get(gsub('mat','meta',data[i]))/ncol(tmp)
    r = rank(m)
    o = order(r)

    heatmap.2(tmp[,o], trace='none',scale='none', col='bluered', ColSideColors = bluered(length(r))[rank(r)][o], Colv=F, Rowv=T,main = data[i])
}

##make plots
for(i in 1:length(data)) {
    tmp = get(paste('std',data[i],sep=''))
    m = get(gsub('mat','meta',data[i]))/ncol(tmp)
    r = rank(m)
    clin = get(gsub('mat','bmi',data[i]))

    boxplot(r~clin$bmi_status,main = data[i])
    plot(clin$bmi, r, pch = 20, main = data[i])
}

######################################################################

Fuentes-Mattei stuff

######################################################################

##look at task5.R and use the first half of the code
##make sure you change some variable names to prevent overwriting.

##post-Fuentes-Mattei processing:
## make new matrices using the genes used for Fuentes-Mattei:
fmgenes = as.vector(tmp$Gene.Symbol)
fmgenes = fmgenes[-(which(fmgenes == " "))] ##remove bad ones
data = c('BLCAmat','CESCmat','COADmat','KIRPmat','LIHCmat','READmat','SKCMmat','UCECmat')
for(i in 1:length(data)) {
    a = get(data[i])
    b = t(a[,which(colnames(a) %in% fmgenes)])
    assign(paste('fmobs',data[i], sep=''), b)
}

##standardise data:
data = paste('fmobs',data,sep='')
for(i in 1:length(data)) {
    x = get(data[i])
    x = log10(x + 1)
    x = t(apply(x, 1, function(x) (x-mean(x))/sd(x)))
    x[x < -3] = -3
    x[x > 3] = 3
    heatmap.2(x, trace='non',scale='none',col='bluered',main=data[i])
    txt = gsub('obs','std',data[i])
    assign(txt,x)
}

## apply transformation matrix to the cancer data:

##make sure that the transformation matrix has the same genes as other cancer types:
fmobegenemat = fmobegenemat[rownames(fmobsBLCAmat),]
fmsvd = svd(fmobegenemat)
fmtransmat = diag(1/fmsvd$d) %*% t(fmsvd$u)
for(i in 1:length(data)) {
    txt = gsub('mat','meta',data[i])
    x = get(gsub('obs','std', data[i]))
    tmpmeta = t(fmtransmat %*% x)
    tmpmeta = tmpmeta[,1]
    tmpmeta = 1-tmpmeta
    assign(txt,tmpmeta)
}

##make heatmaps with the metagene:
for(i in 1:length(data)) {
    tmp = get(gsub('obs','std', data[i]))
    m = get(gsub('mat','meta',data[i]))/ncol(tmp)
    r = rank(m)
    o = order(r)

    heatmap.2(tmp[,o], trace='none',scale='none', col='bluered', ColSideColors = bluered(length(r))[rank(r)][o], Colv=F, Rowv=T,main = data[i])
}

##make plots
for(i in 1:length(data)) {
    tmp = get(gsub('obs','std', data[i]))
    m = get(gsub('mat','meta',data[i]))/ncol(tmp)
    r = rank(m)
    clin = gsub('mat','bmi',data[i])
    clin = get(gsub('fmobs','',clin))

    boxplot(r~clin$bmi_status,main = data[i])
    plot(clin$bmi, r, pch = 20, main = data[i])
}


