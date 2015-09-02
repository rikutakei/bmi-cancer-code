######################################################################

                        ##R code for task4

######################################################################

## Re-do the analysis with logged data

######################################################################
## pre-analysis stuff
setwd('~/Documents/masters/data/task04/')

library(limma)
library(edgeR)
library(DESeq)
library(dplyr)
library(gplots)
######################################################################

raw = dget('raw_read_count.txt')
ucecbmi = dget('hwdata.txt')

## function to strip ICGC IDs into TCGA ID:
icgc_to_tcga = function(x) {
    tcgaid = gsub(pattern = '-[[:alnum:]]{3}-[[:alnum:]]{3}-[[:alnum:]]{4}-[[:alnum:]]{2}',
                                   replacement = '', rownames(x))
    rownames(x) = tcgaid
    return(x)
}

raw = icgc_to_tcga(raw)

## get TCGA IDs that have height/weight data.
x = which(rownames(ucecbmi) %in% rownames(raw))
ucecbmi = ucecbmi[x,]
raw = raw[as.vector(rownames(ucecbmi)),]

## any duplicates?
dim(ucecbmi) ## 482 samples
dim(raw) ## 482 samples, so no duplicates.
unique(rownames(raw)) %>% length ## double checking

## Quality control - remove 'bad' smaples:
finalsample = t(raw)
finalsample = finalsample + 1 ## to account for the 0 values
finalsample = log10(finalsample) ## LOG THE DATA!!!
## repeat this a few times, if low correlation (obviously with different values in cutree, etc.)
c = cor(finalsample, use = 'all.obs', method = 'pearson')
heatmap.2(c, trace = 'none', scale = 'none', col = 'bluered')
hc = hclust(dist(c))
cut = cutree(hc, k = 2)
table(cut)
ind = which(cut == 2)
finalsample = finalsample[,-ind]

######################################################################

## Creighton SVD, transform to endometrial data

######################################################################
## need to get metagene data from the Creighton et al dataset, so you
## can use it to transform the endometrial data based on the obesity 
## metagene.

library(hgu133a.db) ## load library for microarray annotation

## read series matrix (microarray) data
gse24185exp = read.table('./GSE24185_series_matrix.txt', header=T,
                          sep='',comment.char='!', row.names=1)

## read in the obesity gene file:
obsgenenames = readLines('obsGenes.txt')

## matrix of 799 obesity gene data:
obsgenedat = as.matrix(gse24185exp[obsgenenames,])

## get gene symbols and map it to the data:
genesym2 = mapIds(hgu133a.db, keys = obsgenenames, column = 'SYMBOL', 
                  keytype = 'PROBEID', multiVals = 'first')
rownames(obsgenedat) = as.vector(genesym2)

## remove the NA gene symbols:
tmp = which(is.na(rownames(obsgenedat)))
obsgenedat = obsgenedat[-tmp,]

## check which genes to use for SVD:

## get obesity gene names:
obsgene = as.vector(rownames(obsgenedat))
obsgene = unique(obsgene)

## which of these obesity genes are present in the endometrial cancer
## data?
## Get a list of genes in the endometrial cancer:
tmp = as.vector(rownames(finalsample))
genelist = which(obsgene %in% tmp)
genelist = obsgene[genelist]

## check how many genes we're using:
length(genelist)

## use the genes identified from the endometrial cancer to do SVD on the
## Creighton et al data:
dat = obsgenedat[genelist,]

## normalise the obesity data:
normobsdat = t(apply(dat, 1, function(x) (x-mean(x))/sd(x)))

obssvd = svd(normobsdat) ## do SVD

## make transformation matrix:
transmatrix = diag(1/obssvd$d) %*% t(obssvd$u)

## apply it to the endometrial cancer data:
## make sample by gene matrix of endometrial cancer data:
ucecdata = finalsample[genelist,]

##(normalise ucecdata)
normucecdata = t(apply(ucecdata, 1, function(x) (x-mean(x))/sd(x)))

## for better colour on heatmap
normucecdata[normucecdata > 3] = 3
normucecdata[normucecdata < -3] = -3

ucecmeta = t(transmatrix %*% normucecdata) # apply transformation matrix
ucecmeta = ucecmeta[,1]

## produce heatmap of the metagene with endometrial data:
heatmap.2(normucecdata, scale = 'none', col = 'bluered', trace = 'none', 
          ColSideColors = bluered(length(ucecmeta))[rank(ucecmeta)])

## reorder the heatmap based on the metagene values:
ord = order(ucecmeta)
heatmap.2(normucecdata[,ord], scale = 'none', col = 'bluered', trace = 'none', 
          ColSideColors = bluered(length(ucecmeta))[rank(ucecmeta)][ord],
          Colv = F, Rowv = T)

## plot metagene with clinical data:
plot(ucecbmi[colnames(finalsample),]$bmi,ucecmeta,pch=20, xlab='BMI')
boxplot(ucecmeta~ucecbmi[colnames(finalsample),]$bmi_status)

######################################################################

## DEG of endometrial, compare with Creighton et al obesity genes

######################################################################






