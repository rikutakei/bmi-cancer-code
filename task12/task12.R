######################################################################

                        ##R code for task 12

######################################################################

## Try pathway enrichment analysis again, but with all the samples
## "on the same scale".

######################################################################
## pre-analysis stuff
setwd('~/Documents/masters/data/task12/')

#source the file to get all the functions from task10.R
source('~/Documents/codes/bmi-cancer-code/task10/func10.R')

library(data.table)
library(gplots)

#libraries for data wrangling
library(WGCNA)
library(hgu133a.db)
library(tidyr)
library(dplyr)

#libraries for DEG analysis
library(limma)
library(edgeR)
library(DESeq)
library(affy)

#libraries for pathway analysis
library(org.Hs.eg.db)
library(KEGG.db)
library(GO.db)
library(reactome.db)

##for venn diagram
source('http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/overLapper.R')

######################################################################

## Re-do the Creighton et al analysis, and make sure that the marker
## isn't affected by other clinical variables.

## Read in Cregihton data:
files = readLines('../raw/creighton/files.txt')
files = paste('../raw/creighton/',files,sep='')
crraw = ReadAffy(filenames = files)
crraw = rma(crraw) ## RMA normalise the data
crraw = exprs(crraw) ## change the format into matrix
colnames(crraw) = gsub('.CEL', '', colnames(crraw))

## get clinical data:
crclin = read.csv('../task05/GSE24185_clinical.csv')

## read in the obesity gene file:
crobsgene = readLines('../task05/obsGenes.txt')

######################################################################
## Do DEG analysis on raw data and see if I get the same genes as they did
######################################################################
group = crclin$bmiStatus
group = ifelse(group == 'obese', 'obese', 'non-obese')
design = model.matrix(~group)
fit = lmFit(crraw, design)
fit = eBayes(fit)
tt = topTable(fit, coef = 2, n = nrow(crraw))

## pull out anything with p < 0.01
myobsgene = tt[which(tt$P.Value < 0.01),]
## how many probes from the creighton probes are in the p < 0.01?
length(which(crobsgene %in% rownames(myobsgene)))## 390

myobsgene = rownames(myobsgene)[1:length(crobsgene)] ## pull out top 799 probes

## pull out anything with FC < 1.2
fcgenes = rownames(tt)[which(2^tt$logFC > 1.2)] #only 132 genes

## do metagene analysis with my 799 gene set
mycrmat = crraw[myobsgene,]

mycrmat = t(apply(mycrmat, 1, function(x) (x-mean(x))/sd(x)))
mycrmat[mycrmat > 3] = 3
mycrmat[mycrmat < -3] = -3

heatmap.2(mycrmat,trace='none',scale='none',col='bluered')

## Apply SVD to mycrmat:
svd = svd(mycrmat)
meta = rank(svd$v[,1])/ncol(mycrmat)
ord = order(meta)

heatmap.2(mycrmat[,ord],trace='none',scale='none',col='bluered', ColSideColors = bluered(length(meta))[rank(meta)][ord],Colv=F,Rowv=T)

boxplot(meta~crclin$bmiStatus)
plot(crclin$bmi, meta, pch = 20)

## Do the same thing using Creighton's obesity genes
crobsmat = crraw[crobsgene,]

crobsmat = t(apply(crobsmat, 1, function(x) (x-mean(x))/sd(x)))
crobsmat[crobsmat > 3] = 3
crobsmat[crobsmat < -3] = -3

heatmap.2(crobsmat,trace='none',scale='none',col='bluered')

## Apply SVD to mycrmat:
crsvd = svd(crobsmat)
crmeta = 1-rank(crsvd$v[,1])/ncol(crobsmat)
crord = order(crmeta)

heatmap.2(crobsmat[,crord],trace='none',scale='none',col='bluered', ColSideColors = bluered(length(crmeta))[rank(crmeta)][crord],Colv=F,Rowv=T)

boxplot(crmeta~crclin$bmiStatus)
plot(crclin$bmi, crmeta, pch = 20)

plot(1-crsvd$v[,1], svd$v[,1], pch = 20)

## I wonder if the overlapping genes between cr and my obesity-related genes have similar effect by itself.
redobsgene = myobsgene[which(myobsgene %in% crobsgene)]
test = crraw[redobsgene,]
test = t(apply(test, 1, function(x) (x-mean(x))/sd(x)))
test[test > 3] = 3
test[test < -3] = -3
heatmap.2(test,trace='none',scale='none',col='bluered')

## Apply SVD to mycrmat:
testsvd = svd(test)
testmeta = rank(testsvd$v[,1])/ncol(test)
testord = order(testmeta)

heatmap.2(test[,testord],trace='none',scale='none',col='bluered', ColSideColors = bluered(length(testmeta))[rank(testmeta)][testord],Colv=F,Rowv=T)

boxplot(testmeta~crclin$bmiStatus)
plot(crclin$bmi,testmeta, pch = 20)

######################################################################
## Do something similar, but using only the residuals (remove the effect of other variables)

## Fit linear model on the data using all the other clinical variables:
residuals = crraw
residuals = t(apply(residuals, 1, function(x) lm(x ~ crclin$age + crclin$race + crclin$menopause + crclin$grade + crclin$LNstatus + crclin$ERstatus + crclin$PRstatus + crclin$HER2status)$residuals))

## Fit linear model on the residuals, using obese vs non-obese design matrix
group = crclin$bmiStatus
group = ifelse(group == 'obese', 'obese', 'non-obese')
design = model.matrix(~group)
fit = lmFit(residuals, design)
fit = eBayes(fit)
tt = topTable(fit, coef=2, n = nrow(residuals))

## pull out top 799 genes:
resobsgene = rownames(tt)[which(tt$P.Value < 0.01)]
resobsgene = resobsgene[1:799]

## check which genes are in residual obesity genes compared to my obesity or Creighton obesity genes
length(which( resobsgene %in% myobsgene )) ## 423
length(which( resobsgene %in% crobsgene )) ## 168
length(which(resobsgene %in% redobsgene)) ## 112

## make venn diagram
setlist = list(residual = resobsgene, my_genes = myobsgene, creighton = crobsgene)
OLlist = overLapper(setlist=setlist, sep='_', type='vennsets')
counts = sapply(OLlist$Venn_List, length)
vennPlot(counts = counts)

## are all probes in crobsgene present in the toptable of residuals?
length(which(crobsgene %in% rownames(tt)[which(tt$P.Value < 0.01)])) ## 216

##store the  216 genes in a variable
rescrobsgene = rownames(tt)[which(tt$P.Value < 0.01)]
rescrobsgene = rescrobsgene[which(rescrobsgene %in% crobsgene)]

setlist = list(residual = resobsgene, my_genes = myobsgene, creighton = crobsgene, cr_residual = rescrobsgene)
OLlist = overLapper(setlist=setlist, sep='_', type='vennsets')
counts = sapply(OLlist$Venn_List, length)
vennPlot(counts = counts)

##resobsgene:
resobsmat = crraw[resobsgene,]

resobsmat = t(apply(resobsmat, 1, function(x) (x-mean(x))/sd(x)))
resobsmat[resobsmat > 3] = 3
resobsmat[resobsmat < -3] = -3

heatmap.2(resobsmat,trace='none',scale='none',col='bluered')

## Apply SVD to myresmat:
ressvd = svd(resobsmat)
resmeta = rank(ressvd$v[,1])/ncol(resobsmat)
resord = order(resmeta)

heatmap.2(resobsmat[,resord],trace='none',scale='none',col='bluered', ColSideColors = bluered(length(resmeta))[rank(resmeta)][resord],Colv=F,Rowv=T)

boxplot(resmeta~crclin$bmiStatus)
plot(crclin$bmi, resmeta, pch = 20)

##rescrobsgene:
rescrobsmat = crraw[rescrobsgene,]

rescrobsmat = t(apply(rescrobsmat, 1, function(x) (x-mean(x))/sd(x)))
rescrobsmat[rescrobsmat > 3] = 3
rescrobsmat[rescrobsmat < -3] = -3

heatmap.2(rescrobsmat,trace='none',scale='none',col='bluered')

## Apply SVD to myrescrmat:
rescrsvd = svd(rescrobsmat)
rescrmeta = rank(rescrsvd$v[,1])/ncol(rescrobsmat)
rescrord = order(rescrmeta)

heatmap.2(rescrobsmat[,rescrord],trace='none',scale='none',col='bluered', ColSideColors = bluered(length(rescrmeta))[rank(rescrmeta)][rescrord],Colv=F,Rowv=T)

boxplot(rescrmeta~crclin$bmiStatus)
plot(crclin$bmi, rescrmeta, pch = 20)

##resredobsgene:
resredobsgene = rownames(tt)[which(tt$P.Value < 0.01)]
resredobsgene = resredobsgene[which(resredobsgene %in% redobsgene)]

resredobsmat = crraw[resredobsgene,]

resredobsmat = t(apply(resredobsmat, 1, function(x) (x-mean(x))/sd(x)))
resredobsmat[resredobsmat > 3] = 3
resredobsmat[resredobsmat < -3] = -3

heatmap.2(resredobsmat,trace='none',scale='none',col='bluered')

## Apply SVD to myresredmat:
resredsvd = svd(resredobsmat)
resredmeta = rank(resredsvd$v[,1])/ncol(resredobsmat)
resredord = order(resredmeta)

heatmap.2(resredobsmat[,resredord],trace='none',scale='none',col='bluered', ColSideColors = bluered(length(resredmeta))[rank(resredmeta)][resredord],Colv=F,Rowv=T)

boxplot(resredmeta~crclin$bmiStatus)
plot(crclin$bmi, resredmeta, pch = 20)

######################################################################
## repeat the analyses above, but only using the caucasian samples.
######################################################################

## pull out relative information for caucasian samples:
camat = crraw
ind = which(crclin$race == 'w')
camat = camat[,ind] ## 77 caucasian samples in total

caclin = crclin[ind,]

group = caclin$bmiStatus
group = ifelse(group == 'obese', 'obese', 'non-obese')
design = model.matrix(~group)
fit = lmFit(camat, design)
fit = eBayes(fit)
tt = topTable(fit, coef = 2, n = nrow(camat))

## pull out anything with p < 0.01
caobsgene = tt[which(tt$P.Value < 0.01),]
## how many probes from the creighton probes are in the p < 0.01?
cacrobsgene = crobsgene[which(crobsgene %in% rownames(caobsgene))]## 271

caobsgene = rownames(caobsgene)[1:length(crobsgene)] ## pull out top 799 probes

## do metagene analysis with 799 genes found from caucasian samples
caobsmat = crraw[caobsgene,]

caobsmat = t(apply(caobsmat, 1, function(x) (x-mean(x))/sd(x)))
caobsmat[caobsmat > 3] = 3
caobsmat[caobsmat < -3] = -3

heatmap.2(caobsmat,trace='none',scale='none',col='bluered')

## Apply SVD to caobsmat:
casvd = svd(caobsmat)
cameta = rank(casvd$v[,1])/ncol(caobsmat)
caord = order(cameta)

heatmap.2(caobsmat[,caord],trace='none',scale='none',col='bluered', ColSideColors = bluered(length(cameta))[rank(cameta)][caord],Colv=F,Rowv=T)

boxplot(cameta~crclin$bmiStatus)
plot(crclin$bmi, cameta, pch = 20)

## do metagene analysis with  271 genes found from caucasian samples that were also in the Creighton gene list
cacrobsmat = crraw[cacrobsgene,]

cacrobsmat = t(apply(cacrobsmat, 1, function(x) (x-mean(x))/sd(x)))
cacrobsmat[cacrobsmat > 3] = 3
cacrobsmat[cacrobsmat < -3] = -3

heatmap.2(cacrobsmat,trace='none',scale='none',col='bluered')

## Apply SVD to cacrobsmat:
cacrsvd = svd(cacrobsmat)
cacrmeta = rank(cacrsvd$v[,1])/ncol(cacrobsmat)
cacrord = order(cacrmeta)

heatmap.2(cacrobsmat[,cacrord],trace='none',scale='none',col='bluered', ColSideColors = bluered(length(cacrmeta))[rank(cacrmeta)][cacrord],Colv=F,Rowv=T)

boxplot(cacrmeta~crclin$bmiStatus)
plot(crclin$bmi, cacrmeta, pch = 20)

## repeat on residual data
camatres = camat
camatres = t(apply(camatres, 1, function(x) lm(x ~ caclin$age + caclin$menopause + caclin$grade + caclin$LNstatus + caclin$ERstatus + caclin$PRstatus + caclin$HER2status)$residuals))

group = caclin$bmiStatus
group = ifelse(group == 'obese', 'obese', 'non-obese')
design = model.matrix(~group)
fit = lmFit(camatres, design)
fit = eBayes(fit)
tt = topTable(fit, coef = 2, n = nrow(camatres))

## pull out anything with p < 0.01
caresobsgene = tt[which(tt$P.Value < 0.01),]
## how many probes from the creighton probes are in the p < 0.01?
cacrresobsgene = crobsgene[which(crobsgene %in% rownames(caresobsgene))]## 166

caresobsgene = rownames(caresobsgene)[1:length(crobsgene)] ## pull out top 799 probes

# How many genes from non-residual caucasian sample data are also in the residual data
which(caobsgene %in% caresobsgene) %>% length #0

## do metagene analysis with 799 genes found from caucasian samples
caresobsmat = crraw[caresobsgene,]

caresobsmat = t(apply(caresobsmat, 1, function(x) (x-mean(x))/sd(x)))
caresobsmat[caresobsmat > 3] = 3
caresobsmat[caresobsmat < -3] = -3

heatmap.2(caresobsmat,trace='none',scale='none',col='bluered')

## Apply SVD to caresobsmat:
caressvd = svd(caresresobsmat)
caresmeta = rank(caressvd$v[,1])/ncol(caresresobsmat)
caresord = order(caresmeta)

heatmap.2(caresobsmat[,caresord],trace='none',scale='none',col='bluered', ColSideColors = bluered(length(caresmeta))[rank(caresmeta)][caresord],Colv=F,Rowv=T)

boxplot(caresmeta~crclin$bmiStatus)
plot(crclin$bmi, caresmeta, pch = 20)

######################################################################
## metagene analysis:
######################################################################

#load all the data:
files = readLines('../task09/files.txt')
files = paste('../task09/',files,sep='')
bmifiles = gsub('mat','bmi',files)

#create variables for cancer data matrix
for (i in 1:length(files)) {
    txt = gsub('t9.txt','',files[i])
    txt = gsub('../task09/','',txt)
    assign(txt, t(dget(files[i]))) #transpose the matrix so the columns are samples
    files[i] = txt #rename files to its variable names

    txt = gsub('t9.txt','',bmifiles[i])
    txt = gsub('../task09/','',txt)
    assign(txt, dget(bmifiles[i]))
    bmifiles[i] = txt #rename files to its variable names
}

# make transformation matrices:
#restrmat    = diag(1/ressvd$d)%*%t(ressvd$u)
#resredtrmat = diag(1/resredsvd$d)%*%t(resredsvd$u)
#rescrtrmat  = diag(1/rescrsvd$d)%*%t(rescrsvd$u)
#catrmat     = diag(1/casvd$d)%*%t(casvd$u)
#cacrtrmat   = diag(1/cacrsvd$d)%*%t(cacrsvd$u)

prefix = c("res","resred","rescr","ca","cacr","cares","cacrres")

## need to re-adjust the transformation matrices so that it's compatible with the cancer data
for (j in 1:length(prefix)){
    ## get relevant genes:
    genes = paste(prefix[j], "obsgene",sep = '')
    genes = get(genes)
    genesymbol = mapIds(hgu133a.db, keys = genes, column = 'SYMBOL', keytype = "PROBEID", multiVals = 'first') ## convert probe IDs into symbols
    ## account for any NA genes:
    if(length(which(is.na(genesymbol))) > 0) {
        genesymbol = genesymbol[-(which(is.na(genesymbol)))] #remove any NAs
    }
    genesymbol = unique(genesymbol) #include only unique genes.
    genesymbol = genesymbol[which(genesymbol %in% rownames(BLCAmat))]
    genes = mapIds(hgu133a.db, keys = genesymbol, column = 'PROBEID', keytype = "SYMBOL", multiVals = 'first') ## convert symbols into probe ids

    ## make corresponding metagene
    tmp = crraw[genes,]
    tmp = t(apply(tmp, 1, function(x) (x-mean(x))/sd(x)))
    tmp[tmp > 3] = 3
    tmp[tmp < -3] = -3

    ## Apply SVD to tmp and make transformation matrix:
    tmpsvd   = svd(tmp)
    tmptrmat = diag(1/tmpsvd$d)%*%t(tmpsvd$u)

    txt = paste(prefix[j],"tcga-res.pdf",sep='')

    pdf(file = txt)
    for(i in 1:length(files)) {
        dat = get(files[i]) ## get data
        dat = dat[genesymbol,] ## pull out relevant genes

        dat = log10(dat + 1)
        dat = t(apply(dat, 1, function(x) (x-mean(x))/sd(x)))
        dat[dat < -3] = -3
        dat[dat > 3] = 3

        ## account for the NaN values (change it to 0):
        if(length(which(is.nan(dat))) > 0) {
            k = arrayInd(which(is.nan(dat)), dim(dat))
            dat[k] = 0
        }

        tmpmeta = t(tmptrmat %*% dat)
        tmpmeta = tmpmeta[,1]
        tmpmeta = rank(tmpmeta)/ncol(dat)
        tmpord = order(tmpmeta)

        txt = paste(files[i], 'with')
        txt = paste( txt, prefix[j])

        ## make results
        heatmap.2(dat,trace='none',scale='none',col='bluered', main = txt)
        heatmap.2(dat[,tmpord],trace='none',scale='none',col='bluered', ColSideColors = bluered(length(tmpmeta))[rank(tmpmeta)][tmpord],Colv=F,Rowv=T, main = txt)

        bmi = get(gsub('mat', 'bmi', files[i]))

        boxplot(tmpmeta~bmi[,4], main = txt)
        plot(bmi[,3], tmpmeta, pch = 20, main = txt)
    }
    dev.off()
}











