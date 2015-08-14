######################################################################

                        ##R code for task3

######################################################################

## Differentially expressed gene analysis of ucec data.

######################################################################
## pre-analysis stuff
setwd('~/Documents/masters/data/task03/')

library(limma)
library(edgeR)
library(DESeq)
######################################################################

raw = dget('raw_read_count.txt')
ucecbmi = dget('hwdata.txt')

## need to identify which samples are present in both the hwdata and the count data:

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

## need to make a model matrix for DEG analysis.
## Quality control - remove 'bad' smaples:
finalsample = t(raw)
## repeat this a few times
c = cor(finalsample, use = 'all.obs', method = 'pearson')
heatmap.2(c, trace = 'none', scale = 'none', col = 'bluered')
hc = hclust(dist(c))
cut = cutree(hc, k = 2)
table(cut)
ind = which(cut == 2)
finalsample = finalsample[,-ind]

## get the BMI data for these samples:
samplenames = as.vector(colnames(finalsample))
finalsamplebmi = ucecbmi[samplenames,]

## split the samples into obese or normal/overweight:
x = finalsamplebmi[,4]
x = ifelse(x == 'obese', 'obese','normal/overweight')
design = model.matrix(~x) ## make design matrix
colnames(design)[2] = 'normal/overweight vs. obese'

## since the bmi data IDs are in alphabetical order, order the sample IDs alphabetically in the count data as well:
ord = order(colnames(finalsample))
finalsample = finalsample[,ord]

#finalsample = t(finalsample) ## transpose the count data for DEG analysis

## DEG analysis:
counts = finalsample[rowSums(cpm(finalsample))>9,] ## remove genes with small count data
norm_factor = calcNormFactors(counts)
y = voom(counts, design, lib.size = colSums(counts)*norm_factor) ## voom normalisation

fit = lmFit(y,design)
fit = eBayes(fit)
top = topTable(fit, coef=2, n = nrow(counts)) ## top-ranked DEGs

#########################################################################
## Mik's stuff

top[1:5,]
ord<-order(design[,2])
heatmap.2(y$E[match(rownames(top)[1:100],rownames(y)),],trace='none',
          ColSideColors=ifelse(design[,2]==1,"black","blue"))

zz<-t(apply(y$E[match(rownames(top)[1:5],rownames(y)),],1,
            function(x) (x-mean(x))/sd(x)))
zz[zz< -3]<- -3
zz[zz>3] <- 3
heatmap.2(zz[,ord],trace='none',
          ColSideColors=ifelse(design[,2]==1,"black","blue")[ord],
          Colv=F,scale='row')
mean(y$E[match(rownames(top)[1],rownames(y)),design[,2]==1])
mean(y$E[match(rownames(top)[1],rownames(y)),design[,2]==0])

top[1,]

#########################################################################

volcanoplot(fit,coef=2)

## get DEG from top table:
ucecdeg = top[top$adj.P.Val < 0.05,]

## Store the DEG identified (218 genes):
degnames = as.vector(rownames(ucecdeg))

## before doing singular value decomposition, remove samples with correlation less than 0.6
finalsample = raw

## repeat this a few times
c = cor(finalsample, use = 'all.obs', method = 'pearson')
heatmap.2(c, trace = 'none', scale = 'none', col = 'bluered')
hc = hclust(dist(c))
cut = cutree(hc, k = 2)
table(cut)
ind = which(cut == 2)
finalsample = finalsample[,-ind]

## SVD:

## make metagene using the DEGs identified above.
degmat = finalsample[degnames,]

normdegmat = t(apply(degmat,1,function(x) (x-mean(x))/sd(x)))
normdegmat[normdegmat > 3] = 3
normdegmat[normdegmat < -3] = -3

ucecmeta = rank(svd(normdegmat)$v[,1])/ncol(degmat)
ord = order(ucecmeta)

heatmap.2(normdegmat[,ord], trace = 'none', col = 'bluered', scale = 'none',
          ColSideColors = bluered(length(ucecmeta))[rank(ucecmeta)][ord],
          Colv = F, Rowv = T)

#########################################################################
## mik's stuff - order by BMI rather than metagene:
mik<-ucecbmi$bmi[match(colnames(normdegmat),rownames(ucecbmi))]
ord = order(mik)
mikCol<-bluered(length(mik))[rank(mik)]
heatmap.2(normdegmat[,ord], trace = 'none', col = 'bluered', scale = 'none',
          ColSideColors = mikCol[ord],
          Colv = F, Rowv = T)

#########################################################################














