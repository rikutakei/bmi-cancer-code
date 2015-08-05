######################################################################

                        ##R code for task3

######################################################################

## Differentially expressed gene analysis of cesc data.

######################################################################
## pre-analysis shit
setwd('~/Documents/masters/data/task03/')

library(limma)
######################################################################

raw = dget('raw_read_count.txt')
cescbmi = dget('hwdata.txt')

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
x = which(rownames(cescbmi) %in% rownames(raw))
cescbmi = cescbmi[x,]
raw = raw[as.vector(rownames(cescbmi)),]

## any duplicates?
dim(cescbmi) ## 482 samples
dim(raw) ## 482 samples, so no duplicates.

## need to make a model matrix for DEG analysis.
## split the samples into obese or normal/overweight:

x = cescbmi[,4]
x = ifelse(x == 'obese', 'obese','normal/overweight')
design = model.matrix(~x) ## make design matrix
colnames(design)[2] = 'normal/overweight vs. obese'

## since the bmi data IDs are in alphabetical order, order the sample IDs alphabetically in the count data as well:
ord = order(rownames(raw))
raw = raw[ord,]

raw = t(raw) ## transpose the count data for DEG analysis

## DEG analysis:
counts = raw[rowSums(cpm(raw))>9,] ## remove genes with small count data
norm_factor = calcNormFactors(counts)
y = voom(counts, design, lib.size = colSums(counts)*norm_factor) ## voom normalisation

fit = lmFit(y,design)
fit = eBayes(fit)
top = topTable(fit, coef=2, n = nrow(counts)) ## top-ranked DEGs

volcanoplot(fit,coef=2)















