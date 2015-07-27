######################################################################

                        ##R code for task2

######################################################################

## set working directory:
setwd('~/Documents/masters/data/task02')

## read the RNA seq data:
seq = read.table('exp_seq.UCEC-US.tsv', sep = '\t', header = T)

## the file format is useless, so use dplyr to manipulate it:
library(dplyr)
library(gplots)

seq = tbl_df(seq) ## change data frame into tbl class

## extract all the raw count data and put it in a matrix:
initsamplenames = filter(seq, gene_id == 'PAN2')$submitted_sample_id
initgenenames = unique(seq$gene_id)
raw = matrix(nrow = 541, ncol = 20502)
rownames(raw) = initsamplenames
colnames(raw) = initgenenames

for (i in 1:ncol(raw)) {
    gene = as.vector(initgenenames[i])
    tmpreadcount = filter(seq, gene_id == gene) %>%
    select(submitted_sample_id, raw_read_count)

    for (j in 1:nrow(tmpreadcount)) {
        sample = as.matrix(tmpreadcount[j,1])
        raw[sample, gene] = as.numeric(tmpreadcount[j, 2])
    }
    print(paste("Finished processing ", gene))
}

## save the matrix:
## dput(raw, file = 'raw_read_count.txt')

## need to SVD the raw read count and validate the Creighton et al results

## load obesity gene names and series matrix data
library(hgu133a.db) ## load library for microarray annotation

## read series matrix (microarray) data
## gse24185exp = read.table('./GSE24185_series_matrix.txt', header=T,
##                          sep='',comment.char='!', row.names=1)

## read in the obesity gene file:
obsgenenames = readLines('obsGenes.txt')

## check the valid columns to get from the hgu133a.db
columns(hgu133a.db)

## select the gene symbols that are assigned to the probes:
genesym = select(hgu133a.db, keys = obsgenenames, columns = 'SYMBOL')

## OR map the probe IDs to the symbols:
genesym2 = mapIds(hgu133a.db, keys = obsgenenames, column = 'SYMBOL', 
                  keytype = 'PROBEID', multiVals = 'first')

## NOTE: this only takes the first gene symbol that the probe ID is
## matched to, if the ID has multiple gene values.

## going to use the genesym2 data from here:

## get unique gene symbols and remove NA:
genesym2 = unique(genesym2)[2:length(genesym2)]
genesym2 = sort(genesym2) ## sort the symbols

## find which columns of the raw data I need:
colind = which(colnames(raw) %in% genesym2)

## extract the columns:
obsgeneraw = raw[,colind]

## there is a row of NAs in the data, so delete that:
which(is.na(obsgeneraw[,1] == T)) ## to find which row has the NAs
obsgeneraw = obsgeneraw[-221,] ## delete the row of NAs

## need to find which samples are actually from the primary tumour site
## load specimen data:
specdata = read.table(file = 'specimen.UCEC-US.tsv', sep = '\t',
                      header = T)
primtum = specdata[specdata$specimen_type == 'Primary tumour - solid tissue',]
primtumid = as.vector(primtum$submitted_specimen_id)

## scale the data so it has mean of 0, standard deviation of 1:
scaledobsgene = scale(obsgeneraw)
scaledobsgene = scaledobsgene[order(rownames(scaledobsgene)),] ## order by rownames

## make the max and min values to be 3 or -3:
scaledobsgene[scaledobsgene >= 3] = 3
scaledobsgene[scaledobsgene <= -3] = -3

## change the rownames so that it matches the TCGA sample names
specimenid = gsub(pattern = '-[[:alnum:]]{3}-[[:alnum:]]{4}-[[:alnum:]]{2}',
                               replacement = '', rownames(scaledobsgene))

## find which of these samples are from the primary tumour site:
ind = which(!(primtumid %in% specimenid))
scaledobsgene = scaledobsgene[-ind,]

## Import the clinical data to match the scaledobsgene.
clin = read.table('ucec_clinical_patient.txt', sep = '\t', skip = 1,
                  header = T)
clin = clin[-1,] ## remove random header

## get height and weight data:
hwdata = clin[,c('weight','height')]
rownames(hwdata) = clin$bcr_patient_barcode ## get rownames

## remove any samples with no height or weight data:
hwdata = hwdata[-(which(hwdata == '[Not Available]')),]
hwdata = hwdata[-(which(hwdata[,2] == '[Not Available]')),]

## convert to numerical matrix
rows = rownames(hwdata)
hwdata = as.matrix((hwdata))
hwdata = apply(hwdata, 2, as.numeric)
rownames(hwdata) = rows

## found a sample with a height of 66 cm (probs should be 166cm)
hwdata[hwdata[,2] < 100, 2] = 100 + hwdata[hwdata[,2] < 100, 2]

## calculate BMI:
bmi = hwdata[,1]/((hwdata[,2] / 100) ^ 2) ## need to convert height into m
bmi = round(bmi, 2)
hwdata = cbind(hwdata, bmi)

 #check which samples have BMI data:
#bmidataind = which(!(rownames(scaledobsgene) %in% rownames(hwdata)))
#scaledobsgene = scaledobsgene[-bmidataind,]

 #select the unique BMI samples:
#hwdata = hwdata[-(which(!(rownames(hwdata) %in% rownames(scaledobsgene)))),]


## need to identify which sample names are repeated:
tcgaid = gsub(pattern = '-[[:alnum:]]{3}-[[:alnum:]]{3}-[[:alnum:]]{4}-[[:alnum:]]{2}',
                               replacement = '', rownames(scaledobsgene))
rownames(scaledobsgene) = tcgaid

tmprow = rownames(scaledobsgene) ## get rownames
repind = table(tmprow) > 1 ## find any names that are repeated
repnames = unique(tmprow)[repind] ## get the sample names that are repeated

## for loop to graph the replicates to decide which one should be
## discarded:
pdf(file = 'replicates.pdf',width = 8) ## output the heatmaps

for (i in 1:length(repnames)) {
    ind = which(rownames(scaledobsgene) == repnames[i])
    heatmap.2(scaledobsgene[ind,], trace = 'none', col = 'bluered',
              scale = 'none', labRow = ind, cexRow = 0.9)
}

dev.off()

## make a vector of the samples to be used in SVD

## These are the samples I've chosen to keep:
tmp = c(47,64,68,124,140,142,152,155,166,168,170,172,286,290,292,295,
        297,299,302,305,307,310,414,416)

tmp2 = vector()
for (i in 1:length(repnames)) {
    ind = which(rownames(scaledobsgene) == repnames[i])
    tmp2 = append(tmp2, ind)
}

tmp3 = tmp2[which(tmp2 %in% tmp == F)] ## the unwanted sample row index

## remove the duplicated samples:
scaledobsgene = scaledobsgene[-tmp3,]

## find all the samples that have BMI data:
tmp = which(rownames(hwdata) %in% rownames(scaledobsgene))
tmp = rownames(hwdata)[tmp]

scaledobsgene = scaledobsgene[tmp,]
hwdata = hwdata[tmp,]

## make colours for BMI:
obsCol = greenred(nrow(scaledobsgene))[rank(hwdata[,3])]

## correlation distance function:
distfun = function(x) { as.dist(1-cor(t(x))) }

## make heatmap using the obesity data:
heatmap.2(t(scaledobsgene), trace = 'none', col = 'bluered',
          scale = 'none', ColSideColors = obsCol, distfun = distfun)

## make heatmap, ordered by BMI:
ord = order(hwdata[,3])
heatmap.2(t(scaledobsgene)[,ord], trace = 'none', col = 'bluered',
          scale = 'none', ColSideColors = obsCol[ord], distfun = distfun,
          Colv = 'none')

## make BMI metagene from the expression data:
## scaled BMI metagene:
scaledsvd = svd(t(scaledobsgene))
scaledmetagene = rank(scaledsvd$v[,1])/nrow(scaledobsgene)

## Unscaled (non-standardised) BMI metagene:
obsdat = obsgeneraw
tcgaid = gsub(pattern = '-[[:alnum:]]{3}-[[:alnum:]]{3}-[[:alnum:]]{4}-[[:alnum:]]{2}',
                               replacement = '', rownames(obsdat))
rownames(obsdat) = tcgaid
obsdat = obsdat[tmp,]

unscaledsvd = svd(t(obsdat))
unscaledmetagene = -unscaledsvd$v[,1]

## compare metagene:
plot(scaledmetagene,unscaledmetagene)

## add metagene to the heatmap:
heatmap.2(t(scaledobsgene), trace = 'none', col = 'bluered',
          scale = 'none', ColSideColors = bluered(length(scaledmetagene))[rank(scaledmetagene)])

## re-order heatmap based on metagene values:
ordmg = order(scaledmetagene)
heatmap.2(t(scaledobsgene)[,ordmg], trace = 'none', col = 'bluered',
          scale = 'none', Colv = F, Rowv = T,
          ColSideColors = bluered(length(scaledmetagene))[rank(scaledmetagene)][ordmg])


######################################################################

## Quality assessment

## make density plots for all the samples on log scale:

rawlog = log2(raw)
d = density(rawlog[1,])
plot(d, ylim = c(0,0.2))

for (i in 2:nrow(rawlog)) {
    d = density(rawlog[i,])
    lines(d)
}


## make correlation matrix for samples across all genes:
c = cor(t(raw), use='all.obs',method='pearson')

heatmap.2(c, trace = 'none', col = 'bluered',scale = 'none')

## get the dendrogram on the heatmap
hc = hclust(dist(c))
plot(hc, cex = 0.1, hang = -1) #show dendrogram

tmpcut = cutree(hc, k = 2) #cut dendrogram into two groups
table(tmpcut) # check which group to remove

ind = which(tmpcut == 2)
finalsample = raw[-ind,] # remove the bad samples

## repeat the above until no bad samples (correlation < 0.6)
## takes about 3 cycles (about 312 samples in the end)

## check for any replicates:
tcgaid = gsub(pattern = '-[[:alnum:]]{3}-[[:alnum:]]{3}-[[:alnum:]]{4}-[[:alnum:]]{2}',
                               replacement = '', rownames(finalsample))

rownames(finalsample) = tcgaid

## scatterplot stuff (?):
plot(x = finalsample[4,], y = finalsample[101,], pch = '.')


## voom normalisation:

library(edgeR)
library(limma)

## need to remove low count data:
count = t(finalsample)
count = count[rowSums(cpm(count)) > 9,]

dge = DGEList(counts = count)
dge = calcNormFactors(dge)


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
tmp = as.vector(rownames(count))
genelist = which(obsgene %in% tmp)
genelist = obsgene[genelist]

## use the genes identified from the endometrial cancer to do SVD on the
## Creighton et al data:
dat = obsgenedat[genelist,]

## normalise the obesity data:
normobsdat = t(apply(dat, 1, function(x) (x-mean(x))/sd(x)))

obssvd = svd(normobsdat) ## do SVD

## make transformation matrix:
transmatrix = diag(1/obssvd$d) %*% t(obssvd$u)

## apply it to the endometrial cancer data:
cescdata = count[genelist,]

##(normalise cescdata??)
normcescdata = t(apply(cescdata, 1, function(x) (x-mean(x))/sd(x)))

## for better colour on heatmap
normcescdata[normcescdata > 3] = 3
normcescdata[normcescdata < -3] = -3

cescmeta = t(transmatrix %*% normcescdata) # apply transformation matrix
cescmeta = cescmeta[,1]

## produce heatmap of the metagene with endometrial data:
heatmap.2(normcescdata, scale = 'none', col = 'bluered', trace = 'none', 
          ColSideColors = bluered(length(cescmeta))[rank(cescmeta)])

## reorder the heatmap based on the metagene values:
ord = order(cescmeta)
heatmap.2(normcescdata[,ord], scale = 'none', col = 'bluered', trace = 'none', 
          ColSideColors = bluered(length(cescmeta))[rank(cescmeta)][ord],
          Colv = F, Rowv = T)

## Need to check whether the metagene corresponds to BMI/obesity:
## get BMI data for the samples used in the cescmeta data
samplenames = names(cescmeta)
ind = which(rownames(hwdata) %in% samplenames)
cescbmi = as.data.frame(hwdata[ind,3])

bmimeta = cescmeta[rownames(cescbmi)]

tmp = cescbmi[,1]

for(i in 1:length(tmp)){
    if (cescbmi[i,1] <= 25) {
        tmp[i] = 'normal'
    } else if ((cescbmi[i,1] > 25) && cescbmi[i,1] <= 30) {
        tmp[i] = 'obese'
    } else {
        tmp[i] = 'overweight'
    }
}

cescbmi = cbind(cescbmi, tmp)
names(cescbmi) = c('bmi_value', 'bmi_status')




