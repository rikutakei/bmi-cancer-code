###############################################################################

                        ## Final R code for the project

###############################################################################

## set working directory:
setwd('/home/riku/Documents/masters/data/task15')

## Load previously saved workspace:
# load(".RData")

## load functions:
source('/home/riku/Documents/codes/bmi-cancer-code/task15/functions_final.R')

## load libraries:
# library(RColorBrewer)
# library(mclust)
library(GMD)
library(colorRamps)
library(data.table)
library(gplots)
library(lattice)
library(sva)

## heatmap stuff
library(devtools)
devtools::install_github("TomKellyGenetics/heatmap.2x", ref="master")
library(heatmap.2x)

## libraries for data wrangling
library(WGCNA)
library(dplyr)
library(hgu133a.db)
library(tidyr)

## libraries for DEG analysis
# library(edgeR)
library(DESeq)
library(affy)
library(limma)

## libraries for pathway analysis
library(GO.db)
library(KEGG.db)
library(org.Hs.eg.db)
library(reactome.db)

## script to make venn diagrams
source('http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/My_R_Scripts/overLapper.R')

###############################################################################
## LOAD DATA

# source("load_data.R")

# Load Creighton et al data:

files = readLines('./raw/creighton/files.txt')
files = paste('./raw/creighton/', files, sep='')
crraw = ReadAffy(filenames = files)

crrma = rma(crraw) ## RMA normalise the data
crrma = exprs(crrma) ## change the format into matrix
colnames(crrma) = gsub('.CEL', '', colnames(crrma))

crmas = mas5(crraw) ## Mas5 normalise the data
crmas = exprs(crmas) ## change the format into matrix
crmas = log2(crmas) ## log2 the data
colnames(crmas) = gsub('.CEL', '', colnames(crmas))

# standardise data (mean = 0, sd = 1):
crstdrma = t(apply(crrma, 1, function(x) (x-mean(x))/sd(x)))
crstdmas = t(apply(crmas, 1, function(x) (x-mean(x))/sd(x)))

# load Creighton's obesity-associated genes:
crobsgenes = read.csv('./obsgenes/crobsgenes.txt', header=F)
crobsgenes = as.vector(crobsgenes[,1])

# load clinical data:
crclin = read.csv('./clindata/crclin.csv', sep=',', header=T)
rownames(crclin) = crclin$geo_accession

###############################################################################
## Fuentes-Mattei et al data:

# Load Fuentes-Mattei et al data:
files = readLines('./raw/fuentes-mattei/files.txt')
files = paste('./raw/fuentes-mattei/', files, sep='')
fmraw = ReadAffy(filenames = files)

fmrma = rma(fmraw) ## RMA normalise the data
fmrma = exprs(fmrma) ## change the format into matrix
colnames(fmrma) = gsub('GSM[^_]*_', '', colnames(fmrma))
colnames(fmrma) = gsub('.CEL', '', colnames(fmrma))

fmmas = mas5(fmraw) ## Mas5 normalise the data
fmmas = exprs(fmmas) ## change the format into matrix
fmmas = log2(fmmas) ## log2 the data
colnames(fmmas) = gsub('GSM[^_]*_', '', colnames(fmmas))
colnames(fmmas) = gsub('.CEL', '', colnames(fmmas))

# standardise data (mean = 0, sd = 1):
fmstdrma = t(apply(fmrma, 1, function(x) (x-mean(x))/sd(x)))
fmstdmas = t(apply(fmmas, 1, function(x) (x-mean(x))/sd(x)))

# load Fuentes-Mattei's obesity-associated genes:
fmobsgenes = read.csv('./obsgenes/fmobsgenes.txt', header=T)

# load clinical data:
fmclin = read.csv('./clindata/fmclin.csv', sep=',', header=T)
rownames(fmclin) = gsub('[ ()]', '_', fmclin$Sample.name)

###############################################################################
## Cris Print's Breast cancer data:

# Load Print et al data:
files = readLines('./raw/cris/files.txt')
files = paste('./raw/cris/', files, sep='')
crisraw = ReadAffy(filenames = files)

crisrma = rma(crisraw) ## RMA normalise the data
crisrma = exprs(crisrma) ## change the format into matrix
colnames(crisrma) = gsub('GSM[^_]*_', '', colnames(crisrma))

crismas = mas5(crisraw) ## Mas5 normalise the data
crismas = exprs(crismas) ## change the format into matrix
crismas = log2(crismas) ## log2 the data
colnames(crisrma) = gsub('GSM[^_]*_', '', colnames(crisrma))

# Load clinical data:
crisclin = read.csv('./clindata/crisclin2.csv', sep=',', header=T)
rownames(crisclin) = crisclin$CelFile

# Remove samples with no BMI information:
ind = which(!is.na(crisclin$BMI))
crisclin = crisclin[ind,]
crisrma = crisrma[,ind]
crismas = crismas[,ind]

# Standardise data (mean = 0, sd = 1):
crisstdrma = t(apply(crisrma, 1, function(x) (x-mean(x))/sd(x)))
crisstdmas = t(apply(crismas, 1, function(x) (x-mean(x))/sd(x)))

# Clean the clinical data so that it's compatible with the functions/other
# parts of code:
colnames(crisclin)[which(names(crisclin) == 'BMI')] = 'bmi'
bmiStatus = crisclin$bmi
bmiStatus[bmiStatus >= 30] = 'obese'
bmiStatus[bmiStatus < 25] = 'normal'
bmiStatus[bmiStatus < 30] = 'overweight'
crisclin = cbind(crisclin, bmiStatus = bmiStatus)

###############################################################################
## Gatza et al. cancer data:

# Import the sample names included in different GSE experiments:
files = readLines('./raw/gatza/sample_names/files.txt')
files = paste('./raw/gatza/sample_names/', files, sep='')

sample_list = list()
for (i in 1:length(files)) {
	tmp = readLines(files[i])
	txt = gsub('.*/', '', files[i])
	txt = gsub('_.*', '', txt)
	sample_list[[i]] = tmp
	names(sample_list)[i] = txt
}
sample_list = lapply(sample_list, function(x) gsub('.gz','', x))

# Import the hgu133a samples:
files = readLines('./raw/gatza/hgu133a/files.txt')
sample_list = lapply(sample_list, function(x) x[which(x %in% files)])

# remove GSE dataset with no HGU133A dataset:
sample_list = sample_list[lapply(sample_list, length) > 0]
sample_list = lapply(sample_list, function(x) paste('./raw/gatza/hgu133a/', x, sep=''))

# import each dataset separately
rawlist = sample_list
rawlist = lapply(rawlist, function(x) ReadAffy(filenames = x))

# RMA normalise each data set independently:
rmalist = lapply(rawlist, function(x) rma(x))
rmalist = lapply(rmalist, function(x) exprs(x))

# MAS5 normalise each data set independently:
mas5list = lapply(rawlist, function(x) mas5(x))
mas5list = lapply(mas5list, function(x) exprs(x))
mas5list = lapply(mas5list, function(x) log2(x))

# combine the independently normalised data into a single matrix:
gatzarma = do.call(cbind, rmalist)
gatzamas = do.call(cbind, mas5list)

dim(gatzarma) # 22283 by 1060
dim(gatzamas) # 22283 by 1060

# generate batch info:
batch = sample_list
batch = lapply(batch, function(x) gsub('./raw/gatza/hgu133a/', '' , x))
batch = as.vector(unlist(batch))
groups = c()
for (i in 1:length(sample_list)) {
	tmp = rep(names(sample_list)[[i]], length(sample_list[[i]]))
	groups = c(groups, tmp)
}
batch = cbind(batch, groups)

# Correct for the batch effect:
# First make a model matrix for the different batches:
mod = model.matrix(~1, data=as.data.frame(batch))

# Adjust for batch effect:
gatzarma = ComBat(dat = gatzarma, batch=batch[,2], mod=mod)
gatzamas = ComBat(dat = gatzamas, batch=batch[,2], mod=mod)

# standardise the data:
gtrmastd = t(apply(gatzarma, 1, function(x) (x-mean(x))/sd(x)))
gtmasstd = t(apply(gatzamas, 1, function(x) (x-mean(x))/sd(x)))

## import pathway gene list from the Gatza paper
files = readLines('./gatzagenelist/pathlist.txt')
files = paste('./gatzagenelist/', files, sep='')
for (i in 1:length(files)) {
	txt = gsub('./gatzagenelist/','',files[i])
	txt = gsub('.txt','',txt)
	genes = readLines(files[i])
	genes = genes[which(genes %in% rownames(gtrmastd))]
	assign(txt, genes)
}

# make a variable with list of the name of the pathways
paths = gsub('./gatzagenelist/', '', files)
paths = gsub('.txt', '', paths)

###############################################################################
## ICGC data:

## (may have to use the server for extra RAM (use the command: ulimit -s 20480))

files = c('BLCA', 'CESC', 'COAD', 'KIRP', 'LIHC', 'READ', 'SKCM', 'UCEC')
cancerfiles = paste('./raw/raw-cancer-data/exp_seq.', files, sep='')
cancerfiles = paste(cancerfiles, '-US.tsv', sep='')
cancerbmi = paste('./clindata/TCGA/', tolower(files), sep='')
cancerbmi = paste(cancerbmi, '_clinical_patient.txt', sep='')
for (i in 1:length(files)) {
	## process sequence data:
	seq = read.table(cancerfiles[i], sep = '\t', header=T)
	seq = tbl_df(seq)
	dup = duplicated(paste(seq$submitted_sample_id,seq$gene_id,sep=''))
	seq = seq[!dup,c('submitted_sample_id','gene_id','raw_read_count')]
	genes = unique(seq$gene_id)
	seq  = spread(seq,submitted_sample_id,raw_read_count)
	rownames(seq) = genes
	seq  = data.matrix(seq)
	seq = seq[-1,-1]
	seq = icgc_to_tcga(t(seq))
    if (length(rownames(seq)) > length(unique(rownames(seq)))) {
        seq = collapseRows(seq, unique(rownames(seq)), unique(rownames(seq)))
        seq = seq$datETcollapsed
    }

	## process clinical data:
	clin = read.csv(cancerbmi[i], sep = '\t', skip=1, header=T)
	clin = clin[-1,]
	rownames(clin) = clin$bcr_patient_barcode

	## Remove samples without height and/or weight data
	if (length(which(clin$height == "[Not Available]")) > 0){
		clin = clin[-which(clin$height == "[Not Available]"),]
	}
	if (length(which(clin$weight == "[Not Available]")) > 0){
		clin = clin[-which(clin$weight == "[Not Available]"),]
	}

	## Remove samples that are not in the sequence data:
	seq = seq[which(rownames(seq) %in% rownames(clin)),]
	clin = clin[which(rownames(clin) %in% rownames(seq)),]

	## Calculate BMI and BMI status:
	bmi = as.numeric(as.vector(clin$weight))/((as.numeric(as.vector(clin$height))/100)^2)
	bmiStatus = bmi
	bmiStatus[bmiStatus >= 30] = "obese"
	bmiStatus[bmiStatus < 25] = "normal"
	bmiStatus[(bmiStatus < 30) & (bmiStatus >= 25)] = "overweight"

	## Add BMI information to the clinical data:
	clin = cbind(clin, bmi, bmiStatus)

	## assign variable names:
	assign(paste(files[i], 'raw', sep = ''), seq)
	assign(paste(files[i], 'clin', sep=''), clin)
}

for (i in 1:length(files)) {
	txt = paste(files[i], 'raw', sep='')
	mat = t(get(txt))
	mat = standardise_data(mat)
	txt = paste(files[i], 'std', sep='')
	assign(txt, t(mat))
}

## TODO: clean the sample names of ICGC/TCGA and Gatza data
###############################################################################
# Make a matrix of gene symbols for checking the direction of specific genes

crsymrma    = make_sym_mat(crrma)
crsymmas    = make_sym_mat(crmas)
fmsymrma    = make_sym_mat(fmrma)
fmsymmas    = make_sym_mat(fmmas)
crissymrma  = make_sym_mat(crisrma)
crissymmas  = make_sym_mat(crismas)
gatzasymrma = make_sym_mat(gatzarma)
gatzasymmas = make_sym_mat(gatzamas)

crsymrmastd = standardise_data(crsymrma, log = F)
crsymmasstd = standardise_data(crsymmas, log = F)
fmsymrmastd = standardise_data(fmsymrma, log = F)
fmsymmasstd = standardise_data(fmsymmas, log = F)
crissymrmastd = standardise_data(crissymrma, log = F)
crissymmasstd = standardise_data(crissymmas, log = F)
gatzasymrmastd = standardise_data(gatzasymrma, log = F)
gatzasymmasstd = standardise_data(gatzasymmas, log = F)

###############################################################################
## Pathway data base:

#Import Human Gene Symbols
SYMBOL.list = as.list(org.Hs.egSYMBOL)

KEGG.list = as.list(org.Hs.egPATH) ##Import KEGG pathways:
names(KEGG.list) = unlist(SYMBOL.list) #name KEGG pathway lists with corresponding gene symbols
keggpath = as.list(KEGGPATHID2NAME) #mapping KEGG path IDs to human read pathway name

GO.list = as.list(org.Hs.egGO) ##Import GO pathways:
tmp = GO.list #... and reformat so matching KEGG.list

# pull out the GOID
for (i in 1:length(tmp)) {
    ind = names(tmp)[i]
    if(class(tmp[[ind]]) == 'list') {
        tmp[[ind]] = names(tmp[[ind]])
    }
}
tmp = tmp[which(!is.na(tmp))]##filter out the NA values
GO.list = tmp
names(GO.list) = SYMBOL.list[names(GO.list)] #name GO pathway lists with corresponding gene symbols
gopath = as.list(GOTERM) #mapping GO IDs to human read pathway name
#... and reformat so matching KEGG.list
tmp = list()
for(i in 1:length(gopath)) tmp[[i]] = gopath[[i]]@Term
names(tmp) = names(gopath)
gopath = tmp

## TODO: Do I really need this bit of code?
#Import Human Reactome pathways
reactome.list = as.list(reactomeEXTID2PATHID)
reactome.list = reactome.list[order(as.numeric(names(reactome.list)))] #Sort the names in the list:
reactome.list = reactome.list[which(names(reactome.list) %in% names(SYMBOL.list))] #get paths that have gene symbols:
names(reactome.list) = unlist(SYMBOL.list[names(reactome.list)]) #rename the entrez gene ID into gene symbol
reactomepath = as.list(reactomePATHID2NAME) #mapping reactome path IDs to human read pathway name
reactomepath = reactomepath[grep('Homo sapiens',reactomepath)] #pull out all human-related pathways
reactomepath = lapply(reactomepath, function(x) gsub('Homo sapiens: ', '', x)) #cut out the 'Homo sapiens: ' bit so it's only the pathway names.
reactomepath = lapply(reactomepath, function(x) return(x[1])) #cut out the 'Homo sapiens: ' bit so it's only the pathway names.

save.image()

###############################################################################
## Creighton's obesity-associated metagene analysis

# source("creighton_metagene.R")

## Check if Creighton's obesity-associated metagenes are actually associated to
## obesity:

## Get the genes that are common from the creighton gene set and the ICGC data set, then use these genes to make the matrix for SVD/metagene analysis
commongenes = mapIds(hgu133a.db, keys = crobsgenes, column = 'SYMBOL', keytype = "PROBEID", multiVals = 'first')
commongenes = commongenes[which(!is.na(commongenes))]
commongenes = unique(commongenes)
commongenes = commongenes[which(commongenes %in% colnames(BLCAraw))]

## Make matrix with just the obesity-associated genes:
## Use RMA normalised data, and use standardised or non-standardised data:
crobsmatraw = crsymrma[commongenes,]
crobsmatstd = crsymrmastd[commongenes,]

## 644 gene symbols common in Creighton's obesity-associated genes and the ICGC data:
dim(crobsmatraw)
dim(crobsmatstd)

## Apply SVD and generate metagene from the Creighton data:
## The metagene from standardised data was flipped, so the metagene was in the same direction as the BMI
crmetaraw = make_metagene(crobsmatraw, flip = F, raw = T)
crmetastd = make_metagene(crobsmatstd, flip = T, raw = T)

# Check if the raw or adjusted metagenes have any differences:
pdf('pdf/creighton_mg_results/check_raw_vs_std_mg.pdf')
plot(crmetastd, crmetaraw, pch=20, main='Creighton metagene comparison (raw values)', ylab='Metagene values from raw data', xlab='Metagene values from standardised data')

# Look for differences between ranked metagene:
crmetaraw = make_metagene(crobsmatraw, flip = F, raw = F)
crmetastd = make_metagene(crobsmatstd, flip = T, raw = F)
plot(crmetastd, crmetaraw, pch=20, main='Creighton metagene comparison (ranked values)', ylab='Metagene values from raw data', xlab='Metagene values from standardised data')
dev.off()

## Check if the raw metagene correlates with BMI:
pdf('pdf/creighton_mg_results/creighton_mg_heatmap1.pdf')
main='Creighton metagene (raw)'
# For the heatmap to show up nicely, standardised data needs to be used:
mgheatmap(crobsmatstd, crmetaraw, main = main)
bmiplot(crclin, crmetaraw, main = main)

# repeat with standardised data:
main='Creighton metagene (standardised)'
mgheatmap(crobsmatstd, crmetastd, main = main)
bmiplot(crclin, crmetastd, main = main)
dev.off()

## Keep the ordering of the genes from Creighton's original obesity-associated
## gene heatmap:
dend = get_dend(crobsmatstd)

## Make transformation matrix in both data:
crtransmatraw = make_trans_mat(crobsmatraw)
crtransmatstd = make_trans_mat(crobsmatstd)

###############################################################################
## Check if Creighton's obesity-associated metagene transfers to ICGC cancer
## types

## Pull out the obesity-associated genes from the ICGC data:
cancertypes = c('BLCA', 'CESC', 'COAD', 'KIRP', 'LIHC', 'READ', 'SKCM', 'UCEC')
for (i in 1:length(cancertypes)) {
	var = paste(cancertypes[i], 'raw', sep='')
	seq = get(var)
	seq = seq[,commongenes]
	seq = t(seq) ## make it genes by samples
	txt = paste(cancertypes[i], 'cr', sep='')
	assign(txt, seq)
}

## check if it has the right number of genes:
dim(BLCAcr) # 644 genes

## log10 and standardise data:
for (i in 1:length(cancertypes)) {
	txt = paste(cancertypes[i], 'cr', sep='')
	seq = get(txt)
	seq = standardise_data(seq)
	txt = paste(txt, 'std', sep='')
	assign(txt, seq)
}

## Check for any differences in the metagenes produced from raw or standardised
## transformation matrix in raw or standardised data:

pdf('pdf/creighton_mg_results/crtcga_raw_vs_std.pdf')
for (i in 1:length(cancertypes)) {
	txt = paste(cancertypes[i], 'cr', sep='')
	mat = get(txt)
	txt = paste(cancertypes[i], 'crstd', sep='')
	matstd = get(txt)
	meta1 = make_transmeta(mat, crtransmatraw, flip = F)
	meta2 = make_transmeta(mat, crtransmatstd, flip = F)
	meta3 = make_transmeta(matstd, crtransmatraw, flip = F)
	meta4 = make_transmeta(matstd, crtransmatstd, flip = F)
	plot(meta1, meta2, pch=20, main="Raw/raw vs. raw/std FM metagene in Creighton's data (ranked values)", xlab='Metagene values from raw/raw TM', ylab='Metagene values from raw/std TM')
	plot(meta1, meta3, pch=20, main="Raw/raw vs. std/raw FM metagene in Creighton's data (ranked values)", xlab='Metagene values from raw/raw TM', ylab='Metagene values from std/raw TM')
	plot(meta2, meta4, pch=20, main="Raw/std vs. std/std FM metagene in Creighton's data (ranked values)", xlab='Metagene values from raw/std TM', ylab='Metagene values from std/std TM')
	plot(meta3, meta4, pch=20, main="Std/raw vs. std/std FM metagene in Creighton's data (ranked values)", xlab='Metagene values from std/raw TM', ylab='Metagene values from std/std TM')
}
dev.off()

## The above results show how the standardised data is the best one to use,
## regardless of the transformation matrix used.

## Transform the ICGC data:
for (i in 1:length(cancertypes)) {
	txt = paste(cancertypes[i], 'crstd', sep='')
	mat = get(txt)
	t = make_transmeta(mat, crtransmatstd, flip = F)
	txt = paste(cancertypes[i], 'crmetastd', sep='')
	assign(txt, t)
}

## Check if the metagene correlates with sample gene expression and/or BMI:
pdf('pdf/creighton_mg_results/crtcga_std.pdf')
for (i in 1:length(cancertypes)) {
	dat  = get(paste(cancertypes[i], 'crstd', sep=''))
	meta = get(paste(cancertypes[i], 'crmetastd', sep=''))
	clin = get(paste(cancertypes[i], 'clin', sep=''))
	main = paste('Creighton metagene (', cancertypes[i], sep = "")
	main = paste(main, ')', sep = "")
	mgheatmap(dat, meta, dend = dend, main = main)
	bmiplot(clin, meta, main = main)
}
dev.off()

###############################################################################
## Check if Creighton's obesity-associated metagene transfers to Cris' breast
## cancer data

## Extract Creighton's obesity-associated genes from Cris' data (RMA normalised
## and standardised):
crisobsmat = crissymrmastd[commongenes,]

## Transform the data:
cristransmetacr = make_transmeta(crisobsmat, crtransmatstd, flip = F)

## Check if the metagene correlates with sample gene expression and/or BMI:
pdf('pdf/creighton_mg_results/cris_cr_trans_meta.pdf')
main = "Creighton metagene in Print's data"
mgheatmap(crisobsmat, cristransmetacr, dend = dend, main = main)
bmiplot(crisclin, cristransmetacr, main = main)
dev.off()

###############################################################################
## Creighton gene expression analysis

# source("creighton_DEG.R")

## Make groups (obesity and non-obesity groups) for DEG analysis:
group = crclin$bmiStatus
group = ifelse(group=='obese', 'Obese', 'Non-obese')

## Find the DEGs in the RMA normalised data (not standardised):
crdeg = get_deg(crrma, group, p = 0.05)
sum(crdeg$P.Value < 0.05)   # 5278 genes
sum(crdeg$adj.P.Val < 0.05) # 9 genes
sum(crdeg$P.Value < 0.01)   # 1781 genes
sum(crdeg$adj.P.Val < 0.01) # 0 genes
length(which(2^crdeg$logFC > 1.2)) # only 61 genes had > 1.2 fold change with p < 0.05

## Look for any overlap of genes between the DEGs I have found and Creighton's
## original obesity-associated genes:

## Take the top 799 DEGs:
rawobsgene = rownames(crdeg[1:799,])
crolgene = rawobsgene[which(rawobsgene %in% crobsgenes)]
length(crolgene) ## 239 genes

## Do the same, but in residual data:

## Fit linear model on the data using all the other clinical variables:
residuals = crrma
residuals = t(apply(residuals, 1, function(x) lm(x ~ crclin$age + crclin$race + crclin$menopause + crclin$grade + crclin$LNstatus + crclin$ERstatus + crclin$PRstatus + crclin$HER2status)$residuals))

## Fit linear model on the residuals, using obese vs non-obese design matrix
crresdeg = get_deg(residuals, group, p = 0.05)
sum(crresdeg$P.Value < 0.05)   # 4371 genes
sum(crresdeg$adj.P.Val < 0.05) # 0 genes
sum(crresdeg$P.Value < 0.01)   # 1104 genes
sum(crresdeg$adj.P.Val < 0.01) # 0 genes

## Take the top 799 DEGs:
resobsgene = rownames(crresdeg[1:799,])
rescrolgene = resobsgene[which(resobsgene %in% crobsgenes)]
length(rescrolgene) ## 168 genes

## Do the same, but in Caucasian-only data:
camat  = crrma
ind    = which(crclin$race == 'w')
camat  = camat[,ind] ## 77 caucasian samples in total
caclin = crclin[ind,]

group = caclin$bmiStatus
group = ifelse(group == 'obese', 'Obese', 'Non-obese')
cadeg = get_deg(camat, group, p = 0.05)
sum(cadeg$P.Value < 0.05)   # 6029 genes
sum(cadeg$adj.P.Val < 0.05) # 0 genes
sum(cadeg$P.Value < 0.01)   # 2129 gene
sum(cadeg$adj.P.Val < 0.01) # 0 gene

## Take the top 799 DEGs:
caobsgene = rownames(cadeg[1:799,])
cacrolgene = caobsgene[which(caobsgene %in% crobsgenes)]
length(cacrolgene) ## 148 genes

## Fit linear model on the Caucasian-only data after removing all the other clinical variables:
caresiduals = camat
caresiduals = t(apply(caresiduals, 1, function(x) lm(x ~ caclin$age + caclin$menopause + caclin$grade + caclin$LNstatus + caclin$ERstatus + caclin$PRstatus + caclin$HER2status)$residuals))

## Fit linear model on the residuals, using obese vs non-obese design matrix
caresdeg = get_deg(caresiduals, group)
sum(caresdeg$P.Value < 0.05)   # 5427 genes
sum(caresdeg$adj.P.Val < 0.05) # 0 genes
sum(caresdeg$P.Value < 0.01)   # 1558 gene
sum(caresdeg$adj.P.Val < 0.01) # 0 gene

## Take the top 799 DEGs:
caresobsgene = rownames(caresdeg[1:799,])
carescrolgene = caresobsgene[which(caresobsgene %in% crobsgenes)]
length(carescrolgene) ## 92 genes

pdf('pdf/creighton_mg_results/deg_venn.pdf')
## Make venn diagram with res, raw, and cr obesity gene:
setlist = list("DEGs from\nresidual data" = resobsgene, "DEGs from\nCreighton's data" = rawobsgene, "Creighton's original\nobesity genes" = crobsgenes)
OLlist = overLapper(setlist=setlist, sep='_', type='vennsets')
counts = sapply(OLlist$Venn_List, length)
vennPlot(counts = counts)

## Make a Venn diagram with ca, cares, and cr obesity gene:
setlist = list("DEGs from Caucasian-only\nresidual data" = caresobsgene , "DEGs from\nCaucasian-only data" = caobsgene, "Creighton's original\nobesity genes" = crobsgenes)
OLlist = overLapper(setlist=setlist, sep='_', type='vennsets')
counts = sapply(OLlist$Venn_List, length)
vennPlot(counts = counts)
dev.off()

###############################################################################
## Check if the obesity-associated genetic signatures correlate with sample
## BMI and BMI status in Creighton's data.

## For each metagene, identify the gene that are in the ICGC data and use these
## gene for validataion
dim(crsymrmastd) #13031 gene
icgcgene = colnames(BLCAraw)

allobsname = c("rawobsgene","crolgene","resobsgene","rescrolgene", "caobsgene","cacrolgene","caresobsgene","carescrolgene")
metalength = matrix(1,8)
rownames(metalength) = allobsname
for (i in 1:length(allobsname)) {
	tmp = get(allobsname[i])
	tmp = mapIds(hgu133a.db, keys = tmp, column = 'SYMBOL', keytype = "PROBEID", multiVals = 'first')
	tmp = unique(tmp)
	tmp = tmp[which(!is.na(tmp))]
	tmp = tmp[which(tmp %in% rownames(crsymrmastd))]
	tmp = tmp[which(tmp %in% icgcgene)]
	assign(allobsname[i], tmp)
	metalength[i,1] = length(tmp)
}

metalength
# rawobsgene     678
# crolgene       199
# resobsgene     655
# rescrolgene    147
# caobsgene      657
# cacrolgene     128
# caresobsgene   651
# carescrolgene  86

# Need to first check if all the metagene are going in the same direction

# Find the genes that are common across all of the metagenes:
namelist = c("Raw metagene","CR overlap metagene (raw)","Residual metagene","CR overlap metagene (residual)","Raw metagene (Caucasian/raw)","CR overlap metagene (Caucasian/raw)","Residual metagene (Caucasian/residual)","CR overlap metagene (Caucasian/residual)")
checkgene= c(rawobsgene, crolgene, resobsgene, rescrolgene, caobsgene, cacrolgene, caresobsgene, carescrolgene)
checkgene = table(checkgene)[which(table(checkgene) == length(allobsname))]
checkgene = names(checkgene)

# Check for the direction of each of the metagene, using the common genes
# identified above
flip = c('resobsgene', 'rescrolgene', 'caresobsgene')
pdf('pdf/creighton_mg_results/cr_meta_direction.pdf')
for (i in 1:length(allobsname)) {
	gene = get(allobsname[i])
	mat = crsymrmastd[gene,]
	tmpmeta = make_metagene(mat, flip = F, raw = F)
	if (allobsname[i] %in% flip) {
		tmpmeta = 1 - tmpmeta
	}
	mat = mat[checkgene,]
	mgheatmap(mat, tmpmeta, main = namelist[i])
}
dev.off()

# Residual, residual overlap, and Caucasian-only residual metagenes must be
# flipped so all the metagenes are going in the same direction

## Validate these metagenes in creighton data:
metalist  = list()
pdf('pdf/creighton_mg_results/cr_deg_meta_vs_clin.pdf')
for (i in 1:length(allobsname)) {
	gene    = get(allobsname[i])
	mat     = crsymrmastd[gene,]
	tmpmeta = make_metagene(mat, flip = F, raw = F)
	if (allobsname[i] %in% flip) {
		tmpmeta = 1-tmpmeta
	}
	metalist[[i]] = tmpmeta
	mgheatmap(mat, tmpmeta, main = namelist[i])
	bmiplot(crclin, tmpmeta, main = namelist[i])

	txt = gsub("gene", "transmat", allobsname[i])
	trmat = make_trans_mat(mat)
	assign(txt, trmat)
}
dev.off()

## Get row dendrogram for each obesity-associated gene signature, and use these
## to order the ICGC data:
dendlist = list()
for (i in 1:length(allobsname)) {
	mat = crsymrmastd[get(allobsname[i]),]
	dendlist[[i]] = get_dend(mat)
}

## See what the correlation between the metagenes are like:
pdf('pdf/creighton_mg_results/cr_meta_cor.pdf')
x = matrix(1,103)
for(i in 1:length(metalist)) {
	x = cbind(x, metalist[[i]])
}
x = x[,-1]
colnames(x) = gsub('gene', '', allobsname)
cormat = cor(x, method = 'pearson')
cormatadj = x[,-which(colnames(x) == "caresobs")]
cormatadj = cor(cormatadj, method = 'pearson')
heatmap.2(cormat, trace = 'none', scale='none', col='bluered', cexRow = 1.0, cexCol = 1.0)
heatmap.2(cormatadj, trace = 'none', scale='none', col='bluered', cexRow = 1.0, cexCol = 1.0)
dev.off()

###############################################################################
## Check if these obesity-associated genetic signatures correlate with sample
## BMI and BMI status in ICGC data.

alltransname = gsub('gene', 'transmat', allobsname)
cancertypes = c('BLCA', 'CESC', 'COAD', 'KIRP', 'LIHC', 'READ', 'SKCM', 'UCEC')

for (i in 1:length(allobsname)) {
	pdfname = paste('pdf/creighton_mg_results/', allobsname[i], sep='')
	pdfname = paste(pdfname, 'ICGC.pdf', sep='_')
	pdf(pdfname)
	gene = get(allobsname[i])
	transmat = get(alltransname[i])
	for (j in 1:length(cancertypes)) {
		txt = paste(cancertypes[j], 'std', sep='')
		cancer = t(get(txt))
		mat = cancer[gene,]
		clin = get(paste(cancertypes[j], 'clin', sep=''))
		cancermeta = make_transmeta(mat, transmat, flip = F)
		if (allobsname[i] %in% flip) {
			cancermeta = 1-cancermeta
		}
		main = paste(namelist[i], '(')
		main = paste(main, cancertypes[j], sep = '')
		main = paste(main, ')', sep = '')
		dend = dendlist[[i]]
		mgheatmap(mat, cancermeta, dend = dend, main = main)
		bmiplot(clin, cancermeta, main = main)
	}
	dev.off()
}

###############################################################################
## Check if Creighton's obesity-associated metagene transfers to Cris' breast
## cancer data

## Check if the metagene correlates with sample gene expression and/or BMI:
pdf('pdf/creighton_mg_results/cris_crdeg_trans_meta.pdf')
for (i in 1:length(allobsname)) {
	mat = crissymrmastd[get(allobsname[i]),]
	transmat = get(alltransname[i])
	meta = make_transmeta(mat, transmat, flip = F)
	if (allobsname[i] %in% flip) {
		meta = 1 - meta
	}
	main = paste(namelist[i], "(Cris' data)")
	dend = dendlist[[i]]
	mgheatmap(mat, meta, dend = dend, main = main)
	bmiplot(crisclin, meta, main = main)
}
dev.off()

###############################################################################
## Continuous BMI genes from Creighton's data

# source('cont_bmi.R')

## Get obesity-associated genes from the correlation of the genes with the BMI
## values, rather than the sample BMI status

# Use RMA normalised Creighton data:
mat = crrma
dim(mat) # 22283 103

# Correlate it with the BMI value of the samples:
bmicor = cor(t(mat), crclin$bmi, method = "spearman")
colnames(bmicor) = "correlation"

pval = apply(mat, 1, function(x) cor.test(crclin$bmi, x)$p.value)
pval = p.adjust(pval, method = 'fdr')

bmicor = cbind(bmicor, p.value = pval)

# Are there any genes significantly correlating with sample BMI?
which(bmicor[,2] <= 0.05) %>% length ## 0 genes significant

###############################################################################
## Compare the top correlating genes from the list and compare it with other
## metagenes

## Take the top 799 BMI correlated genes from the list:
abscor = abs(bmicor[,1])
abscor = rev(sort(abscor))
genes = names(abscor)[1:799]

## See if the genes are similar to Creighton's or FM's gene signatures:
length(which(genes %in% crobsgenes)) ## 218 genes common
length(which(genes %in% as.vector(fmobsgenes$Probe))) ## 4 genes in common

## Make gene symbol version of the genes:
contgenes = mapIds(hgu133a.db, keys = genes, column = 'SYMBOL', keytype = "PROBEID", multiVals = 'first')
contgenes = unique(contgenes)
contgenes = contgenes[which(!is.na(contgenes))]
contgenes = contgenes[which(contgenes %in% rownames(crsymrmastd))]
contgenes = contgenes[which(contgenes %in% colnames(BLCAraw))]
length(contgenes) ## 665 genes

## See if the genes correlate with BMI status and/or BMI in Creighton's data:
pdf('pdf/creighton_mg_results/contbmi_cr.pdf')
crcontmat = crsymrmastd[contgenes,]
meta = make_metagene(crcontmat, flip = F, raw = F)
main = "Continuous BMI genes (Creighton's data)"
mgheatmap(crcontmat, meta, main = main)
bmiplot(crclin, meta, main = main)
dev.off()

dend = get_dend(crcontmat)

contbmitransmat = make_trans_mat(crcontmat)

## See if the genes correlate with BMI status and/or BMI in Cris' data:
pdf('pdf/creighton_mg_results/contbmi_cris.pdf')
criscontmat = crissymrmastd[contgenes,]
meta = make_transmeta(criscontmat, contbmitransmat)
main = "Continuous BMI genes (Cris' data)"
mgheatmap(criscontmat, meta, main = main)
bmiplot(crisclin, meta, main = main)
dev.off()

## See if the genes correlate with BMI status and/or BMI in ICGC data:
pdf('pdf/creighton_mg_results/contbmi_ICGC.pdf')
for (i in 1:length(cancertypes)) {
	mat = t(get(paste(cancertypes[i], 'std', sep='')))
	mat = mat[contgenes,]
	clin = get(paste(cancertypes[i], 'clin', sep=''))
	meta = make_transmeta(mat, contbmitransmat)
	main = paste("Continuous BMI genes (", cancertypes[i], sep='')
	main = paste(main, ')', sep = '')
	mgheatmap(mat, meta, main = main)
	bmiplot(clin, meta, main = main)
}
dev.off()

###############################################################################
## Try top 100 genes, instead of top 799 genes

abscor = abs(bmicor[,1])
abscor = rev(sort(abscor))
genes = names(abscor)[1:100]

## See if the genes are similar to Creighton's or FM's gene signatures:
length(which(genes %in% crobsgenes)) ## 38 genes common
length(which(genes %in% as.vector(fmobsgenes$Probe))) ## 1 gene in common

## Make gene symbol version of the genes:
contgenes = mapIds(hgu133a.db, keys = genes, column = 'SYMBOL', keytype = "PROBEID", multiVals = 'first')
contgenes = unique(contgenes)
contgenes = contgenes[which(!is.na(contgenes))]
contgenes = contgenes[which(contgenes %in% rownames(crsymrmastd))]
contgenes = contgenes[which(contgenes %in% colnames(BLCAraw))]
length(contgenes) ## 83 genes

## See if the genes correlate with BMI status and/or BMI in Creighton's data:
pdf('pdf/creighton_mg_results/contbmi100_cr.pdf')
crcontmat = crsymrmastd[contgenes,]
meta = make_metagene(crcontmat, flip = F, raw = F)
main = "Continuous BMI genes (Creighton's data)"
mgheatmap(crcontmat, meta, main = main)
bmiplot(crclin, meta, main = main)
dev.off()

dend = get_dend(crcontmat)

contbmitransmat = make_trans_mat(crcontmat)

## See if the genes correlate with BMI status and/or BMI in Cris' data:
pdf('pdf/creighton_mg_results/contbmi100_cris.pdf')
criscontmat = crissymrmastd[contgenes,]
meta = make_transmeta(criscontmat, contbmitransmat)
main = "Continuous BMI genes (Cris' data)"
mgheatmap(criscontmat, meta, main = main)
bmiplot(crisclin, meta, main = main)
dev.off()

## See if the genes correlate with BMI status and/or BMI in ICGC data:
pdf('pdf/creighton_mg_results/contbmi100_ICGC.pdf')
for (i in 1:length(cancertypes)) {
	mat = t(get(paste(cancertypes[i], 'std', sep='')))
	mat = mat[contgenes,]
	clin = get(paste(cancertypes[i], 'clin', sep=''))
	meta = make_transmeta(mat, contbmitransmat)
	main = paste("Continuous BMI genes (", cancertypes[i], sep='')
	main = paste(main, ')', sep = '')
	mgheatmap(mat, meta, main = main)
	bmiplot(clin, meta, main = main)
}
dev.off()

###############################################################################
## Fuentes-Mattei's obesity-associated metagene anlysis

# source("fm_metagene.R")

## Get the genes that are common from the creighton data set and the ICGC data
## set and use these genes to make the FM metagene:
tmpgenes = mapIds(hgu133a.db, keys = as.vector(fmobsgenes$Probe), column = 'SYMBOL', keytype = "PROBEID", multiVals = 'first')
tmpgenes = tmpgenes[which(!is.na(tmpgenes))]
tmpgenes = unique(tmpgenes)
tmpgenes = tmpgenes[which(tmpgenes %in% colnames(BLCAraw))]
tmpgenes = tmpgenes[which(tmpgenes %in% rownames(fmsymrma))]

## Make gene symbol matrix with FM's obesity-associated genes:
fmobsmat    = fmsymrma[tmpgenes,]
fmobsmatstd = fmsymrmastd[tmpgenes,]
dim(fmobsmat) ## 106 obesity-associated genes
dim(fmobsmatstd) ## 106 genes by 278 samples

## Check the metagene in FM's data first:
fmmeta = make_metagene(fmobsmat, flip = F, raw = T)
fmmetastd = make_metagene(fmobsmatstd, flip = F, raw = T)

# Check if the raw or standardised metagenes have any difference:
pdf('pdf/fm_mg_results/fm_meta_check.pdf')
plot(fmmetastd, fmmeta, pch=20, main='Fuentes-Mattei metagene\ncomparison (raw values)', ylab='Metagene values from raw data', xlab='Metagene values from standardised data')

# Look at the ranked metagenes:
fmmeta = make_metagene(fmobsmat, flip = F, raw = F)
fmmetastd = make_metagene(fmobsmatstd, flip = F, raw = F)
plot(fmmetastd, fmmeta, pch=20, main='Fuentes-Mattei metagene\ncomparison (ranked values)', ylab='Metagene values from raw data', xlab='Metagene values from standardised data')

# For the heatmap to show up nicely, standardised data needs to be used:
main='FM metagene (raw)'
mgheatmap(fmobsmatstd, fmmeta, main = main)
# Repeat with standardised metagene:
main='FM metagene (standardised)'
mgheatmap(fmobsmatstd, fmmetastd, main = main)
dev.off()

## The metagenes seem to be going all over the place...
## Maybe the FM obesity-associated genes aren't good

## Get the row dendrogram of the original heatmap:
dend = get_dend(fmobsmatstd)

## Make transformation matrix:
fmstdtransmat = make_trans_mat(fmobsmatstd)

###############################################################################
## Get FM metagene from Creighton data using the transformation matrix
##
## NOTE: Transformation matrix from FM data was used to obtain the metagene in
## Creighton data, rather than using SVD with FM's obesity-associated genes.
## This was to ensure the obesity-associated metagene from FM's study was
## transferrable to other data sets, and not dependent on the data's expression
## values (I think).

## Make a matrix in Creighton's data (standardised):
crfmobsmat = crsymrmastd[tmpgenes,]
dim(crfmobsmat) ## 106 genes

crfmtransmeta = make_transmeta(crfmobsmat, fmstdtransmat, flip = F)

pdf('pdf/fm_mg_results/cr_fm_meta.pdf')

main="FM metagene from Creighton's data"
mgheatmap(crfmobsmat, crfmtransmeta, dend = dend, main = main)
bmiplot(crclin, crfmtransmeta, main = main)

dev.off()

###############################################################################
## Try it on ICGC data

## Pull out all the genes common in Creighton obesity genes and ICGC data:
cancertypes = c('BLCA', 'CESC', 'COAD', 'KIRP', 'LIHC', 'READ', 'SKCM', 'UCEC')
for (i in 1:length(cancertypes)) {
	txt = paste(cancertypes[i], 'std', sep='')
	seq = t(get(txt))
	seq = seq[tmpgenes,]
	txt = paste(cancertypes[i], 'fm', sep='')
	assign(txt, seq)
}

## Transform the ICGC data:
for (i in 1:length(cancertypes)) {
	txt = paste(cancertypes[i], 'fm', sep='')
	mat = get(txt)
	t   = make_transmeta(mat, fmstdtransmat, flip = T)
	txt = paste(cancertypes[i], 'fmmeta', sep='')
	assign(txt, t)
}

## Check if the metagene correlates with sample gene expression and/or BMI in
## ICGC data:
pdf('pdf/fm_mg_results/fm_meta_ICGC.pdf')
for (i in 1:length(cancertypes)) {
	dat  = get(paste(cancertypes[i], 'fm', sep=''))
	meta = get(paste(cancertypes[i], 'fmmeta', sep=''))
	clin = get(paste(cancertypes[i], 'clin', sep=''))
	main = paste('FM metagene (', cancertypes[i], sep = "")
	main = paste(main, ')', sep = "")
	mgheatmap(dat, meta, dend = dend, main = main)
	bmiplot(clin, meta, main = main)
}
dev.off()

###############################################################################
## Try it on Cris's data

## Make a matrix in Creighton's data (standardised):
crisfmobsmat = crissymrmastd[tmpgenes,]
dim(crisfmobsmat) ## 106 genes

crisfmtransmeta = make_transmeta(crisfmobsmat, fmstdtransmat, flip = F)

pdf('pdf/fm_mg_results/cris_fm_meta.pdf')

main="FM metagene from Print's data"
mgheatmap(crisfmobsmat, crisfmtransmeta, dend = dend, main = main)
bmiplot(crisclin, crisfmtransmeta, main = main)

dev.off()

###############################################################################
## Gene expression analysis on ICGC data to see if there are any genes common
## across multiple cancer types

# source("icgc_deg.R")

## Get DEGs for each ICGC cancer type:
## Use obese vs. non-obese as the group design
icgc_deg_list = list()
for (i in 1:length(cancertypes)) {
	tmpmat     = get(paste(cancertypes[i], 'raw', sep=''))
	tmpclin    = get(paste(cancertypes[i], 'clin', sep=''))
	group      = ifelse(tmpclin$bmiStatus == 'obese', 'obese', 'non-obese')
	tmpmatnorm = voom_norm(t(tmpmat), group)
	tmpdeg     = get_deg(tmpmatnorm, group)
	icgc_deg_list[[i]] = tmpdeg
}

## How many DEGs were found in each cancer type?
for (i in 1:length(cancertypes)) {
	print(nrow(icgc_deg_list[[i]]))
}
# BLCA 679  genes
# CESC 1229 genes
# COAD 974  genes
# KIRP 687  genes
# LIHC 3340 genes
# READ 796  genes
# SKCM 1137 genes
# UCEC 2934 genes

## How many genes are common to multiple cancer types?
degs = vector()
for (i in 1:length(cancertypes)) {
	degs = c(degs, rownames(icgc_deg_list[[i]]))
}
deg_table = table(degs)
table(deg_table)
#    1    2    3    4    5  total genes
# 7093 1886  266   27    1		9273

## There were 9273 unique DEGs found from all 8 cancer types.
## Of those, 7093 genes were unique to one cancer type, 1886 genes were common
## to two cancer types, 266 genes were common to three cancer types, 27 genes
## common to four cancer types, and 1 gene was common to five cancer types.

###############################################################################
## Validate the statistical significance by repeating it 1000 times with
## randomly chosen BMI status.

set.seed(1) # set seed for reproducibility

## Simulate 1000 times:
cancerdegres = matrix(0, 8, 0) ## empty matrix for the results
for (i in 1:1000) {
	tmplist = list()
	for (j in 1:length(cancertypes)) {
		tmpmat     = get(paste(cancertypes[j], 'raw', sep=''))
		tmpclin    = get(paste(cancertypes[j], 'clin', sep=''))
		group      = sample(tmpclin$bmiStatus, length(tmpclin$bmiStatus))
		group      = ifelse(group == 'obese', 'obese', 'non-obese')
		tmpmatnorm = voom_norm(t(tmpmat), group)
		tmpdeg     = get_deg(tmpmatnorm, group)
		tmplist[[j]] = rownames(tmpdeg)
	}
	## get table of common genes
	x = unlist(tmplist)
	x = table(table(x))

	## make an empty vector of length 8 and store the result:
	v = c(rep(0, 8))
	for (k in 1:length(x)) {
		v[k] = x[k]
	}

	if (i %% 50 == 0) {
		print(paste("Done", i))
	}

	## store data in matrix:
	cancerdegres = cbind(cancerdegres, v)
}
m = apply(cancerdegres, 1, mean)
percentile = apply(cancerdegres, 1, function(x) quantile(x, 0.95))

## TODO: re-run the simulation...

## TODO: tabulate the simulated results
# No. cancer      & 1        & 2       & 3      & 4     & 5     & 6     & 7     & 8
# Mean            & 5263.821 & 848.395 & 77.461 & 4.546 & 0.166 & 0.003 & 0.000 & 0.000
# 95th percentile &  & & & & & & &

###############################################################################
## Pathway enrichment analyses on all of the DEGs found so far (Creighton, FM,
## ICGC data)

# Make a gene-by-pathway matrix to use it in the pathway enrichment analysis

# KEGG matrix
pathnames = unique(unlist(keggpath)) #Get all the pathways in the keggpath list
genenames = colnames(BLCAraw) #Get all the gene names in the cancer data
tmp = matrix(nrow = length(genenames), ncol = length(pathnames)) #make an empty matrix to fill
colnames(tmp) = pathnames
rownames(tmp) = genenames

for(i in 1:length(rownames(tmp))) {
    v = colnames(tmp) %in% lookup(rownames(tmp)[i], "KEGG")
    if(length(v) > 0) {
        tmp[rownames(tmp)[i],] = v
    }
}

keggTFmat = ifelse(tmp == F, 0, 1)

# Reactome matrix
pathnames = unique(unlist(reactomepath))
tmp = matrix(nrow = length(genenames), ncol = length(pathnames))
colnames(tmp) = pathnames
rownames(tmp) = genenames

for(i in 1:length(rownames(tmp))) {
    v = colnames(tmp) %in% lookup(rownames(tmp)[i], "reactome")
    if(length(v) > 0) {
        tmp[rownames(tmp)[i],] = v
    }
}

reactomeTFmat = ifelse(tmp == F, 0, 1)

# GO matrix
# For GO, there are too many pathways in total, so will have to make a matrix
# using the pathways relevant to the genes in the cancer data.
goentrez = GO.list[genenames] ## pull out pathways that are relevant
goentrez = goentrez[which(!is.na(names(goentrez)))] #remove NAs
genenames = names(goentrez) #update genenames
goentrez = unique(unlist(goentrez)) #get unique entrezID

pathnames = gopath[goentrez] ## pull out necessary paths using the entrezID
pathnames = unlist(pathnames)

tmp = matrix(nrow = length(genenames), ncol = length(pathnames))
rownames(tmp) = genenames
colnames(tmp) = pathnames

for(i in 1:length(rownames(tmp))) {
    v = colnames(tmp) %in% lookup(rownames(tmp)[i], "GO")
    if(length(v) > 0) {
        tmp[rownames(tmp)[i],] = v
    }
}

goTFmat = ifelse(tmp == F, 0, 1)

##remove pathways that have no genes involved (i.e. 0s)
keggTFmat = keggTFmat[,!apply(keggTFmat==0, 2, all)]
goTFmat = goTFmat[,!apply(goTFmat ==0, 2, all)]
reactomeTFmat = reactomeTFmat[,!apply(reactomeTFmat==0, 2, all)]

##save the files as txt
dput(keggTFmat, file = 'kegggenebypath.txt')
dput(reactomeTFmat, file = 'reactomegenebypath.txt')
dput(goTFmat, file = 'gogenebypath.txt')

# Pathway enrichment analysis function
# It should take in DEGs and pull out relevant data from a pre-made
# gene-by-pathway matrix.
# deg should be in the top table format (containing only significant genes)
# db is the type of database you want to use in string format (i.e. 'KEGG',
# 'reactome', or 'GO')

pathenrich = function(deg, db, gn = genenames) {
    #check which database you're using
    if (db == "KEGG") {
        db = keggTFmat
    } else if (db == "GO") {
        db = goTFmat
    } else {
        db = reactomeTFmat
    }

    ##pull out relevant genes from the TF mat of the database
    #mat = db[which(rownames(db) %in% rownames(deg)),]
    ##remove columns with 0s (i.e. pathways with no DEGs)
    #mat = mat[,!apply(mat==0, 2, all)]

    # matrix to dump the p-values of the pathway enrichment
    v = matrix(nrow = ncol(db), ncol = 1, dimnames = list(colnames(db),'p.value'))
    x = gn %in% rownames(deg)
    for (i in 1:ncol(db)){
        y = gn %in% rownames(db)[db[,i]==1]
        if(length(table(y)) > 1 ) { ##do Fisher's test only if there is a '1' inthere
            f = fisher.test(table(x,y))
            v[i,1] = f$p.value
        } else {
            v[i,1] = NA
        }
    }

    return(v)
}

#set genenames to all 20501 genes in the  cancer data
genenames = rownames(BLCAmat)

##This function does path enrichment analysis n times on a single cancer type
# dat is the cancer data, sample_dat is the sample distribution in original data (normal/overweight/obese), gn is a vector containing all the genes in the cancer data (i.e. rownames(BLCAmat)), n is the number of analyses to be carried out.
npathenrich = function(dat, sample_dat, gn = genenames, n = 100, db = "KEGG"){
    for(i in 1:n) {
        #pick samples randomly
        group = sample(seq(sum(sample_dat)),sample_dat[3], replace = F)
        #make a vector out of the randomly chosen samples
        group = seq(sum(sample_dat)) %in% group
        #make model matrix using the randomly chosen samples.
        design = model.matrix(~group)
        #normalise the data
        x = normVoom(dat, design)
        #make top table from the data
        tt = make_tt(x, design)
        # pull out DEGs
        deg = pull_deg(tt)

        tmppath = pathenrich(deg, db, gn = genenames)

        if(i == 1) {
            mat = tmppath
        } else if (i > 1){
            mat = cbind(mat, tmppath)
        }
    }
    return(mat)
}


#############################################################################
# start pathway enrichment analysis:
#############################################################################

#begin by doing the pathway enrichment analysis on the ob vs lean/ov on each
# cancer type

files = gsub('mat','bmi',files)

originalpathlist = originalList

## do pathway enrichment analysis on each cancer type
originalpathlist = lapply(seq_along(originalpathlist), function(x) {
                       group = get(files[x])[,4]
                       group = ifelse(group == 'obese', 'obese', 'normal/overweight')
                       model = model.matrix(~group)
                       dat = originalpathlist[[x]]
                       dat = normVoom(dat, model)
                       tt = make_tt(dat, model)
                       deg = pull_deg(tt)
                       path = pathenrich(deg, db = "GO", gn = rownames(dat))
})

## adjust the p-values using FDR:
originalpathlist = lapply(originalpathlist, function(x) {
                        p.adj = p.adjust(x, method = 'BH', n = nrow(x))
                        x = cbind(x, p.adj)
})

## pull out the names of enriched paths:
enrpaths = lapply(originalpathlist, function(x) {
                  ind = which(x[,2] <= 0.05)
                  rownames(x)[ind]
})

names(enrpaths) = gsub('bmi', '', files)

## ready to do the pathway enrichment analysis:

set.seed(1) ##set the seed for reproducibility

#make copy just in case
ntestrun = originalList

# repeat the test run, but with npathenrich function instead of pathenrich:
ntestrun = lapply(seq_along(ntestrun ), function(x) {
                          dat = ntestrun [[x]]
                          path = npathenrich( dat, sample_dat = samples[x,] ,n = 10, db = "GO", gn = rownames(dat))
})

ntestrun2 = ntestrun

ntestrun2 = lapply(ntestrun2 , function(x) {
                   tmp = apply(x, 2, function(y) p.adjust(y, method = 'BH', n = length(y)))
                   colnames(tmp) = paste('adj.', colnames(tmp), sep = '')
                   return(tmp)
})

## So, the iteration for the Fisher's test-based enrichment analysis is working
## Now I'll have to look into a rank-based enrichment analysis

###############################################################################
# Rank-based enrichment analysis
###############################################################################

## test out the camera function on the BLCA dataset
test = originalList[[1]]

## Make a design matrix:
group = BLCAbmi[,4]
group = ifelse(group == 'obese', 'obese', 'normal/overweight')
design = model.matrix(~group)

## normalise the data with normVoom()
test = normVoom(test, design)$E

## need to make a list of indices (that describe what genes are in each pathway) for the camera function:
goind = split(goTFmat, rep(1:ncol(goTFmat), each = nrow(goTFmat)))
names(goind) = colnames(goTFmat)

## lapply the list so that it contains the genenames involved in the pathway
## This list will be used as the index in the camera function
genenames = rownames(test)
goind = lapply(goind, function(x) genenames[which(x == 1)])

## use the camera function in limma package
test = camera(test, goind, design, use.ranks=T)

## use the roast function to see if it does any better:
test2 = roast(test, goind, design)

###############################################################################
# try camera on all the cancer types
###############################################################################

files = gsub('mat', 'bmi', files)

test = originalList

# default starting index list for lapply
goind = split(goTFmat, rep(1:ncol(goTFmat), each = nrow(goTFmat)))
names(goind) = colnames(goTFmat)

test = lapply(seq_along(test), function(x) {
              # make model matrix from bmi info
              group = get(files[x])
              group = group[,4]
              group = ifelse(group == 'obese', 'obese', 'normal/overweight')
              design = model.matrix(~group)

              # normalise RNA-seq data:
              y = test[[x]]
              #print(dim(y))
              #print(count)

              y = normVoom(y, design)$E

              genenames = rownames(y)

              # make index list for pathway genes:
              ind = goind
              ind = lapply(ind, function(z) genenames[which(z == 1)])

              # pathway enrichment:
              y = camera(y, goind, design, contrast = 2, use.ranks=T)
})


###############################################################################
# try pathway enrichment analysis with all the samples from different cancer
# types on the 'same scale'.
###############################################################################

## normalise the samples and put them on the same scale:

## Combine all the samples from different cancer types:
allsamples = cbind(BLCAmat, CESCmat, COADmat, KIRPmat, LIHCmat, READmat, SKCMmat, UCECmat)
allbmi = rbind(BLCAbmi, CESCbmi, COADbmi, KIRPbmi, LIHCbmi, READbmi, SKCMbmi, UCECbmi)
allbmi = allbmi[,4]

## normalise using the bmi data:
group = ifelse(allbmi == 'obese', 'obese', 'normal/overweight')
design = model.matrix(~group)

allsamplesnorm = normVoom(allsamples, design)$E

## split the normalised data back into separate cancer types:
files = gsub('bmi', 'mat', files)
varnames = paste('scaled',files, sep='')
for(i in 1:length(varnames)) {
    sample = colnames(get(files[i]))
    assign(varnames[i],allsamplesnorm[,sample])
}

scaledtest = list(scaledBLCAmat, scaledCESCmat, scaledCOADmat, scaledKIRPmat, scaledLIHCmat, scaledREADmat, scaledSKCMmat, scaledUCECmat)

files = gsub('mat', 'bmi', files)

scaledtest = lapply(seq_along(scaledtest), function(x) {
              # make model matrix from bmi info
              group = get(files[x])
              group = group[,4]
              group = ifelse(group == 'obese', 'obese', 'normal/overweight')
              design = model.matrix(~group)

              # normalise RNA-seq data:
              y = scaledtest[[x]]

              genenames = rownames(y)

              # make index list for pathway genes:
              ind = goind
              ind = lapply(ind, function(z) genenames[which(z == 1)])

              # pathway enrichment:
              y = camera(y, goind, design, contrast = 2, use.ranks=T)
})














###############################################################################
## Gatza pathway metagene analysis

## NOTE: From previous results on transformation matrix and SVD, standardised
## data is used.

# List of genes related to/representing each of the Gatza pathway:
checkgene = c('AKT1', 'CTNNB1', 'E2F1', 'EGFR', 'ESR1', 'ERBB2', 'IFNA1',
			  'IFNG', 'MYC', 'TP53', 'TP63', 'PIK3CA', 'PGR', 'HRAS', 'SRC',
			  'STAT3', 'TGFB1', 'TNF')

## Get rank-based and probit metagene from RMA normalised data:
flip = rep(F, 18) # Don't flip any metagene
gtrmameta       = make_multi_meta(gtrmastd, paths, flip = flip, raw = T)
gtrmarankmeta   = apply(gtrmameta, 1, function(x) rank(x)/length(x))
gtrmaprobitmeta = apply(gtrmameta, 1, function(x) pnorm(scale(x)))

## Repeat in MAS5 normalised data:
gtmasmeta = make_multi_meta(gtmasstd, paths, flip = flip, raw = T)
gtmasrankmeta   = apply(gtmasmeta, 1, function(x) rank(x)/length(x))
gtmasprobitmeta = apply(gtmasmeta, 1, function(x) pnorm(scale(x)))

## Look at the difference between RMA vs. MAS5 and rank-based vs. probit:
pdf('pdf/gt_mg_results/gt_meta_comparison.pdf')
for (i in 1:ncol(gtrmarankmeta)) {
	rank1   = gtrmarankmeta[,i]
	rank2   = gtmasrankmeta[,i]
	probit1 = gtrmaprobitmeta[,i]
	probit2 = gtmasprobitmeta[,i]
	txt = colnames(gtrmarankmeta)[i]
	txt = toupper(gsub('_probes', '', txt))
	main = paste("RMA vs. MAS5 (ranked metagene)", txt, sep = ' - ')
	plot(rank1, rank2, pch=20, main=main, xlab='Metagene values from RMA data', ylab='Metagene values from MAS5 data')
	main = gsub('ranked', 'probit', main)
	plot(probit1, probit2, pch=20, main=main, xlab='Metagene values from RMA data', ylab='Metagene values from MAS5 data')
	main = paste("Ranked vs. probit metagenes (RMA data)", txt, sep = ' - ')
	plot(rank1, probit1, pch=20, main=main, xlab='Ranked metagene values', ylab='Probit metagene values')
	main = gsub('RMA', 'MAS5', main)
	plot(rank2, probit2, pch=20, main=main, xlab='Ranked metagene values', ylab='Probit metagene values')
}
dev.off()

















## TODO: work from here
print()

###############################################################################
## Gatza and BMI metagene in ICGC

## make metagene with Gatza pathways in ICGC data samples:
gatzametalist = list()
for (i in 1:length(paths)) {
	pdfname = gsub('probes', 'meta.pdf', paths[i])
	pdfname = paste('pdf/', pdfname, sep='')
	pdf(pdfname)
	gene = get(paths[i])
	mat = crsymmat[gene,]
	mat = t(apply(mat, 1, function(x) (x-mean(x))/sd(x)))
	mat[mat < -3] = -3
	mat[mat > 3] = 3
	tmpsvd = svd(mat)
	trmat = diag(1/tmpsvd$d) %*% t(tmpsvd$u)
	tmpmeta = vector()
	for (j in 1:length(cancertypes)) {
		cancer = t(get(cancertypes[j]))
		mat = cancer[gene,]
		mat = standardisedata(mat)
		bmi = get(paste(cancertypes[j], 'bmi', sep=''))
		tmpsvd = t(trmat %*% mat)
		tmpsvd = tmpsvd[,1]
		tmpsvd = rank(tmpsvd)/length(tmpsvd)
		if (paths[i] %in% mgflip) {
			tmpsvd = 1-tmpsvd
		}
		main = gsub('probes', '', paths[i])
		main = toupper(main)
		main = paste(main, 'metagene (')
		main = paste(main, cancertypes[j], sep = '')
		main = paste(main, ')', sep = '')
		metaplot3(mat, tmpsvd, bmi, name = main)
		tmpmeta = c(tmpmeta, tmpsvd)
	}
	gatzametalist[[i]] = tmpmeta
	dev.off()
}
gatzametalist = as.data.frame(gatzametalist)
colnames(gatzametalist) = paths

## stick the BMI metagene and Gatza metagene together:
allmetagene = cbind(bmimetalist, gatzametalist)
allmetagene = t(data.matrix(allmetagene))

## need to make a combined clinical data for all the cancer types
x = c('bcrpatientbarcode', 'tumortissuesite', 'weight', 'height')
allclin = BLCAclin[,x]
for(i in 2:length(cancertypes)) {
	txt = paste(cancertypes[i], 'clin', sep='')
	clin = get(txt)
	clin = clin[,x]
	allclin = rbind(allclin, clin)
}
rownames(allclin) = allclin$bcrpatientbarcode
allclin = allclin[colnames(allmetagene),-1]
allclin$weight = as.numeric(as.character(allclin$weight))
allclin$height = as.numeric(as.character(allclin$height))
allclin[,4] = allclin$weight/((allclin$height/100)^2)
colnames(allclin)[4] = 'bmi'
x = allclin[,4]
for (i in 1:length(x)) {
	if(x[i] >= 30){
		x[i] = "obese"
	} else if(x[i] < 25){
		x[i] = "normal"
	} else {
		x[i] = "overweight"
	}
}
allclin[,5] = x
colnames(allclin)[5] = 'bmiStatus'

## prepare colours for cancer types:
brewercol = brewer.pal(8, "Set1")
cancertypecol = c()
type = c()
for (i in 1:length(cancertypes)) {
	type = c(type, rep(cancertypes[i], nrow(get(cancertypes[i]))))
	cancertypecol = c(cancertypecol, rep(brewercol[i], nrow(get(cancertypes[i]))))
}
allclin[,1] = type

## prepare colours for BMI status:
brewercol = brewer.pal(3, "Set1")
bmicol = allclin[,4]
for (i in 1:length(bmicol)) {
	if(bmicol[i] >= 30){
		bmicol[i] = brewercol[1]
	} else if(bmicol[i] < 25){
		bmicol[i] = brewercol[2]
	} else {
		bmicol[i] = brewercol[3]
	}
}

#bmicol = ifelse(allclin[,5] == "obese", brewercol[1], brewercol[2])

## make a vector of the order of the samples (by BMI values), per cancer type.
bmiord = vector()
bmivalcol = vector()
count = 0
for (i in 1:length(cancertypes)) {
	type = cancertypes[i]
	cl = allclin[which(allclin$tumortissuesite == type),]
	cl = cl$bmi
	ord = order(rank(cl)) + count
	col = bluered(length(ord))[rank(cl)]

	bmiord = c(bmiord, ord)
	bmivalcol = c(bmivalcol, col)
	count = length(bmiord)
}

allcol = rbind(BMI=bmivalcol, BMIStatus=bmicol, CancerType=cancertypecol)

pdf(file='pdf/gatzabmimeta.pdf',width=14, height=7)
heatmap.2x(allmetagene, scale='none', trace='none', col=bluered(2000), ColSideColors=allcol)
heatmap.2x(allmetagene[, 1:1872], scale='none', trace='none', col=bluered(2000), ColSideColors=allcol, Colv=F)
heatmap.2x(allmetagene[, bmiord], scale='none', trace='none', col=bluered(2000), ColSideColors=allcol[,bmiord], Colv=F)
dev.off()

pdf('pdf/allmetacor.pdf')
x = cor(t(allmetagene), method='pearson')
heatmap.2(x, trace = 'none', scale='none', col='bluered', main="Pearson correlation")
x = cor(t(allmetagene), method='spearman')
heatmap.2(x, trace = 'none', scale='none', col='bluered', main='Spearman correlation')
x = cor(gatzametalist, method='spearman')
heatmap.2(x, trace = 'none', scale='none', col='bluered', main='Spearman correlation')
dev.off()

# Plot heatmaps of all the metagene for each cancer type
pdf(file='pdf/cancersepallmeta.pdf', width=7, height=7)
for (i in 1:length(cancertypes)) {
	type = cancertypes[i]
	ind = which(allclin$tumortissuesite == type)
	samples = rownames(allclin)[ind]
	mat = allmetagene[,samples]
	heatmap.2(x=mat, scale='none', trace='none', col='bluered', main=cancertypes[i])
	heatmap.2x(allmetagene[, bmiord[ind]], scale='none', trace='none', col=bluered(2000), ColSideColors=allcol[c("BMI", "BMIStatus"),bmiord[ind]], Colv=F)
}
dev.off()








print('hello world!') ## command to stop nvimr to go down to the bottom






























































