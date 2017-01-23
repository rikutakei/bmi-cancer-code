###############################################################################

                        ## Final R code for the project

###############################################################################

## set working directory:
setwd('/home/riku/Documents/masters/data/task15')

## Load previously saved workspace:
# load(".RData")

## load functions:
# source('/home/riku/Documents/codes/bmi-cancer-code/task15/functions_final.R')
source('~/scratch/task15/functions_final.R')

source('https://bioconductor.org/biocLite.R')
install.packages(GMD)
install.packages(colorRamps)
install.packages(data.table)
install.packages(gplots)
install.packages(lattice)
install.packages(sva)

## heatmap stuff
install.packages(devtools)
devtools::install_github("TomKellyGenetics/heatmap.2x", ref="master")

## libraries for data wrangling
install.packages(WGCNA)
install.packages(dplyr)
install.packages(hgu133a.db)
install.packages(tidyr)

## libraries for DEG analysis
install.packages(edgeR)
install.packages(DESeq)
install.packages(affy)
install.packages(limma)

## libraries for pathway analysis
install.packages(GO.db)
install.packages(KEGG.db)
install.packages(org.Hs.eg.db)
install.packages(reactome.db)


## load libraries:
# library(mclust)
library(RColorBrewer)
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
library(edgeR)
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
# source("shortcut.R")

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

# Save data:
save(crrma, crmas, crmstdrma, crstdmas, crclin, file = 'shortcutData/cr_data.txt')

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

# Save data:
save(fmrma, fmmas, fmstdrma, fmstdmas, fmclin, 'shortcutData/fm_data.txt')

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

# Save data:
save(crisrma, crismas, crisstdrma, crisstdmas, crisclin, 'shortcutData/cris_data.txt')

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

# Save data:
save(gatzarma, gatzamas, gtrmastd, gtmasstd, 'shortcutData/gatza_data.txt')

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
	mat = read.table(cancerfiles[i], sep = '\t', header=T)
	mat = tbl_df(mat)
	dup = duplicated(paste(mat$submitted_sample_id,mat$gene_id,sep=''))
	mat = mat[!dup,c('submitted_sample_id','gene_id','raw_read_count')]
	genes = unique(mat$gene_id)
	mat  = spread(mat,submitted_sample_id,raw_read_count)
	rownames(mat) = genes
	mat  = data.matrix(mat)
	mat = mat[-1,-1]
	mat = icgc_to_tcga(t(mat))
    if (length(rownames(mat)) > length(unique(rownames(mat)))) {
        mat = collapseRows(mat, unique(rownames(mat)), unique(rownames(mat)))
        mat = mat$datETcollapsed
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

	## Remove samples that are not in the matuence data:
	mat = mat[which(rownames(mat) %in% rownames(clin)),]
	clin = clin[which(rownames(clin) %in% rownames(mat)),]

	## Calculate BMI and BMI status:
	bmi = as.numeric(as.vector(clin$weight))/((as.numeric(as.vector(clin$height))/100)^2)
	bmiStatus = bmi
	bmiStatus[bmiStatus >= 30] = "obese"
	bmiStatus[bmiStatus < 25] = "normal"
	bmiStatus[(bmiStatus < 30) & (bmiStatus >= 25)] = "overweight"

	## Add BMI information to the clinical data:
	clin = cbind(clin, bmi, bmiStatus)

	## assign variable names:
	assign(paste(files[i], 'raw', sep = ''), mat)
	assign(paste(files[i], 'clin', sep=''), clin)
}

for (i in 1:length(files)) {
	txt = paste(files[i], 'raw', sep='')
	mat = t(get(txt))
	mat = standardise_data(mat)
	txt = paste(files[i], 'std', sep='')
	assign(txt, t(mat))
}

# Save data:
save(BLCAraw, BLCAstd, BLCAclin, CESCraw, CESCstd, CESCclin, COADraw, COADstd, COADclin, KIRPraw, KIRPstd, KIRPclin, LIHCraw, LIHCstd, LIHCclin, READraw, READstd, READclin, SKCMraw, SKCMstd, SKCMclin, UCECraw, UCECstd, UCECclin, 'shortcutData/icgc_data.txt')

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

## Save the data:
save(crsymrma, crsymmas, crsymrmastd, crsymmasstd, fmsymrma, fmsymmas, fmsymrmastd, fmsymmasstd, crissymrma, crissymmas, crissymrmastd, crissymmasstd, gatzasymrma, gatzasymmas, gatzasymrmastd, gatzasymmasstd, 'shortcutData/symmat_data.txt')

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

#Import Human Reactome pathways
reactome.list = as.list(reactomeEXTID2PATHID)
reactome.list = reactome.list[order(as.numeric(names(reactome.list)))] #Sort the names in the list:
reactome.list = reactome.list[which(names(reactome.list) %in% names(SYMBOL.list))] #get paths that have gene symbols:
names(reactome.list) = unlist(SYMBOL.list[names(reactome.list)]) #rename the entrez gene ID into gene symbol
reactomepath = as.list(reactomePATHID2NAME) #mapping reactome path IDs to human read pathway name
reactomepath = reactomepath[grep('Homo sapiens',reactomepath)] #pull out all human-related pathways
reactomepath = lapply(reactomepath, function(x) gsub('Homo sapiens: ', '', x)) #cut out the 'Homo sapiens: ' bit so it's only the pathway names.
reactomepath = lapply(reactomepath, function(x) return(x[1])) #cut out the 'Homo sapiens: ' bit so it's only the pathway names.

## Save pathway databases (raw):
save(keggpath, reactomepath, gopath, file = 'rawpath_data.txt')

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
	txt = paste(cancertypes[i], 'raw', sep='')
	mat = get(txt)
	mat = mat[,commongenes]
	mat = t(mat) ## make it genes by samples
	txt = paste(cancertypes[i], 'cr', sep='')
	assign(txt, mat)
}

## check if it has the right number of genes:
dim(BLCAcr) # 644 genes

## log10 and standardise data:
for (i in 1:length(cancertypes)) {
	txt = paste(cancertypes[i], 'cr', sep='')
	mat = get(txt)
	mat = standardise_data(mat)
	txt = paste(txt, 'std', sep='')
	assign(txt, mat)
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
	mat = t(get(txt))
	mat = mat[tmpgenes,]
	txt = paste(cancertypes[i], 'fm', sep='')
	assign(txt, mat)
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

# source('path_enrich.R')

# Make a gene-by-pathway matrix to use it in the pathway enrichment analysis

# KEGG matrix
genenames = colnames(BLCAraw) #Get all the gene names in the cancer data
save(genenames, 'shortcutData/icgcgenenames.txt') # save the genenames in separate data
pathnames = unique(unlist(keggpath)) #Get all the pathways in the keggpath list
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

## Save pathway databases:
save(keggTFmat, reactomeTFmat, goTFmat, file = 'tfmat_data.txt')

###############################################################################
# Rank-based enrichment analysis (each cancer type separately)

## Use the path_enrich() function to find the enriched pathways between the
## obese and non-obese samples in the ICGC cancer data:
reslist = list()
for (i in 1:length(cancertypes)) {
	mat = t(get(paste(cancertypes[i], 'raw', sep = '')))
	clin = get(paste(cancertypes[i], 'clin', sep = ''))
	tmplist = list()
	tmplist[[1]] = path_enrich(mat, clin, keggTFmat)
	tmplist[[2]] = path_enrich(mat, clin, reactomeTFmat)
	tmplist[[3]] = path_enrich(mat, clin, goTFmat)
	names(tmplist) = c('KEGG', 'Reactome', 'GO')
	reslist[[i]] = tmplist
}
names(reslist) = cancertypes

###############################################################################
# Rank-based enrichment analysis (all cancer types combined)

## First try with Voom normalising individual cancer types separately, then
## batch correct with ComBat:
pdf(file='pdf/voom_pdf/per_cancer_voom.pdf', width=7, height=7)
# pdf(file='pdf/voom_pdf/per_cancer_voom(nullGroup).pdf', width=7, height=7)
# pdf(file='pdf/voom_pdf/per_cancer_voom(5percent).pdf', width=7, height=7)
samplename = rownames(BLCAraw)
batchinfo = rep('BLCA', nrow(BLCAraw))
group = ifelse(BLCAclin$bmiStatus == 'obese', 'obese', 'non-obese')
allcancer = voom_norm(t(BLCAraw), group, plot = T)
# allcancer = voom_norm(t(BLCAraw), group = NULL, plot = T)
# allcancer = voom_norm(t(BLCAraw), group, dropout = 0.05, plot = T)
allcancer = tbl_df(t(allcancer))
for (i in 2:length(cancertypes)) {
	mat = get(paste(cancertypes[i], 'raw', sep=''))
	clin = get(paste(cancertypes[i], 'clin', sep=''))
	samplename = c(samplename, rownames(mat))
	batchinfo = c(batchinfo, rep(cancertypes[i], nrow(mat)))
	group = ifelse(clin$bmiStatus == 'obese', 'obese', 'non-obese')
	mat = voom_norm(t(mat), group, plot = T)
	# mat = voom_norm(t(mat), group = NULL, plot = T)
	# mat = voom_norm(t(mat), group, dropout = 0.05, plot = T)
	mat = tbl_df(t(mat))
	allcancer = full_join(allcancer, mat)
}
dev.off()
allcancer = as.matrix(allcancer)
rownames(allcancer) = samplename
allcancer = t(allcancer)
allcancer = allcancer[which(complete.cases(allcancer)),]
dim(allcancer) ## 20494 genes by 1872 samples
# dim(allcancer) ## 13630 genes by 1872 samples after 5% dropout

batchinfo = cbind(samplename, batchinfo)
mod = model.matrix(~1, data = as.data.frame(batchinfo))
allcancer = ComBat(dat = allcancer, batch = batchinfo[,2], mod = mod)

## Try Voom normalising all of the cancer types after merging, then batch
## correct with ComBat:
samplename = rownames(BLCAraw)
batchinfo = rep('BLCA', nrow(BLCAraw))
bmidata = as.character(BLCAclin$bmiStatus)
allcancer2 = tbl_df(BLCAraw)
for (i in 2:length(cancertypes)) {
	mat = get(paste(cancertypes[i], 'raw', sep=''))
	clin = get(paste(cancertypes[i], 'clin', sep=''))
	samplename = c(samplename, rownames(mat))
	batchinfo = c(batchinfo, rep(cancertypes[i], nrow(mat)))
	bmidata = c(bmidata, as.character(clin$bmiStatus))
	mat = tbl_df(mat)
	allcancer2 = full_join(allcancer2, mat)
}
allcancer2 = as.matrix(allcancer2)
rownames(allcancer2) = samplename
allcancer2 = t(allcancer2)
allcancer2 = allcancer2[which(complete.cases(allcancer2)),]
dim(allcancer2) ## 20494 genes by 1872 samples

group = ifelse(bmidata == 'obese', 'obese', 'non-obese')
pdf(file='pdf/voom_pdf/all_cancer_voom.pdf', width=7, height=7)
# pdf(file='pdf/voom_pdf/all_cancer_voom(nullGroup).pdf', width=7, height=7)
# pdf(file='pdf/voom_pdf/all_cancer_voom(5percent).pdf', width=7, height=7)
allcancer2 = voom_norm(allcancer2, group, plot = T)
# allcancer2 = voom_norm(allcancer2, group = NULL, plot = T)
# allcancer2 = voom_norm(allcancer2, group, dropout = 0.05, plot = T)
dev.off()
# dim(allcancer2) ## 19469 genes by 1872 samples after 5% dropout

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
batchinfo = cbind(samplename, batchinfo)
mod = model.matrix(~1, data = as.data.frame(batchinfo))
allcancer2 = ComBat(dat = allcancer2, batch = batchinfo[,2], mod = mod)

###############################################################################
## See what the difference in the metagene are between per-cancer normalised
## data set and the all-cancer normalised data set

## Use Creighton's metagene to test the correlation between the two data sets
commongenes = mapIds(hgu133a.db, keys = crobsgenes, column = 'SYMBOL', keytype = "PROBEID", multiVals = 'first')
commongenes = commongenes[which(!is.na(commongenes))]
commongenes = unique(commongenes)
commongenes = commongenes[which(commongenes %in% rownames(allcancer))]

## Get Creighton metagene from each cancertype (voom normalised individually):
percancermeta = list()
percancermat = list()
for (i in 1:length(cancertypes)) {
	mat = get(paste(cancertypes[i], 'raw', sep=''))
	clin = get(paste(cancertypes[i], 'clin', sep=''))
	group = ifelse(clin$bmiStatus == 'obese', 'obese', 'non-obese')
	mat = voom_norm(t(mat), group)
	percancermat[[i]] = mat
	mat = mat[commongenes,]
	meta = make_metagene(mat, flip = F, raw = T)
	percancermeta[[i]] = meta
}

## Get Creighton metagene from the combined cancer data set:
samplename = rownames(BLCAraw)
bmidata = as.character(BLCAclin$bmiStatus)
tmpallcancer = tbl_df(BLCAraw)
for (i in 2:length(cancertypes)) {
	mat = get(paste(cancertypes[i], 'raw', sep=''))
	clin = get(paste(cancertypes[i], 'clin', sep=''))
	samplename = c(samplename, rownames(mat))
	bmidata = c(bmidata, as.character(clin$bmiStatus))
	mat = tbl_df(mat)
	tmpallcancer = full_join(tmpallcancer, mat)
}
tmpallcancer = as.matrix(tmpallcancer)
rownames(tmpallcancer) = samplename
tmpallcancer = t(tmpallcancer)
tmpallcancer = tmpallcancer[which(complete.cases(tmpallcancer)),]
dim(tmpallcancer) ## 20494 genes by 1872 samples

group = ifelse(bmidata == 'obese', 'obese', 'non-obese')
tmpallcancer = voom_norm(tmpallcancer, group)

crtmpallcancer = tmpallcancer[commongenes,]
allcancermeta = make_metagene(crtmpallcancer, flip = F, raw = T)
names(allcancermeta) = samplename

## Split the metagene into 8 different cancertypes:
tmp = list()
for (i in 1:length(cancertypes)) {
	mat = get(paste(cancertypes[i], 'raw', sep=''))
	tmp[[i]] = allcancermeta[which(names(allcancermeta) %in% rownames(mat))]
}
allcancermeta = tmp

## Correlation of the raw metagenes:
tmp = percancermeta
tmp2 = allcancermeta

rawmetacor = list()
for (i in 1:length(tmp)) {
	correlation = cor(tmp[[i]], tmp2[[i]], method = 'spearman')
	rawmetacor[[i]] = correlation
}

## Correlation of the ranked metagenes:
tmp = lapply(percancermeta, function(x) rank(x)/length(x))
tmp2 = lapply(allcancermeta, function(x) rank(x)/length(x))

rankmetacor = list()
for (i in 1:length(tmp)) {
	correlation = cor(tmp[[i]], tmp2[[i]], method = 'spearman')
	rankmetacor[[i]] = correlation
}

###############################################################################
## Try the metagene comparison between per-cancer and combined cancer with IFNG
## signature, and split the samples based on the IFNG signature (high vs. low)

## Convert gene probe IDs into symbols:

commongenes = mapIds(hgu133a.db, keys = ifng_probes, column = 'SYMBOL', keytype = "PROBEID", multiVals = 'first')
commongenes = commongenes[which(!is.na(commongenes))]
commongenes = unique(commongenes)
commongenes = commongenes[which(commongenes %in% rownames(allcancer))]

## Get Creighton metagene from each cancertype (voom normalised individually):
percancermeta = list()
percancermat = list()
for (i in 1:length(cancertypes)) {
	mat = get(paste(cancertypes[i], 'raw', sep=''))
	mat = voom_norm(t(mat), group = NULL)
	percancermat[[i]] = mat
	mat = mat[commongenes,]
	meta = make_metagene(mat, flip = F, raw = T)
	percancermeta[[i]] = meta
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



## Get Creighton metagene from the combined cancer data set:
samplename = rownames(BLCAraw)
bmidata = as.character(BLCAclin$bmiStatus)
tmpallcancer = tbl_df(BLCAraw)
for (i in 2:length(cancertypes)) {
	mat = get(paste(cancertypes[i], 'raw', sep=''))
	samplename = c(samplename, rownames(mat))
	mat = tbl_df(mat)
	tmpallcancer = full_join(tmpallcancer, mat)
}
tmpallcancer = as.matrix(tmpallcancer)
rownames(tmpallcancer) = samplename
tmpallcancer = t(tmpallcancer)
tmpallcancer = tmpallcancer[which(complete.cases(tmpallcancer)),]
dim(tmpallcancer) ## 20494 genes by 1872 samples

tmpallcancer = voom_norm(tmpallcancer, group = NULL)

ifngtmpallcancer = tmpallcancer[commongenes,]
allcancermeta = make_metagene(ifngtmpallcancer, flip = F, raw = T)
names(allcancermeta) = samplename

## Split the metagene into 8 different cancertypes:
tmp = list()
for (i in 1:length(cancertypes)) {
	mat = get(paste(cancertypes[i], 'raw', sep=''))
	tmp[[i]] = allcancermeta[which(names(allcancermeta) %in% rownames(mat))]
}
allcancermeta = tmp

## Correlation of the raw metagenes:
tmp = percancermeta
tmp2 = allcancermeta

rawmetacor = list()
for (i in 1:length(tmp)) {
	correlation = cor(tmp[[i]], tmp2[[i]], method = 'spearman')
	rawmetacor[[i]] = correlation
}

## Correlation of the ranked metagenes:
tmp = lapply(percancermeta, function(x) rank(x)/length(x))
tmp2 = lapply(allcancermeta, function(x) rank(x)/length(x))

rankmetacor = list()
for (i in 1:length(tmp)) {
	correlation = cor(tmp[[i]], tmp2[[i]], method = 'spearman')
	rankmetacor[[i]] = correlation
}

###############################################################################
## Pathway enrichment analysis with per-cancer and combined-cancer normalised
## data:

## Load shortcut data:
# load('shortcutData/allcancer.txt')

## Make clinical data for all cancer types:
allbmi = c()
for (i in 1:length(cancertypes)) {
	clin = get(paste(cancertypes[i], 'clin', sep=''))
	bmi = as.character(clin$bmiStatus)
	allbmi = c(allbmi, bmi)
}
names(allbmi) = colnames(allcancer)
allbmi = as.data.frame(allbmi)
colnames(allbmi) = 'bmiStatus'

## Path enrichment analysis using camera:
cancermat = c('allcancer', 'allcancer2')
db = c('keggTFmat', 'reactomeTFmat', 'goTFmat')
enrichreslist = list()
for (i in 1:length(cancermat)) {
	mat = get(cancermat[i])
	tmplist = list()
	for (j in 1:length(db)) {
		tmpdb = get(db[j])
		tmplist[[j]] = path_enrich(mat, allbmi, tmpdb, norm = F)
	}
	names(tmplist) = db
	enrichreslist[[i]] = tmplist
}
names(enrichreslist) = c('per-cancer', 'combined-cancer')

###############################################################################
## Gatza pathway analysis
###############################################################################

## First, check the directionality of the Gatza pathway metagenes so that they
## produce similar results as their paper

# List of genes representing/related to the pathway:
checkgene = c('AKT1', 'CTNNB1', 'E2F1', 'EGFR', 'ESR1', 'ERBB2', 'IFNA1',
			  'IFNG', 'MYC', 'TP53', 'TP63', 'PIK3CA', 'PGR', 'HRAS', 'SRC',
			  'STAT3', 'TGFB1', 'TNF')

path_symbol = paste(paths, '_sym', sep='')
for (i in 1:length(paths)) {
	gene = get(paths[i])
	gene = mapIds(hgu133a.db, keys = gene, column = 'SYMBOL', keytype = "PROBEID", multiVals = 'first')
	gene = c(gene[which(!is.na(gene))], checkgene[i])
	gene = unique(gene)
	gene = gene[which(gene %in% rownames(gatzasymrma))]
	gene = gene[which(gene %in% rownames(gatzasymmas))]
	assign(path_symbol[i], gene)
}

## Make metagenes for each pathway genes (gene probe IDs) in data sets, then
## visualise it using the symbol matrix and pathway genes (gene symbols):

## NOTE: standardised data were used (experience from Creighton metagene
## results).

## RMA normalised data:

# Ranked metagenes:
mgfliprma = c(
			   # 'akt_probes',
			   'bcat_probes',
			   'e2f1_probes',
			   'egfr_probes',
			   # 'er_probes',
			   'her2_probes',
			   'ifna_probes',
			   'ifng_probes',
			   # 'myc_probes',
			   'p53_probes',
			   'p63_probes',
			   'pi3k_probes',
			   'pr_probes',
			   # 'ras_probes',
			   # 'src_probes',
			   'stat3_probes',
			   'tgfb_probes',
			   'tnfa_probes'
			  )
flip = paths %in% mgfliprma
meta = make_multi_meta(gtrmastd, paths, flip = flip, raw = T)
rownames(meta) = gsub('_probes', '', rownames(meta))
rmarankmeta = t(apply(meta, 1, function(x) convert_raw_meta(x, method = 'rank')))
corrank = cor(t(rmarankmeta), method = 'pearson')

pdf('pdf/gatza_pdf/gtrmastdrank.pdf')
recreate_gatza(gatzasymrmastd, rmarankmeta, path_symbol, checkgene, type = 'ranked')
dev.off()

# Probit metagenes:
mgflipprobit = c(
				 # 'akt_probes',
				 'bcat_probes',
				 'e2f1_probes',
				 'egfr_probes',
				 # 'er_probes',
				 'her2_probes',
				 'ifna_probes',
				 'ifng_probes',
				 # 'myc_probes',
				 'p53_probes',
				 'p63_probes',
				 'pi3k_probes',
				 'pr_probes',
				 # 'ras_probes',
				 # 'src_probes',
				 'stat3_probes',
				 'tgfb_probes',
				 'tnfa_probes'
				 )
flip = paths %in% mgflipprobit
meta = make_multi_meta(gtrmastd, paths, flip = flip, raw = T)
rownames(meta) = gsub('_probes', '', rownames(meta))
rmaprobitmeta = t(apply(meta, 1, function(x) convert_raw_meta(x, method = 'probit')))
corprobit = cor(t(probitmeta), method = 'pearson')
pdf('pdf/gatza_pdf/gtrmastdprobit.pdf')
recreate_gatza(gatzasymrmastd, rmaprobitmeta, path_symbol, checkgene, type = 'probit')
dev.off()

## Check the consistency between ranked and probit metagenes:
pdf(file='pdf/gatza_pdf/gatza_rma_meta_rank_vs_probit.pdf', width=7, height=7)
x1 = rmarankmeta
x2 = rmaprobitmeta
for (i in 1:nrow(x1)) {
	main = rownames(x1)[i]
	plot(x1[i,], x2[i,], main=main, xlab="Ranked metagene", ylab="Probit metagene", pch = 20)
}
dev.off()

## Repeat in MAS5 normalised data:

# Ranked metagenes:
mgflipmas = c(
			  'akt_probes',
			  'bcat_probes',
			  # 'e2f1_probes',
			  # 'egfr_probes',
			  # 'er_probes',
			  # 'her2_probes',
			  'ifna_probes',
			  'ifng_probes',
			  # 'myc_probes',
			  'p53_probes',
			  'p63_probes',
			  'pi3k_probes',
			  'pr_probes',
			  # 'ras_probes',
			  'src_probes',
			  'stat3_probes',
			  'tgfb_probes',
			  'tnfa_probes'
			  )
flip = paths %in% mgflipmas
meta = make_multi_meta(gtmasstd, paths, flip = flip, raw = T)
rownames(meta) = gsub('_probes', '', rownames(meta))
masrankmeta = t(apply(meta, 1, function(x) convert_raw_meta(x, method = 'rank')))
corrank = cor(t(masrankmeta), method = 'pearson')

pdf('pdf/gatza_pdf/gtmasstdrank.pdf')
recreate_gatza(gatzasymmasstd, masrankmeta, path_symbol, checkgene, type = 'ranked')
dev.off()

# Probit metagene:
mgflipmas = c(
			  'akt_probes',
			  'bcat_probes',
			  # 'e2f1_probes',
			  # 'egfr_probes',
			  # 'er_probes',
			  # 'her2_probes',
			  'ifna_probes',
			  'ifng_probes',
			  # 'myc_probes',
			  'p53_probes',
			  'p63_probes',
			  'pi3k_probes',
			  'pr_probes',
			  # 'ras_probes',
			  'src_probes',
			  'stat3_probes',
			  'tgfb_probes',
			  'tnfa_probes'
			  )
flip = paths %in% mgflipmas
meta = make_multi_meta(gtmasstd, paths, flip = flip, raw = T)
rownames(meta) = gsub('_probes', '', rownames(meta))
masprobitmeta = t(apply(meta, 1, function(x) convert_raw_meta(x, method = 'probit')))
corrank = cor(t(masprobitmeta), method = 'pearson')

pdf('pdf/gatza_pdf/gtmasstdprobit.pdf')
recreate_gatza(gatzasymmasstd, masprobitmeta, path_symbol, checkgene, type = 'probit')
dev.off()

## Check the consistency between ranked and probit metagenes:
pdf(file='pdf/gatza_pdf/gatza_mas_meta_rank_vs_probit.pdf', width=7, height=7)
x1 = masrankmeta
x2 = masprobitmeta
for (i in 1:nrow(x1)) {
	main = rownames(x1)[i]
	plot(x1[i,], x2[i,], main=main, xlab="Ranked metagene", ylab="Probit metagene", pch = 20)
}
dev.off()

## Compare the metagenes between RMA and MAS5 normalised data:
pdf(file='pdf/gatza_pdf/combmetaplot.pdf', width=7, height=7)
x1 = rmarankmeta
x2 = rmaprobitmeta
x3 = masrankmeta
x4 = masprobitmeta
for (i in 1:nrow(x1)) {
	txt = rownames(x1)[i]
	main = paste(txt, ' (rank)', sep='')
	plot(x1[i,], x3[i,], main=main, xlab="RMA", ylab="MAS5")
	main = paste(txt, ' (probit)', sep='')
	plot(x2[i,], x4[i,], main=main, xlab="RMA", ylab="MAS5")
}
dev.off()

###############################################################################
# Make Gatza pathway transformation matrix in Gatza data:

gtrmatransmat = list()
for (i in 1:length(paths)) {
	gene = get(paths[i])
	mat = gtrmastd[gene,]
	trans = make_trans_mat(mat)
	gtrmatransmat[[i]] = trans
}
names(gtrmatransmat) = paths

gtmastransmat = list()
for (i in 1:length(paths)) {
	gene = get(paths[i])
	mat = gtmasstd[gene,]
	trans = make_trans_mat(mat)
	gtmastransmat[[i]] = trans
}
names(gtmastransmat) = paths

###############################################################################
# Which data set is the best to make transformation matrices (both for Gatza pathways and obesity metagenes)?
# Look at the correlation between the metagene produced from transformation matrix and SVD

# Check if the function works ("Control"):
gt1 = compare_svd_trans(gtrmastd, paths, gtrmatransmat)
gt2 = compare_svd_trans(gtmasstd, paths, gtmastransmat)
gt3 = compare_svd_trans(gtmasstd, paths, gtrmatransmat)
gt4 = compare_svd_trans(gtrmastd, paths, gtmastransmat)
# Correlation should be 1 for both test1 and test2, since the transformation matrices were derived in this dataset
gt1$correlation
gt2$correlation

# Correlation is slightly off for test3 and test4 - probably due to the difference in the normalisation method
gt3$correlation
gt4$correlation

# Repeat for other data sets:
# Creighton data:
cr1 = compare_svd_trans(crstdrma, paths, gtrmatransmat)
cr2 = compare_svd_trans(crstdmas, paths, gtmastransmat)
cr3 = compare_svd_trans(crstdmas, paths, gtrmatransmat)
cr4 = compare_svd_trans(crstdrma, paths, gtmastransmat)

# FM data:
fm1 = compare_svd_trans(fmstdrma, paths, gtrmatransmat)
fm2 = compare_svd_trans(fmstdmas, paths, gtmastransmat)
fm3 = compare_svd_trans(fmstdmas, paths, gtrmatransmat)
fm4 = compare_svd_trans(fmstdrma, paths, gtmastransmat)

# Cris data:
cris1 = compare_svd_trans(crisstdrma, paths, gtrmatransmat)
cris2 = compare_svd_trans(crisstdmas, paths, gtmastransmat)
cris3 = compare_svd_trans(crisstdmas, paths, gtrmatransmat)
cris4 = compare_svd_trans(crisstdrma, paths, gtmastransmat)

# Summarize the correlation data (RMA with RMA transmat, MAS with MAS transmat):
txt = c('Gatza','FM', 'Creighton', 'Cris')
gtrmatransres = cbind(gt1$correlation,fm1$correlation, cr1$correlation, cris1$correlation)
colnames(gtrmatransres) = txt
gtmastransres = cbind(gt2$correlation,fm2$correlation, cr2$correlation, cris2$correlation)
colnames(gtmastransres) = txt

# Summarize the correlation data (MAS with RMA transmat, RMA with MAS transmat):
txt = c('Gatza', 'FM', 'Creighton', 'Cris')
gtopptransres1 = cbind(gt3$correlation, fm3$correlation, cr3$correlation, cris3$correlation)
colnames(gtopptransres1) = txt
gtopptransres2 = cbind(gt4$correlation,fm4$correlation, cr4$correlation, cris4$correlation)
colnames(gtopptransres2) = txt

## Plot a bar graph to compare each metagenes with one another:
pdf('pdf/gatza_pdf/gatza_cor_barplot.pdf')
col = brewer.pal(4, 'Paired')
for(i in 1:nrow(gtrmatransres)) {
	for(j in 1:length(txt)) {
		main = paste(txt[j], '(')
		main = paste(main, rownames(gtrmatransres)[i], sep='')
		main = paste(main, ')', sep='')
		x1 = abs(gtrmatransres[i, j])
		x2 = abs(gtopptransres2[i, j])
		x3 = abs(gtmastransres[i, j])
		x4 = abs(gtopptransres1[i, j])
		dat = c(x1, x2, x3, x4)
		names(dat) = c('RMA/RMA','RMA/MAS5','MAS5/MAS5','MAS5/RMA')
		barplot(dat, ylim = c(0,1), main = main, col = col, ylab = 'Correlation (Absolute value)', xlab = 'Data/Transformation matrix')
	}
}
dev.off()

## From these barplots, some of the pathways do have variability between the
## SVD metagene and the transformation matrix metagene, where the transformation
## matrices are from the Gatza data.

###############################################################################
# Repeat the same thing, but with the BMI metagenes:
# (Make sure to re-run the Creighton DEG analysis to get the gene probe IDs for
# the obesity genes)

tmp = fmobsgenes
fmobsgenes = as.character(tmp$Probe)
obsname = c(allobsname, 'crobsgenes', 'fmobsgenes')

# make transformation matrix in Creighton's data:
obsrmatransmat = list()
for (i in 1:length(obsname)) {
	gene = get(obsname[i])
	mat = crstdrma[gene,]
	svd = svd(mat)
	trans = diag(1/svd$d) %*% t(svd$u)
	obsrmatransmat [[i]] = trans
	names(obsrmatransmat)[i] = obsname[i]
}

obsmastransmat = list()
for (i in 1:length(obsname)) {
	gene = get(obsname[i])
	mat = crstdmas[gene,]
	svd = svd(mat)
	trans = diag(1/svd$d) %*% t(svd$u)
	obsmastransmat [[i]] = trans
	names(obsmastransmat)[i] = obsname[i]
}

# Check if the function works ("Control"):
cr1 = compare_svd_trans(crstdrma, obsname, obsrmatransmat)
cr2 = compare_svd_trans(crstdmas, obsname, obsmastransmat)
cr3 = compare_svd_trans(crstdmas, obsname, obsrmatransmat)
cr4 = compare_svd_trans(crstdrma, obsname, obsmastransmat)

# Correlation should be 1 for both test1 and test2, since the transformation matrices were derived in this dataset
cr1$correlation
cr2$correlation

# Correlation is slightly off for test3 and test4 - probably due to the difference in the normalisation method
cr3$correlation
cr4$correlation

# Repeat for other data sets:
# Gatza data:
gt1 = compare_svd_trans(gtrmastd, obsname, obsrmatransmat)
gt2 = compare_svd_trans(gtmasstd, obsname, obsmastransmat)
gt3 = compare_svd_trans(gtmasstd, obsname, obsrmatransmat)
gt4 = compare_svd_trans(gtrmastd, obsname, obsmastransmat)

# FM data:
fm1 = compare_svd_trans(fmstdrma, obsname, obsrmatransmat)
fm2 = compare_svd_trans(fmstdmas, obsname, obsmastransmat)
fm3 = compare_svd_trans(fmstdmas, obsname, obsrmatransmat)
fm4 = compare_svd_trans(fmstdrma, obsname, obsmastransmat)

# Cris data:
cris1 = compare_svd_trans(crisstdrma, obsname, obsrmatransmat)
cris2 = compare_svd_trans(crisstdmas, obsname, obsmastransmat)
cris3 = compare_svd_trans(crisstdmas, obsname, obsrmatransmat)
cris4 = compare_svd_trans(crisstdrma, obsname, obsmastransmat)

# Summarize the correlation data (RMA with RMA transmat, MAS with MAS transmat):
txt = c('Creighton','Gatza', 'FM', 'Cris')
obsrmatransres = cbind(cr1$correlation,gt1$correlation, fm1$correlation, cris1$correlation)
colnames(obsrmatransres) = txt
obsmastransres = cbind(cr2$correlation,gt2$correlation, fm2$correlation, cris2$correlation)
colnames(obsmastransres) = txt

# Summarize the correlation data (MAS with RMA transmat, RMA with MAS transmat):
txt = c('Creighton','Gatza', 'FM', 'Cris')
obsopptransres1 = cbind(cr3$correlation,gt3$correlation, fm3$correlation, cris3$correlation)
colnames(obsopptransres1) = txt
obsopptransres2 = cbind(cr4$correlation,gt4$correlation, fm4$correlation, cris4$correlation)
colnames(obsopptransres2) = txt

## Plot a bar graph to compare each metagenes with one another:
pdf('pdf/gatza_pdf/obs_cor_barplot.pdf')
col = brewer.pal(4, 'Paired')
for(i in 1:nrow(obsrmatransres)) {
	for(j in 1:length(txt)) {
		main = paste(txt[j], '(')
		main = paste(main, rownames(obsrmatransres)[i], sep='')
		main = paste(main, ')', sep='')
		x1 = abs(obsrmatransres[i, j])
		x2 = abs(obsopptransres2[i, j])
		x3 = abs(obsmastransres[i, j])
		x4 = abs(obsopptransres1[i, j])
		dat = c(x1, x2, x3, x4)
		names(dat) = c('RMA/RMA','RMA/MAS5','MAS5/MAS5','MAS5/RMA')
		barplot(dat, ylim = c(0,1), main = main, col = col, ylab = 'Correlation (Absolute value)', xlab = 'Data/Transformation matrix')
	}
}
dev.off()

## Unlike the Gatza metagenes, it doesn't matter whether you make the metagene
## from the transformation matrix (from Creighton's data) or from SVD.
## Therefore, the transformation matrix for the obesity metagenes can be made in
## any dataset.

###############################################################################
## Does the normalisation method of the raw data and/or the normalisation
## method of the data in which the transformation matrix was derived from
## affect the metagene?

## TODO:










###############################################################################
## Quick metagene direction check for the obesity associated genes in Gatza data:

commongenes = table(c(rawobsgene,crolgene,resobsgene,rescrolgene, caobsgene,cacrolgene,caresobsgene,carescrolgene, crobsgenes))
commongenes = commongenes[commongenes==9]
commongenes = names(commongenes)

tmpname = c(allobsname, 'crobsgenes')

matheat = gtrmastd[commongenes,]
matheat[matheat < -3] = -3
matheat[matheat > 3] = 3

pdf('pdf/gatza_pdf/gtobrma.pdf')
for (i in 1:length(tmpname)) {
	genes = get(tmpname[i])
	tmp = gtrmastd[genes,]
	meta = make_metagene(tmp, flip = F, raw = F)
	ord = order(meta)
	heatmap.2(matheat[,ord], trace='none',scale='none', col='bluered', ColSideColors = bluered(length(meta))[rank(meta)][ord], Colv=NA, main=tmpname[i], cexRow=1.0)
}
dev.off()

## Obesity metagenes to flip in RMA normalised data:
obsfliprma = c('resobsgene')

matheat = gtmasstd[commongenes,]
matheat[matheat < -3] = -3
matheat[matheat > 3] = 3

pdf('pdf/gatza_pdf/gtobmas.pdf')
for (i in 1:length(tmpname)) {
	genes = get(tmpname[i])
	tmp = gtmasstd[genes,]
	meta = make_metagene(tmp, flip = F, raw = F)
	ord = order(meta)
	heatmap.2(matheat[,ord], trace='none',scale='none', col='bluered', ColSideColors = bluered(length(meta))[rank(meta)][ord], Colv=NA, main=tmpname[i], cexRow=1.0)
}
dev.off()

## Obesity metagenes to flip in RMA normalised data:
obsflipmas = c('rawobsgene', 'caobsgene', 'caresobsgene')

###############################################################################
## Apply Gatza metagene transformation matrix in other data sets and check if
## the clustering is similar/same as in the Gatza data:

datanames = c('crstdzzz','crisstdzzz','fmstdzzz','gtzzzstd')
maintxt= c('zzz-normalised Creighton data','zzz-normalised Cris data','zzz-normalised FM data','zzz-normalised Gatza data')

pdf(file='pdf/gatza_pdf/gatza_meta_trans.pdf', width=7, height=7)

tmpvar = gsub('zzz', 'rma', datanames)
tmptxt = gsub('zzz', 'RMA', maintxt)
for (i in 1:length(tmpvar)) {
	mat = get(tmpvar[i])
	flip = paths %in% mgfliprma
	tmpmeta = make_multi_transmeta(mat, paths, gtrmatransmat, flip = flip, raw = T)
	rownames(tmpmeta) = gsub('_probes', '', rownames(tmpmeta))
	tmpmeta = t(apply(tmpmeta, 1, function(x) convert_raw_meta(x, method = 'rank')))
	gatza_heat(tmpmeta, type = 'ranked', main = tmptxt[i])
}

tmpvar = gsub('zzz', 'mas', datanames)
tmptxt = gsub('zzz', 'MAS5', maintxt)
for (i in 1:length(tmpvar)) {
	mat = get(tmpvar[i])
	flip = paths %in% mgflipmas
	tmpmeta = make_multi_transmeta(mat, paths, gtmastransmat, flip = flip, raw = T)
	rownames(tmpmeta) = gsub('_probes', '', rownames(tmpmeta))
	tmpmeta = t(apply(tmpmeta, 1, function(x) convert_raw_meta(x, method = 'rank')))
	gatza_heat(tmpmeta, type = 'ranked', main = tmptxt[i])
}

dev.off()

## The heatmaps show similar clustering in other data sets as the clustering in the Gatza data set.

###############################################################################
## Apply all metagene transformation matrix (all derived from Gatza data) in
## other data sets and check which pathways cluster together with the obesity
## metagenes

# make all transformation matrices in Gatza data
allmeta = c(paths, obsname)
varname = allmeta
varname = gsub('_probes', '', varname)
varname = gsub('genes', '', varname)
varname = gsub('gene', '', varname)

allrmatransmat = list()
for (i in 1:length(allmeta)) {
	gene = get(allmeta[i])
	mat = gtrmastd[gene,]
	trans = make_trans_mat(mat)
	allrmatransmat[[i]] = trans
}
names(allrmatransmat) = varname

allmastransmat = list()
for (i in 1:length(allmeta)) {
	gene = get(allmeta[i])
	mat = gtmasstd[gene,]
	trans = make_trans_mat(mat)
	allmastransmat[[i]] = trans
}
names(allmastransmat) = varname

mgfliprmaall = c(mgfliprma, obsfliprma)
mgflipmasall = c(mgflipmas, obsflipmas)

pdf(file='pdf/gatza_pdf/all_meta_trans.pdf', width=7, height=7)

tmpvar = gsub('zzz', 'rma', datanames)
tmptxt = gsub('zzz', 'RMA', maintxt)
for (i in 1:length(tmpvar)) {
	mat = get(tmpvar[i])
	flip = allmeta %in% mgfliprmaall
	tmpmeta = make_multi_transmeta(mat, allmeta, allrmatransmat, flip = flip, raw = T)
	rownames(tmpmeta) = gsub('_probes', '', rownames(tmpmeta))
	rownames(tmpmeta) = gsub('genes', '', rownames(tmpmeta))
	rownames(tmpmeta) = gsub('gene', '', rownames(tmpmeta))
	tmpmeta = t(apply(tmpmeta, 1, function(x) convert_raw_meta(x, method = 'rank')))
	gatza_heat(tmpmeta, type = 'ranked', main = tmptxt[i])
}

## It seems like the obesity metagenes do not cluster together with any other pathways. FM obesity gene clustered together with BCAT/E2F1/... group, but didn't look like it had a significant correlation with them.

tmpvar = gsub('zzz', 'mas', datanames)
tmptxt = gsub('zzz', 'MAS5', maintxt)
for (i in 1:length(tmpvar)) {
	mat = get(tmpvar[i])
	flip = allmeta %in% mgflipmasall
	tmpmeta = make_multi_transmeta(mat, allmeta, allmastransmat, flip = flip, raw = T)
	rownames(tmpmeta) = gsub('_probes', '', rownames(tmpmeta))
	rownames(tmpmeta) = gsub('genes', '', rownames(tmpmeta))
	rownames(tmpmeta) = gsub('gene', '', rownames(tmpmeta))
	tmpmeta = t(apply(tmpmeta, 1, function(x) convert_raw_meta(x, method = 'rank')))
	gatza_heat(tmpmeta, type = 'ranked', main = tmptxt[i])
}

## In the MAS5-normalised data, it seems like the obesity metagenes from Creighton's data split up into two groups: the Creighton overlap metagenes, and non-overlap metagenes. The overlap groups formed a cluster on its own, but the non-overlap group clustered together with TGFB pathway strongly. FM metagene again clustered together with BCAT/E2F1 group (slightly stronger than in RMA normalised data).

dev.off()

###############################################################################
# Begin on the linear model stuff - predicting BMI metagenes using pathway
# metagenes (in Cris' data)

# First, drop some pathways that didn't look good over different data sets, using previous results
# There were six 'decent' pathways to use:
usepath = c("bcat_probes","er_probes","ifna_probes","ifng_probes","myc_probes","pr_probes")

# Make the data for linear model:

# metagene results from RMA (transmat and data) (Cris's data)
gtmeta = compare_svd_trans(crisstdrma, paths, gtrmatransmat)
gtmeta = gtmeta$trans[,usepath]
obsmeta = compare_svd_trans(crisstdrma, obsname, obsrmatransmat)
obsmeta = obsmeta$trans
crisdata = cbind(crisclin, gtmeta, obsmeta)

# Create linear model:
# Linear model with just BMI:
fit = lm(crobsgenes ~ bmi, crisdata)
summary(fit)
# Linear model with just BMI status:
fit = lm(crobsgenes ~ bmiStatus, crisdata)
summary(fit)
# Linear model with both BMI and BMI status:
fit = lm(crobsgenes ~ bmi + bmiStatus, crisdata)
summary(fit)

# Linear model with other pathway metagenes:
# Linear model with just BMI:
fit2 = lm(crobsgenes ~ bmi + bcat_probes + er_probes + ifna_probes + ifng_probes + myc_probes + pr_probes, crisdata)
summary(fit2)
# Linear model with just BMI status:
fit2 = lm(crobsgenes ~ bmiStatus + bcat_probes + er_probes + ifna_probes + ifng_probes + myc_probes + pr_probes, crisdata)
summary(fit2)
# Linear model with both BMI and BMI status:
fit2 = lm(crobsgenes ~ bmi + bmiStatus + bcat_probes + er_probes + ifna_probes + ifng_probes + myc_probes + pr_probes, crisdata)
summary(fit2)

# Linear model with pathway metagenes only:
fit3 = lm(crobsgenes ~ bcat_probes + er_probes + ifna_probes + ifng_probes + myc_probes + pr_probes, crisdata)
summary(fit3)

# Linear model with PR pathway metagene only:
fit4 = lm(crobsgenes ~ pr_probes, crisdata)
summary(fit4)

###############################################################################
# Predict obesity metagene in Cris' data, using the linear model that was fitted in Cris' data:

prediction1 = predict(fit, crisdata, se.fit = T)
cor(prediction1$fit, crisdata$crobsgenes, method = 'spearman')
cor(prediction1$fit, crisdata$crobsgenes, method = 'pearson')

prediction2 = predict(fit2, crisdata, se.fit = T)
cor(prediction2$fit, crisdata$crobsgenes, method = 'spearman')
cor(prediction2$fit, crisdata$crobsgenes, method = 'pearson')

prediction3 = predict(fit3, crisdata, se.fit = T)
cor(prediction3$fit, crisdata$crobsgenes, method = 'spearman')
cor(prediction3$fit, crisdata$crobsgenes, method = 'pearson')

prediction4 = predict(fit4, crisdata, se.fit = T)
cor(prediction4$fit, crisdata$crobsgenes, method = 'spearman')
cor(prediction4$fit, crisdata$crobsgenes, method = 'pearson')

pdf(file='pdf/gatza_pdf/prediction_cris_with_cris.pdf', width=7, height=7)
plot(crisdata$crobsgenes, prediction1$fit, pch=20)
plot(crisdata$crobsgenes, prediction2$fit, pch=20)
plot(crisdata$crobsgenes, prediction3$fit, pch=20)
plot(crisdata$crobsgenes, prediction4$fit, pch=20)
dev.off()

###############################################################################
# Predict obesity metagene in Creighton's data, using the linear model that was fitted in Cris' data:

# Make the data:

# metagene results from RMA (transmat and data) (Creighton's data)
gtmeta = compare_svd_trans(crstdrma, paths, gtrmatransmat)
gtmeta = gtmeta$trans[,usepath]
obsmeta = compare_svd_trans(crstdrma, obsname, obsrmatransmat)
obsmeta = obsmeta$trans
crdata = cbind(crclin, gtmeta, obsmeta)

# Predict the BMI metagene based on the value of the pathway metagenes:
prediction1 = predict(fit, crdata, se.fit = T)
cor(prediction1$fit, crdata$crobsgenes, method = 'spearman')
cor(prediction1$fit, crdata$crobsgenes, method = 'pearson')

prediction2 = predict(fit2, crdata, se.fit = T)
cor(prediction2$fit, crdata$crobsgenes, method = 'spearman')
cor(prediction2$fit, crdata$crobsgenes, method = 'pearson')

prediction3 = predict(fit3, crdata, se.fit = T)
cor(prediction3$fit, crdata$crobsgenes, method = 'spearman')
cor(prediction3$fit, crdata$crobsgenes, method = 'pearson')

prediction4 = predict(fit4, crdata, se.fit = T)
cor(prediction4$fit, crdata$crobsgenes, method = 'spearman')
cor(prediction4$fit, crdata$crobsgenes, method = 'pearson')

pdf(file='pdf/gatza_pdf/prediction_cr_with_cris.pdf', width=7, height=7)
plot(crdata$crobsgenes, prediction1$fit, pch=20, ylab='Predicted metagene', xlab='True metagene', main='BMI-only')
plot(crdata$crobsgenes, prediction2$fit, pch=20, ylab='Predicted metagene', xlab='True metagene', main='BMI and pathways')
plot(crdata$crobsgenes, prediction3$fit, pch=20, ylab='Predicted metagene', xlab='True metagene', main='Pathways-only')
plot(crdata$crobsgenes, prediction4$fit, pch=20, ylab='Predicted metagene', xlab='True metagene', main='PR-only')
dev.off()

###############################################################################
# Use all the pathways, but drop them out one by one and find the best set of pathways that predict the obesity metagenes








###############################################################################
## TODO: work from here
print('hello world!') ## command to stop nvimr to go down to the bottom


