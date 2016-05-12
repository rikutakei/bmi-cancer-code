###############################################################################

                        ##R code for task 13

###############################################################################

## Combine tasks 1 to 12 into this document

###############################################################################
## pre-analysis stuff
setwd('~/Documents/masters/data/task13/')

#source the file to get all the functions from task10.R
source('~/Documents/codes/bmi-cancer-code/task13/functions.R')

# install packages
source('http://bioconductor.org/biocLite.R')
install.packages("data.table")
install.packages("gplots")
biocLite()
a
y
biocLite("WGCNA")
a
y
biocLite("hgu133a.db")
a
y
biocLite("tidyr")
a
y
biocLite("dplyr")
a
y
biocLite("limma")
a
y
biocLite("edgeR")
a
y
biocLite("DESeq")
a
y
biocLite("affy")
a
y
biocLite("org.Hs.eg.db")
a
y
biocLite("KEGG.db")
a
y
biocLite("GO.db")
a
y
biocLite("reactome.db")
a
y

# load libraries:

library(data.table)
library(gplots)
library(lattice)

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

###############################################################################

## First load all of the necessary data into R:

###############################################################################
## Creighton et al data:

files = readLines('./raw/creighton/files.txt')
files = paste('./raw/creighton/', files, sep='')
cr_raw = ReadAffy(filenames = files)
cr_raw = rma(cr_raw) ## RMA normalise the data
cr_raw = exprs(cr_raw) ## change the format into matrix

cr_obsgene = read.csv('./obsgenes/crobsgenes.txt', header=F)
cr_obsgene = as.vector(cr_obsgene[,1])

crclin = read.csv('./clindata/crclin.csv', sep=',', header=T)

###############################################################################
## Fuentes-Mattei et al data:

files = readLines('./raw/fuentes-mattei/files.txt')
files = paste('./raw/fuentes-mattei/', files, sep='')
fm_raw = ReadAffy(filenames = files)
fm_raw = rma(fm_raw) ## RMA normalise the data
fm_raw = exprs(fm_raw) ## change the format into matrix

fm_obsgene = read.csv('./obsgenes/fmobsgenes.txt', header=T)

fmclin = read.csv('./clindata/fmclin.csv', sep=',', header=T)

###############################################################################
## ICGC data:

## (may have to use the server for extra RAM (ulimit -s 20480))

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

	## process bmi data:
	bmi = clin[,c('height', 'weight')]
	rownames(bmi) = clin$bcr_patient_barcode
	if (length(which(bmi[,1] == "[Not Available]")) > 0){
		bmi = bmi[-which(bmi[,1] == "[Not Available]"),]
	}
	if (length(which(bmi[,2] == "[Not Available]")) > 0){
		bmi = bmi[-which(bmi[,2] == "[Not Available]"),]
	}

	seq = seq[which(rownames(seq) %in% rownames(bmi)),]
	bmi = bmi[which(rownames(bmi) %in% rownames(seq)),]

	## assign variable names:
	assign(files[i], seq)
	assign(paste(files[i], 'clin', sep=''), clin)
	assign(paste(files[i], 'bmi', sep=''), bmi)
}

###############################################################################
## Pathway data base:

#Import Human Gene Symbols
SYMBOL.list<-as.list(org.Hs.egSYMBOL)

KEGG.list<-as.list(org.Hs.egPATH) ##Import KEGG pathways:
names(KEGG.list)<-unlist(SYMBOL.list) #name KEGG pathway lists with corresponding gene symbols
keggpath<-as.list(KEGGPATHID2NAME) #mapping KEGG path IDs to human read pathway name

GO.list<-as.list(org.Hs.egGO) ##Import GO pathways:
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
names(GO.list)<-SYMBOL.list[names(GO.list)] #name GO pathway lists with corresponding gene symbols
gopath<-as.list(GOTERM) #mapping GO IDs to human read pathway name
#... and reformat so matching KEGG.list
tmp = list()
for(i in 1:length(gopath)) tmp[[i]]<-gopath[[i]]@Term
names(tmp) = names(gopath)
gopath = tmp

#Import Human Reactome pathways
reactome.list <-as.list(reactomeEXTID2PATHID)
reactome.list = reactome.list[order(as.numeric(names(reactome.list)))] #Sort the names in the list:
reactome.list = reactome.list[which(names(reactome.list) %in% names(SYMBOL.list))] #get paths that have gene symbols:
names(reactome.list) = unlist(SYMBOL.list[names(reactome.list)]) #rename the entrez gene ID into gene symbol
reactomepath <- as.list(reactomePATHID2NAME) #mapping readtome path IDs to human read pathway name
reactomepath = reactomepath[grep('Homo sapiens',reactomepath)] #pull out all human-related pathways
reactomepath = lapply(reactomepath, function(x) gsub('Homo sapiens: ', '', x)) #cut out the 'Homo sapiens: ' bit so it's only the pathway names.

###############################################################################
## Creighton metagene stuff

## make matrix with just the obesity-associated genes:
cr_obsmat = cr_raw[cr_obsgene,]

## map the gene probe names back to the gene symbols:
cr_genesymbol = mapIds(hgu133a.db, keys = rownames(cr_obsmat), column = 'SYMBOL', keytype = "PROBEID", multiVals = 'first')
rownames(cr_obsmat) = as.vector(cr_genesymbol)

## remove NAs and duplicated genes
cr_obsmat = cr_obsmat[-(which(is.na(rownames(cr_obsmat)))),]
cr_obsmat = collapseRows(cr_obsmat, rowGroup = unique(rownames(cr_obsmat)), rowID = unique(rownames(cr_obsmat)))
cr_obsmat = cr_obsmat$datETcollapsed

## get genes that are common from the creighton gene sets and the ICGC data set and use these genes to make the matrix
genenames = colnames(BLCA)
common_genes = which(rownames(cr_obsmat) %in% genenames)
cr_obsmat = cr_obsmat[common_genes,]
dim(cr_obsmat) ## 644 genes

## TODO: Do the heatmap stuff and add p-values to the plots
cr_rawsvd = svd(cr_obsmat)
cr_rawmeta  = cr_rawsvd$v[,1] # first principle component
cr_rawmeta2 = cr_rawsvd$v[,2] # second principle component
cr_rawmeta3 = cr_rawsvd$v[,3] # third principle component
##cr_rawmeta = rank(cr_rawmeta$v[,1])/ncol(cr_obsmat)
##cr_rawmeta = 1-cr_rawmeta
cr_raword = order(cr_rawmeta)

cr_obsmat_adj = t(apply(cr_obsmat, 1, function(x) (x-mean(x))/sd(x)))
cr_obsmat_adj[cr_obsmat_adj > 3] = 3
cr_obsmat_adj[cr_obsmat_adj < -3] = -3

cr_adjsvd = svd(cr_obsmat_adj)
cr_adjmeta  = cr_adjsvd$v[,1] # first principle component
cr_adjmeta2 = cr_adjsvd$v[,2] # second principle component
cr_adjmeta3 = cr_adjsvd$v[,3] # third principle component
cr_adjmeta = rank(cr_adjmeta)/ncol(cr_obsmat)
cr_adjmeta = 1-cr_adjmeta
cr_adjord = order(cr_adjmeta)

# see if it the raw or adjusted metagenes have any difference
plot(1-cr_adjmeta, -cr_rawmeta, pch=20, main='Creighton metagene comparison', ylab='Raw metagene value', xlab='Adjusted metagene value')

## Check if the metagene correlates with BMI:
heatmap.2(cr_obsmat_adj, trace='none',scale='none', col='bluered')
heatmap.2(cr_obsmat_adj, trace='none',scale='none', col='bluered', ColSideColors = bluered(length(cr_rawmeta))[rank(cr_rawmeta)])
heatmap.2(cr_obsmat_adj[,cr_raword], trace='none',scale='none', col='bluered', ColSideColors = bluered(length(cr_rawmeta))[rank(cr_rawmeta)][cr_raword], Colv=F)
bmifactor = factor(crclin$bmiStatus, levels=c("normal","overweight", "obese"))
boxplot(cr_rawmeta~bmifactor)
plot(crclin$bmi, cr_rawmeta, pch = 20)
fit = lm(cr_rawmeta~crclin$bmi)
abline(fit)
txt = summary(fit)$adj.r.squared ## r squared = 0.145
txt2 = summary(fit)$coef[2, 4] ## p value  = 4.30e-05
legend('bottomright', bty='n', legend = c(as.expression(bquote(R^2 == .(format(txt, digits = 4)))), as.expression(bquote(p == .(format(txt2, digits=4))))))

# repeat with adjusted data:
heatmap.2(cr_obsmat_adj, trace='none',scale='none', col='bluered')
heatmap.2(cr_obsmat_adj, trace='none',scale='none', col='bluered', ColSideColors = bluered(length(cr_adjmeta))[rank(cr_adjmeta)])
heatmap.2(cr_obsmat_adj[,cr_adjord], trace='none',scale='none', col='bluered', ColSideColors = bluered(length(cr_adjmeta))[rank(cr_adjmeta)][cr_adjord], Colv=F)
boxplot(cr_adjmeta~bmifactor)
plot(crclin$bmi, cr_adjmeta, pch = 20)
fit = lm(cr_adjmeta~crclin$bmi)
abline(fit, col="red")
txt = summary(fit)$adj.r.squared ## r squared = 0.105
txt2 = summary(fit)$coef[2, 4] ## p value  = 4.90e-04
legend('bottomright', bty='n', legend = c(as.expression(bquote(R^2 == .(format(txt, digits = 4)))), as.expression(bquote(p == .(format(txt2, digits=4))))))
dend = hclust(dist(cr_obsmat_adj))

## Make transformation matrix:
cr_rawtransmat = diag(1/cr_rawsvd$d) %*% t(cr_rawsvd$u)
cr_adjtransmat = diag(1/cr_adjsvd$d) %*% t(cr_adjsvd$u)

## Pull out all the genes common in Creighton obesity genes and ICGC data:
cancertypes = c('BLCA', 'CESC', 'COAD', 'KIRP', 'LIHC', 'READ', 'SKCM', 'UCEC')
for (i in 1:length(cancertypes)) {
	seq = get(cancertypes[i])
	seq = seq[,common_genes]
	seq = t(seq) ## make it genes by samples for voomICGC function
	txt = paste('cr', cancertypes[i], sep='')
	assign(txt, seq)
}

## check if it has the right number of genes:
dim(crBLCA) # 644 genes

## calculate bmi and assign bmistatus
for (i in 1:length(cancertypes)) {
	txt = paste(cancertypes[i], 'bmi', sep='')
	tmp = get(txt)
	tmp[,1] = as.numeric(as.character(tmp[,1]))
	tmp[,2] = as.numeric(as.character(tmp[,2]))
	x = tmp$weight/((tmp$height/100)^2)
	tmp[,3] = x
	for (j in 1:length(x)) {
		if(x[j] >= 30){
			x[j] = "obese"
		} else if(x[j] < 25){
			x[j] = "normal"
		} else {
			x[j] = "overweight"
		}
	}
	tmp[,4] = x
	colnames(tmp)[3] = "bmi"
	colnames(tmp)[4] = "bmiStatus"
	assign(txt, tmp)
}

## log10 and standardise data:
for (i in 1:length(cancertypes)) {
	txt = paste('cr', cancertypes[i], sep='')
	seq = get(txt)
	seq = standardise_data(seq)
	assign(txt, seq)
}

## Transform the ICGC data:
for (i in 1:length(cancertypes)) {
	t = t(cr_adjtransmat %*% get(paste('cr', cancertypes[i], sep='')))
	t = t[,1]
	assign(paste(cancertypes[i], 'cradjmeta',sep=''), t)
}

## Check if the metagene correlates with sample gene expression and/or BMI:
for (i in 1:length(cancertypes)) {
	dat = get(paste('cr', cancertypes[i], sep=''))
	meta = get(paste(cancertypes[i], 'cradjmeta', sep=''))
	bmi = get(paste(cancertypes[i], 'bmi', sep=''))
	metaplot(dat, meta, bmi, dend, name = paste(cancertypes[i], '- Creighton metagene'))
}

###############################################################################
## Fuentes-Mattei metagene stuff

## make matrix with just the obesity-associated genes:
fm_obsmat = fm_raw[fm_obsgene$Probe,]

## map the gene probe names back to the gene symbols:
fm_genesymbol = mapIds(hgu133a.db, keys = rownames(fm_obsmat), column = 'SYMBOL', keytype = "PROBEID", multiVals = 'first')
rownames(fm_obsmat) = as.vector(fm_genesymbol)

## remove duplicated genes
fm_obsmat = collapseRows(fm_obsmat, rowGroup = unique(rownames(fm_obsmat)), rowID = unique(rownames(fm_obsmat)))
fm_obsmat = fm_obsmat$datETcollapsed

## get genes that are common from the creighton gene sets and the ICGC data set and use these genes to make the matrix
genenames = colnames(BLCA)
crgenenames = rownames(cr_raw)
crgenenames = mapIds(hgu133a.db, keys = crgenenames, column = 'SYMBOL', keytype = "PROBEID", multiVals = 'first')
crgenenames = unique(crgenenames)
common_genes = which(rownames(fm_obsmat) %in% genenames)
common_genes = rownames(fm_obsmat)[common_genes]
common_genes = common_genes[which(common_genes %in% crgenenames)]
fm_obsmat = fm_obsmat[common_genes,]
dim(fm_obsmat) ## 116 genes

## TODO: Do the heatmap stuff and add p-values to the plots
fm_rawsvd = svd(fm_obsmat)
fm_rawmeta  = fm_rawsvd$v[,1] # first principle component
fm_raword = order(fm_rawmeta)

fm_obsmat_adj = t(apply(fm_obsmat, 1, function(x) (x-mean(x))/sd(x)))
fm_obsmat_adj[fm_obsmat_adj > 3] = 3
fm_obsmat_adj[fm_obsmat_adj < -3] = -3

fm_adjsvd = svd(fm_obsmat_adj)
fm_adjmeta  = fm_adjsvd$v[,1] # first principle component
fm_adjmeta = rank(fm_adjmeta)/ncol(fm_obsmat)
fm_adjmeta = 1-fm_adjmeta
fm_adjord = order(fm_adjmeta)

# see if it the raw or adjusted metagenes have any difference
plot(1-fm_adjmeta, -fm_rawmeta, pch=20, main='FM metagene comparison', ylab='Raw metagene value', xlab='Adjusted metagene value')

heatmap.2(fm_obsmat_adj, trace='none',scale='none', col='bluered')
heatmap.2(fm_obsmat_adj, trace='none',scale='none', col='bluered', ColSideColors = bluered(length(-fm_rawmeta))[rank(fm_rawmeta)])
heatmap.2(fm_obsmat_adj[,fm_raword], trace='none',scale='none', col='bluered', ColSideColors = bluered(length(fm_rawmeta))[rank(-fm_rawmeta)][fm_raword], Colv=F)

# repeat with adjusted metagene:
heatmap.2(fm_obsmat_adj, trace='none',scale='none', col='bluered')
heatmap.2(fm_obsmat_adj, trace='none',scale='none', col='bluered', ColSideColors = bluered(length(fm_adjmeta))[rank(fm_adjmeta)])
heatmap.2(fm_obsmat_adj[,fm_adjord], trace='none',scale='none', col='bluered', ColSideColors = bluered(length(fm_adjmeta))[rank(fm_adjmeta)][fm_adjord], Colv=F)

## Make transformation matrix:
fm_rawtransmat = diag(1/fm_rawsvd$d) %*% t(fm_rawsvd$u)
fm_adjtransmat = diag(1/fm_adjsvd$d) %*% t(fm_adjsvd$u)

###############################################################################
## Get metagene using the FM probe set, rather than transforming it with
## transformation matrix

cr_fmobsmat = cr_raw[fm_obsgene$Probe,]

tmp = mapIds(hgu133a.db, keys = rownames(cr_fmobsmat), column = 'SYMBOL', keytype = "PROBEID", multiVals = 'first')
rownames(cr_fmobsmat) = as.vector(tmp)
cr_fmobsmat = cr_fmobsmat[common_genes,]
dim(cr_fmobsmat) ## 116 genes

cr_fmrawsvd = svd(cr_fmobsmat)
cr_fmrawmeta  = cr_fmrawsvd$v[,1] # first principle component
cr_fmraword = order(cr_fmrawmeta)

cr_fmobsmat_adj = t(apply(cr_fmobsmat, 1, function(x) (x-mean(x))/sd(x)))
cr_fmobsmat_adj[cr_fmobsmat_adj > 3] = 3
cr_fmobsmat_adj[cr_fmobsmat_adj < -3] = -3

cr_fmadjsvd = svd(cr_fmobsmat_adj)
cr_fmadjmeta  = cr_fmadjsvd$v[,1] # first principle component
cr_fmadjmeta = rank(cr_fmadjmeta)/ncol(cr_fmobsmat)
cr_fmadjmeta = 1-cr_fmadjmeta
cr_fmadjord = order(cr_fmadjmeta)

# see if it the raw or adjusted metagenes have any difference
plot(1-cr_fmadjmeta, -cr_fmrawmeta, pch=20, main='FM metagene comparison', ylab='Raw metagene value', xlab='Adjusted metagene value')
cr_fmdend= hclust(dist(cr_fmobsmat_adj))

metaplot(cr_fmobsmat_adj, cr_fmadjmeta, crclin, cr_fmdend)

###############################################################################
## Get metagene using the transformation matrix

cr_fmobsmat = cr_raw[fm_obsgene$Probe,]

tmp = mapIds(hgu133a.db, keys = rownames(cr_fmobsmat), column = 'SYMBOL', keytype = "PROBEID", multiVals = 'first')
rownames(cr_fmobsmat) = as.vector(tmp)
cr_fmobsmat = cr_fmobsmat[common_genes,]
dim(cr_fmobsmat) ## 116 genes

cr_fmtransmeta = t(fm_rawtransmat %*% cr_fmobsmat)
cr_fmtransmeta = cr_fmtransmeta[,1]
cr_fmtransmeta = rank(cr_fmtransmeta)/ncol(cr_fmobsmat)
cr_fmtransmeta = 1-cr_fmtransmeta

metaplot(cr_fmobsmat_adj, cr_fmtransmeta, crclin, cr_fmdend)

###############################################################################
## Try it on ICGC data

## Pull out all the genes common in Creighton obesity genes and ICGC data:
cancertypes = c('BLCA', 'CESC', 'COAD', 'KIRP', 'LIHC', 'READ', 'SKCM', 'UCEC')
for (i in 1:length(cancertypes)) {
	seq = get(cancertypes[i])
	seq = seq[,common_genes]
	seq = t(seq) ## make it genes by samples for voomICGC function
	txt = paste('fm', cancertypes[i], sep='')
	assign(txt, seq)
}

## log10 and standardise data:
for (i in 1:length(cancertypes)) {
	txt = paste('fm', cancertypes[i], sep='')
	seq = get(txt)
	seq = standardise_data(seq)
	assign(txt, seq)
}

## Transform the ICGC data:
for (i in 1:length(cancertypes)) {
	t = t(fm_adjtransmat %*% get(paste('fm', cancertypes[i], sep='')))
	t = t[,1]
	assign(paste(cancertypes[i], 'fmadjmeta',sep=''), t)
}

## Check if the metagene correlates with sample gene expression and/or BMI:
for (i in 1:length(cancertypes)) {
	dat = get(paste('fm', cancertypes[i], sep=''))
	meta = get(paste(cancertypes[i], 'fmadjmeta', sep=''))
	bmi = get(paste(cancertypes[i], 'bmi', sep=''))
	metaplot2(dat, meta, bmi, name = paste(cancertypes[i], '- FM metagene'))
}

###############################################################################
## Creighton gene expression analysis stuff

## Make group:
group = crclin$bmiStatus
group = ifelse(group=='obese', 'obese', 'normal/overweight')
design = model.matrix(~group)
fit = lmFit(cr_raw, design)
fit = eBayes(fit)
tt = topTable(fit, coef = 2, adjust = 'BH', n=nrow(cr_raw))
sum(tt$P.Value < 0.05)   # 5278 genes
sum(tt$P.Value < 0.01)   # 1781 genes
sum(tt$adj.P.Val < 0.05) # 9 genes
sum(tt$adj.P.Val < 0.01) # 0 genes
length(which(2^tt$logFC > 1.2)) #only 132 genes

## How many Creighton genes are p < 0.01
rawobsgenes = tt[tt$P.Value < 0.01,]
length(which( cr_obsgene %in% rownames(rawobsgenes) )) ## 390
## What about in top 799 genes?
rawobsgenes = rownames(rawobsgenes[1:799,])
crolgenes = rawobsgenes[which(rawobsgenes %in% cr_obsgene)]
length(crolgenes) ## 239

## Do the same, but in residual data:

## Fit linear model on the data using all the other clinical variables:
residuals = cr_raw
residuals = t(apply(residuals, 1, function(x) lm(x ~ crclin$age + crclin$race + crclin$menopause + crclin$grade + crclin$LNstatus + crclin$ERstatus + crclin$PRstatus + crclin$HER2status)$residuals))

## Fit linear model on the residuals, using obese vs non-obese design matrix
resfit = lmFit(residuals, design)
resfit = eBayes(resfit)
restt = topTable(resfit, coef=2, n = nrow(residuals))
sum(restt$P.Value < 0.01)   # 1104 genes
sum(restt$adj.P.Val < 0.01) # 0 genes

## How many Creighton genes are p < 0.01
resobsgene = restt[restt$P.Value < 0.01,]
length(which( cr_obsgene %in% rownames(resobsgene ) )) ## 216
## What about in top 799 genes?
resobsgene = rownames(resobsgene[1:799,])
rescrolgenes = resobsgene[which(resobsgene %in% cr_obsgene)]
length(rescrolgenes) ## 168

## make venn diagram with res, raw, and cr obesity genes:
setlist = list(Residual = resobsgene, My_genes = rawobsgenes, Creighton = cr_obsgene)
OLlist = overLapper(setlist=setlist, sep='_', type='vennsets')
counts = sapply(OLlist$Venn_List, length)
vennPlot(counts = counts)

## Do the same, but in Caucasian-only data:

camat =cr_raw
ind = which(crclin$race == 'w')
camat = camat[,ind] ## 77 caucasian samples in total
caclin = crclin[ind,]

group = caclin$bmiStatus
group = ifelse(group == 'obese', 'obese', 'normal/overweight')
design = model.matrix(~group)
cafit = lmFit(camat, design)
cafit = eBayes(cafit)
catt = topTable(cafit, coef = 2, n = nrow(camat))
sum(catt$P.Value < 0.01)   # 2129 genes
sum(catt$adj.P.Val < 0.01) # 0 genes

## How many Creighton genes are p < 0.01
caobsgene = catt[catt$P.Value < 0.01,]
length(which( cr_obsgene %in% rownames(caobsgene ) )) ## 271
## What about in top 799 genes?
caobsgene = rownames(caobsgene[1:799,])
cacrolgenes = caobsgene[which(caobsgene %in% cr_obsgene)]
length(cacrolgenes) ## 148

## Fit linear model on the Caucasian-only data  after removing all the other clinical variables:
caresiduals = camat
caresiduals = t(apply(caresiduals, 1, function(x) lm(x ~ caclin$age + caclin$menopause + caclin$grade + caclin$LNstatus + caclin$ERstatus + caclin$PRstatus + caclin$HER2status)$residuals))

## Fit linear model on the residuals, using obese vs non-obese design matrix
caresfit = lmFit(caresiduals, design)
caresfit = eBayes(caresfit)
carestt = topTable(caresfit, coef=2, n = nrow(caresiduals))
sum(carestt$P.Value < 0.01)   # 1558 genes
sum(carestt$adj.P.Val < 0.01) # 0 genes

## How many Creighton genes are p < 0.01
caresobsgene = carestt[carestt$P.Value < 0.01,]
length(which( cr_obsgene %in% rownames(caresobsgene ) )) ## 166
## What about in top 799 genes?
caresobsgene = rownames(caresobsgene[1:799,])
carescrolgenes = caresobsgene[which(caresobsgene %in% cr_obsgene)]
length(carescrolgenes) ## 92

## make a Venn diagram with ca, cares, and cr obesity genes:
setlist = list(Caucasian_residual = caresobsgene , Caucasian= caobsgene, Creighton = cr_obsgene)
OLlist = overLapper(setlist=setlist, sep='_', type='vennsets')
counts = sapply(OLlist$Venn_List, length)
vennPlot(counts = counts)

## validate these metagenes in creighton data:
allobsname = c("rawobsgenes","crolgenes","resobsgene","rescrolgenes", "caobsgene","cacrolgenes","caresobsgene","carescrolgenes")
metalist = list()
for (i in 1:length(allobsname)) {
	gene = get(allobsname[i])
	mat = cr_raw[gene,]
	mat = t(apply(mat, 1, function(x) (x-mean(x))/sd(x)))
	mat[mat < -3] = -3
	mat[mat > 3] = 3
	tmpsvd = svd(mat)
	tmpsvd = rank(tmpsvd$v[,1])/ncol(mat)
	metaplot2(mat, tmpsvd, crclin, name = allobsname[i])
	metalist[[i]] = tmpsvd
}

## make correlation matrix
x = matrix(1,103)
for(i in 1:length(metalist)) {
	x = cbind(x, metalist[[i]])
}
x = x[,-1]
levelplot(cor(x, method = 'pearson'))

## check if this works on ICGC data:





###############################################################################
## Common genes across ICGC data stuff

## use codes from task11?



###############################################################################
## Pathway enrichment stuff

## TODO: do the enrichment analysis with reactome and kegg













