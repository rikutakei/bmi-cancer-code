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
install.packages("RColorBrewer")
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
biocLite("GMD")
a
y

# load libraries:

library(data.table)
library(gplots)
library(lattice)
library(GMD)
library(RColorBrewer)

#heatmap stuff
library("devtools")
devtools::install_github("TomKellyGenetics/heatmap.2x", ref="master")
library("heatmap.2x")

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
cr_rawmeta = -cr_rawmeta
cr_rawmeta = 1-cr_rawmeta
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

pdf('crmeta1.pdf')
# see if it the raw or adjusted metagenes have any difference
plot(cr_adjmeta, cr_rawmeta, pch=20, main='Creighton metagene comparison', ylab='Raw metagene value', xlab='Adjusted metagene value')

main='Creighton metagene (raw)'
## Check if the metagene correlates with BMI:
heatmap.2(cr_obsmat_adj, trace='none',scale='none', col='bluered', main=main)
heatmap.2(cr_obsmat_adj, trace='none',scale='none', col='bluered', ColSideColors = bluered(length(cr_rawmeta))[rank(cr_rawmeta)], main=main)
heatmap.2(cr_obsmat_adj[,cr_raword], trace='none',scale='none', col='bluered', ColSideColors = bluered(length(cr_rawmeta))[rank(cr_rawmeta)][cr_raword], Colv=F, main=main)
bmifactor = factor(crclin$bmiStatus, levels=c("normal","overweight", "obese"))
boxplot(cr_rawmeta~bmifactor, ylab = "Raw metagene", xlab = "BMI Status", main=paste(main, ' vs. BMI Status', sep=''))

normind = which('normal' == bmifactor)
ovind = which('overweight' == bmifactor)
obind = which('obese' == bmifactor)
txt = t.test(cr_rawmeta[c(normalind,ovind)]~bmifactor[c(normalind,ovind)], alternative= 'two.sided')$p.value ## p = 0.3450
txt2 = t.test(cr_rawmeta[c(normalind,obind)]~bmifactor[c(normalind,obind)], alternative= 'two.sided')$p.value ## p = 6.061e-05
txt3 = summary(aov(cr_rawmeta~bmifactor))[[1]]$Pr[1]
legend(x = 1.6, y = 0.91, bty='n', legend = as.expression(bquote(p == .(format(txt, digits=4)))))
legend(x = 2.6, y = 0.91, bty='n', legend = as.expression(bquote(p == .(format(txt2, digits=4)))))
legend('bottomright', bty='n', legend = as.expression(bquote("ANOVA p" == .(format(txt3, digits=4)))))

plot(crclin$bmi, cr_rawmeta, pch = 20, xlab = "Raw metagene", ylab = "BMI", main=paste(main, ' vs. BMI', sep=''))
fit = lm(cr_rawmeta~crclin$bmi)
abline(fit, col='red')
txt = summary(fit)$adj.r.squared ## r squared = 0.145
txt2 = summary(fit)$coef[2, 4] ## p value  = 4.30e-05
legend('bottomright', bty='n', legend = c(as.expression(bquote(R^2 == .(format(txt, digits = 4)))), as.expression(bquote(p == .(format(txt2, digits=4))))))

main='Creighton metagene (adjusted)'
# repeat with adjusted data:
heatmap.2(cr_obsmat_adj, trace='none',scale='none', col='bluered', main=main)
heatmap.2(cr_obsmat_adj, trace='none',scale='none', col='bluered', ColSideColors = bluered(length(cr_adjmeta))[rank(cr_adjmeta)], main=main)
heatmap.2(cr_obsmat_adj[,cr_adjord], trace='none',scale='none', col='bluered', ColSideColors = bluered(length(cr_adjmeta))[rank(cr_adjmeta)][cr_adjord], Colv=F, main=main)
boxplot(cr_adjmeta~bmifactor, xlab = "Adjusted metagene", ylab = "BMI Status", main=paste(main, ' vs. BMI Status', sep=''), ylim = c(-0.05, 1.1))

normind = which('normal' == bmifactor)
ovind = which('overweight' == bmifactor)
obind = which('obese' == bmifactor)
txt = t.test(cr_adjmeta[c(normalind,ovind)]~bmifactor[c(normalind,ovind)], alternative= 'two.sided')$p.value ## p = 0.3450
txt2 = t.test(cr_adjmeta[c(normalind,obind)]~bmifactor[c(normalind,obind)], alternative= 'two.sided')$p.value ## p = 6.061e-05
txt3 = summary(aov(cr_adjmeta~bmifactor))[[1]]$Pr[1]
legend(x = 1.6, y = 1.08, bty='n', legend = as.expression(bquote(p == .(format(txt, digits=4)))))
legend(x = 2.6, y = 1.08, bty='n', legend = as.expression(bquote(p == .(format(txt2, digits=4)))))
legend('bottomright', bty='n', legend = as.expression(bquote("ANOVA p" == .(format(txt3, digits=4)))))

plot(crclin$bmi, cr_adjmeta, pch = 20, xlab = "Adjusted metagene", ylab = "BMI", main=paste(main, ' vs. BMI', sep=''))
fit = lm(cr_adjmeta~crclin$bmi)
abline(fit, col="red")
txt = summary(fit)$adj.r.squared ## r squared = 0.105
txt2 = summary(fit)$coef[2, 4] ## p value  = 4.90e-04
legend('bottomright', bty='n', legend = c(as.expression(bquote(R^2 == .(format(txt, digits = 4)))), as.expression(bquote(p == .(format(txt2, digits=4))))))
dend = hclust(dist(cr_obsmat_adj))

dev.off()

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
	t = rank(t)/length(t)
	t = 1-t
	assign(paste(cancertypes[i], 'cradjmeta',sep=''), t)
}

pdf('crtcga.pdf')

## Check if the metagene correlates with sample gene expression and/or BMI:
for (i in 1:length(cancertypes)) {
	dat = get(paste('cr', cancertypes[i], sep=''))
	meta = get(paste(cancertypes[i], 'cradjmeta', sep=''))
	bmi = get(paste(cancertypes[i], 'bmi', sep=''))
	txt = paste('Creighton metagene (', cancertypes[i], sep = "")
	txt = paste(txt, ')', sep = "")
	metaplot3(dat, meta, bmi, name = txt)
}

dev.off()

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

fm_rawsvd = svd(fm_obsmat)
fm_rawmeta  = fm_rawsvd$v[,1] # first principle component
fm_rawmeta = rank(fm_rawmeta$v[,1])/ncol(fm_obsmat)
fm_rawmeta = -fm_rawmeta
fm_rawmeta = 1-fm_rawmeta
fm_raword = order(fm_rawmeta)

fm_obsmat_adj = t(apply(fm_obsmat, 1, function(x) (x-mean(x))/sd(x)))
fm_obsmat_adj[fm_obsmat_adj > 3] = 3
fm_obsmat_adj[fm_obsmat_adj < -3] = -3

fm_adjsvd = svd(fm_obsmat_adj)
fm_adjmeta  = fm_adjsvd$v[,1] # first principle component
fm_adjmeta = rank(fm_adjmeta)/ncol(fm_obsmat)
fm_adjmeta = 1-fm_adjmeta
fm_adjord = order(fm_adjmeta)

pdf('fmmeta1.pdf')

# see if it the raw or adjusted metagenes have any difference
plot(1-fm_adjmeta, -fm_rawmeta, pch=20, main='FM metagene comparison', ylab='Raw metagene value', xlab='Adjusted metagene value')

main='FM metagene (raw)'
heatmap.2(fm_obsmat_adj, trace='none',scale='none', col='bluered', main=main)
heatmap.2(fm_obsmat_adj, trace='none',scale='none', col='bluered', ColSideColors = bluered(length(-fm_rawmeta))[rank(fm_rawmeta)], main=main)
heatmap.2(fm_obsmat_adj[,fm_raword], trace='none',scale='none', col='bluered', ColSideColors = bluered(length(fm_rawmeta))[rank(-fm_rawmeta)][fm_raword], Colv=F, main=main)

main='FM metagene (Adjusted)'
# repeat with adjusted metagene:
heatmap.2(fm_obsmat_adj, trace='none',scale='none', col='bluered', main=main)
heatmap.2(fm_obsmat_adj, trace='none',scale='none', col='bluered', ColSideColors = bluered(length(fm_adjmeta))[rank(fm_adjmeta)], main=main)
heatmap.2(fm_obsmat_adj[,fm_adjord], trace='none',scale='none', col='bluered', ColSideColors = bluered(length(fm_adjmeta))[rank(fm_adjmeta)][fm_adjord], Colv=F, main=main)

dev.off()

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
cr_fmrawmeta = -cr_fmrawmeta
cr_fmrawmeta = 1-cr_fmrawmeta
cr_fmraword = order(cr_fmrawmeta)

cr_fmobsmat_adj = t(apply(cr_fmobsmat, 1, function(x) (x-mean(x))/sd(x)))
cr_fmobsmat_adj[cr_fmobsmat_adj > 3] = 3
cr_fmobsmat_adj[cr_fmobsmat_adj < -3] = -3

cr_fmadjsvd = svd(cr_fmobsmat_adj)
cr_fmadjmeta  = cr_fmadjsvd$v[,1] # first principle component
cr_fmadjmeta = rank(cr_fmadjmeta)/ncol(cr_fmobsmat)
cr_fmadjmeta = 1-cr_fmadjmeta
cr_fmadjord = order(cr_fmadjmeta)

pdf('fmmeta2.pdf')
# see if it the raw or adjusted metagenes have any difference
plot(1-cr_fmadjmeta, -cr_fmrawmeta, pch=20, main='FM metagene comparison (Creighton data)', ylab='Raw metagene value', xlab='Adjusted metagene value')
cr_fmdend= hclust(dist(cr_fmobsmat_adj))

main = "FM metagene (Creighton data)"
metaplot3(cr_fmobsmat_adj, cr_fmadjmeta, crclin, name = main)

dev.off()

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

pdf('fmmeta3.pdf')

main='FM metagene (Creighton transformed)'
metaplot3(cr_fmobsmat_adj, cr_fmtransmeta, crclin, name = main)

dev.off()

###############################################################################
## Try it on ICGC data

## Pull out all the genes common in Creighton obesity genes and ICGC data:
cancertypes = c('BLCA', 'CESC', 'COAD', 'KIRP', 'LIHC', 'READ', 'SKCM', 'UCEC')
for (i in 1:length(cancertypes)) {
	seq = get(cancertypes[i])
	seq = seq[,common_genes]
	seq = t(seq) ## make it genes by samples for standardisation
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
	t = rank(t)/length(t)
	t = 1-t
	assign(paste(cancertypes[i], 'fmadjmeta',sep=''), t)
}

pdf('fmtcga.pdf')
## Check if the metagene correlates with sample gene expression and/or BMI:
for (i in 1:length(cancertypes)) {
	dat = get(paste('fm', cancertypes[i], sep=''))
	meta = get(paste(cancertypes[i], 'fmadjmeta', sep=''))
	bmi = get(paste(cancertypes[i], 'bmi', sep=''))
	txt = paste('FM metagene (', cancertypes[i], sep = "")
	txt = paste(txt, ')', sep = "")
	metaplot3(dat, meta, bmi, name = txt)
}
dev.off()

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
resobsgenes = restt[restt$P.Value < 0.01,]
length(which( cr_obsgene %in% rownames(resobsgenes ) )) ## 216
## What about in top 799 genes?
resobsgenes = rownames(resobsgenes[1:799,])
rescrolgenes = resobsgenes[which(resobsgenes %in% cr_obsgene)]
length(rescrolgenes) ## 168

pdf('venn1.pdf')
## make venn diagram with res, raw, and cr obesity genes:
setlist = list(Residual = resobsgenes, My_genes = rawobsgenes, Creighton = cr_obsgene)
OLlist = overLapper(setlist=setlist, sep='_', type='vennsets')
counts = sapply(OLlist$Venn_List, length)
vennPlot(counts = counts)
dev.off()

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
caobsgenes = catt[catt$P.Value < 0.01,]
length(which( cr_obsgene %in% rownames(caobsgenes ) )) ## 271
## What about in top 799 genes?
caobsgenes = rownames(caobsgenes[1:799,])
cacrolgenes = caobsgenes[which(caobsgenes %in% cr_obsgene)]
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
caresobsgenes = carestt[carestt$P.Value < 0.01,]
length(which( cr_obsgene %in% rownames(caresobsgenes ) )) ## 166
## What about in top 799 genes?
caresobsgenes = rownames(caresobsgenes[1:799,])
carescrolgenes = caresobsgenes[which(caresobsgenes %in% cr_obsgene)]
length(carescrolgenes) ## 92

pdf('venn2.pdf')
## make a Venn diagram with ca, cares, and cr obesity genes:
setlist = list(Caucasian_residual = caresobsgenes , Caucasian= caobsgenes, Creighton = cr_obsgene)
OLlist = overLapper(setlist=setlist, sep='_', type='vennsets')
counts = sapply(OLlist$Venn_List, length)
vennPlot(counts = counts)
dev.off()

## for each metagenes, identify the genes that are in the ICGC data and use these genes for validataion
cr_symmat = cr_raw
tmpgenes = mapIds(hgu133a.db, keys = rownames(cr_symmat), column = 'SYMBOL', keytype = "PROBEID", multiVals = 'first')
rownames(cr_symmat) = tmpgenes
cr_symmat = cr_symmat[-(which(is.na(rownames(cr_symmat)))),]
cr_symmat = collapseRows(cr_symmat, unique(rownames(cr_symmat)), unique(rownames(cr_symmat)))
cr_symmat = cr_symmat$datETcollapsed
dim(cr_symmat) #13031 genes
icgcgenes = colnames(BLCA)

allobsname = c("rawobsgenes","crolgenes","resobsgenes","rescrolgenes", "caobsgenes","cacrolgenes","caresobsgenes","carescrolgenes")
metalength = matrix(1,8)
rownames(metalength) = allobsname
for (i in 1:length(allobsname)) {
	tmp = get(allobsname[i])
	tmp = mapIds(hgu133a.db, keys = tmp, column = 'SYMBOL', keytype = "PROBEID", multiVals = 'first')
	tmp = unique(tmp)
	if (length(which(is.na(tmp))) > 0){
		tmp = tmp[-which(is.na(tmp))]
	}
	tmp = tmp[which(tmp %in% rownames(cr_symmat))]
	tmp = tmp[which(tmp %in% icgcgenes)]
	assign(allobsname[i], tmp)
	metalength[i,1] = length(tmp)
}

metalength
##rawobsgenes     678
##crolgenes       199
##resobsgenes     655
##rescrolgenes    147
##caobsgenes      657
##cacrolgenes     128
##caresobsgenes   651
##carescrolgenes  86

# need to check if the metagenes are going in the same direction

# find the genes that are common across all of the metagenes:
checkgenes= c(rawobsgenes,crolgenes,resobsgenes,rescrolgenes, caobsgenes,cacrolgenes,caresobsgenes,carescrolgenes)
checkgenes = table(checkgenes)[which(table(checkgenes) == 8)]
checkgenes = names(checkgenes)

# check for the direction of each of the metagene, using the common genes
# identified above
pdf('metadirection.pdf')
for (i in 1:length(allobsname)) {
	gene = get(allobsname[i])
	mat = cr_symmat[gene,]
	mat = t(apply(mat, 1, function(x) (x-mean(x))/sd(x)))
	mat[mat < -3] = -3
	mat[mat > 3] = 3
	tmpsvd = svd(mat)
	tmpmeta = rank(tmpsvd$v[,1])/ncol(mat)
	tmpmeta = 1-tmpmeta
	mat = mat[checkgenes,]
	ord = order(tmpmeta)
	heatmap.2(mat[,ord], trace='none',scale='none', col='bluered', ColSideColors = bluered(length(tmpmeta))[rank(tmpmeta)][ord], Colv=NA, main=namelist[i])
}
dev.off()

# the seventh metagene (caucasian+residual) must be flipped so all the metagenes are going in the same direction

## validate these metagenes in creighton data:
namelist = c("Raw metagene","CR overlap metagene (raw)","Residual metagene","CR overlap metagene (residual)","Raw metagene (Caucasian/raw)","CR overlap metagene (Caucasian/raw)","Residual metagene (Caucasian/residual)","CR overlap metagene (Caucasian/residual)")
metalist = list()
pdf('degmetacr.pdf')
for (i in 1:length(allobsname)) {
	gene = get(allobsname[i])
	mat = cr_symmat[gene,]
	mat = t(apply(mat, 1, function(x) (x-mean(x))/sd(x)))
	mat[mat < -3] = -3
	mat[mat > 3] = 3
	tmpsvd = svd(mat)
	tmpmeta = rank(tmpsvd$v[,1])/ncol(mat)
	if (i == 7) {
		tmpmeta = 1-tmpmeta
	}
	metaplot3(mat, tmpmeta, crclin, name = namelist[i])
	metalist[[i]] = tmpmeta

	txt = gsub("genes", "transmat", allobsname[i])
	trmat = diag(1/tmpsvd$d) %*% t(tmpsvd$u)
	assign(txt, trmat)
}
dev.off()

pdf('cor.pdf')
## make correlation matrix
x = matrix(1,103)
for(i in 1:length(metalist)) {
	x = cbind(x, metalist[[i]])
}
x = x[,-1]
colnames(x) = gsub('genes', '', allobsname)
cormat = cor(x, method = 'pearson')
cormatadj = x[,-which(colnames(x) == "caresobs")]
cormatadj = cor(cormatadj, method = 'pearson')
heatmap.2(cormat, trace = 'none', scale='none', col='bluered', cexRow = 1.0, cexCol = 1.0)
heatmap.2(cormatadj, trace = 'none', scale='none', col='bluered', cexRow = 1.0, cexCol = 1.0)
dev.off()

## check if this works on ICGC data:
alltransname = gsub('genes', 'transmat', allobsname)
bmimetalist = list()

for (i in 1:length(allobsname)) {
	pdfname = gsub('genes', 'meta.pdf', allobsname[i])
	pdf(pdfname)
	gene = get(allobsname[i])
	transmat = get(alltransname[i])
	tmpmeta = vector()
	for (j in 1:length(cancertypes)) {
		cancer = t(get(cancertypes[j]))
		mat = cancer[gene,]
		mat = standardise_data(mat)
		bmi = get(paste(cancertypes[j], 'bmi', sep=''))
		tmpsvd = t(transmat %*% mat)
		#tmpsvd = rank(tmpsvd[,1])/ncol(mat)
		tmpsvd = tmpsvd[,1]
		tmpsvd = rank(tmpsvd)/length(tmpsvd)
		if (i == 7) {
			tmpsvd = 1-tmpsvd
		}
		main = paste(namelist[i], '(')
		main = paste(main, cancertypes[j], sep = '')
		main = paste(main, ')', sep = '')
		metaplot3(mat, tmpsvd, bmi, name = main)
		tmpmeta = c(tmpmeta, tmpsvd)
	}
	bmimetalist[[i]] = tmpmeta
	dev.off()
}
bmimetalist = as.data.frame(bmimetalist)
colnames(bmimetalist) = gsub('genes', '', allobsname)

###############################################################################
## Gatza pathway stuff

## import pathway gene list from the Gatza paper
files = readLines('./gatzagenelist/pathlist.txt')
files = paste('./gatzagenelist/', files, sep='')
for (i in 1:length(files)) {
	txt = gsub('./gatzagenelist/','',files[i])
	txt = gsub('.txt','',txt)
	genes = read.csv(files[i])
	genes = as.vector(genes[,1])
	assign(txt, genes)
}
paths = gsub('./gatzagenelist/','',files)
paths = gsub('.txt','',paths)

## convert the gene probes into gene symbols:
pathlength = matrix(1,length(paths))
rownames(pathlength) = paths
for (i in 1:length(paths)) {
	pathgenes = get(paths[i])
	pathgenes = mapIds(hgu133a.db, keys = pathgenes, column = 'SYMBOL', keytype = "PROBEID", multiVals = 'first')
	pathgenes = unique(pathgenes)
	if (length(which(is.na(pathgenes))) > 0){
		pathgenes = pathgenes[-which(is.na(pathgenes))]
	}
	pathgenes = pathgenes[which(pathgenes %in% colnames(BLCA))]
	assign(paths[i], pathgenes)
	pathlength[i,1] = length(pathgenes)
}
pathlength
#			no. of probes	no. common probes
#akt_probes    206				191
#bcat_probes    77  			 75
#e2f1_probes   128  			120
#egfr_probes   419  			399
#er_probes     102  			 97
#her2_probes   212  			199
#ifna_probes    82  			 81
#ifng_probes    88  			 84
#myc_probes    425  			394
#p53_probes    218  			206
#p63_probes    277  			254
#pi3k_probes   220  			209
#pr_probes     212  			202
#ras_probes    300  			281
#src_probes     81  			 77
#stat3_probes  107  			 98
#tgfb_probes   101  			 93
#tnfa_probes    95  			 90

# need to check if the metagenes are going in the same direction

# find the genes that are common across all of the metagenes:
checkgzgenes= c(akt_probes, bcat_probes, e2f1_probes, egfr_probes, er_probes, her2_probes, ifna_probes, ifng_probes, myc_probes, p53_probes, p63_probes, pi3k_probes, pr_probes, ras_probes, src_probes, stat3_probes, tgfb_probes, tnfa_probes)
checkgzgenes = table(checkgzgenes)[which(table(checkgzgenes) == 8)]
checkgzgenes = names(checkgzgenes)

# check for the direction of each of the gatza metagene

# list of genes related to/representing the pathway:
checkgene = c('AKT2', 'CCNE1', 'MYB', 'EGFR', 'ESR1', 'ERBB2', 'IRF7', 'IRF1', 'MYC', 'CASP3', 'TP63', 'TP63', 'TGFA', 'HRAS', 'SRC', 'STAT3', 'PDGFB', 'TNFAIP2')

# make data matrix for the heatmap:
matheat = cr_symmat[checkgzgenes,]
matheat = t(apply(matheat, 1, function(x) (x-mean(x))/sd(x)))
matheat[matheat < -3] = -3
matheat[matheat > 3] = 3

pdf('gatzametadirection.pdf')
for (i in 1:length(paths)) {
	gene = get(paths[i])
	mat = cr_symmat[gene,]
	mat = t(apply(mat, 1, function(x) (x-mean(x))/sd(x)))
	mat[mat < -3] = -3
	mat[mat > 3] = 3
	tmpsvd = svd(mat)
	tmpmeta = rank(tmpsvd$v[,1])/ncol(mat)
	tmpmeta = 1-tmpmeta
	ord = order(tmpmeta)
	ind = which(gene %in% checkgene[i])
	matheat= cr_symmat[gene,]
	matheat = matheat[ind:(ind+1),]
	matheat = t(apply(matheat, 1, function(x) (x-mean(x))/sd(x)))
	matheat[matheat < -3] = -3
	matheat[matheat > 3] = 3
	heatmap.2(matheat[,ord], trace='none',scale='none', col='bluered', ColSideColors = bluered(length(tmpmeta))[rank(tmpmeta)][ord], Colv=NA, main=paths[i], cexRow=1.0)
}
dev.off()

# check to see if the gatza metagenes and BMI metagenes are going in the same
# direction:
#pdf('bmiandgatzametadirection.pdf')
#for (i in 1:length(allobsname)) {
#	gene = get(allobsname[i])
#	mat = cr_symmat[gene,]
#	mat = t(apply(mat, 1, function(x) (x-mean(x))/sd(x)))
#	mat[mat < -3] = -3
#	mat[mat > 3] = 3
#	tmpsvd = svd(mat)
#	tmpmeta = rank(tmpsvd$v[,1])/ncol(mat)
#	if (i == 7) {
#		tmpmeta = 1-tmpmeta
#	}
#	ord = order(tmpmeta)
#	heatmap.2(matheat[,ord], trace='none',scale='none', col='bluered', ColSideColors = bluered(length(tmpmeta))[rank(tmpmeta)][ord], Colv=NA, main=namelist[i])
#}
#dev.off()

# BMI metagenes were all going in the same direction as the gatza metagenes,
# and these metagenes from Gatza paper needs to be flipped:
mgflip = c('akt_probes', 'e2f1_probes', 'her2_probes', 'ifna_probes', 'ifng_probes', 'p63_probes', 'tgfb_probes')

# Check if it works:
pdf('gatzametadirectioncheck.pdf')
for (i in 1:length(paths)) {
	gene = get(paths[i])
	mat = cr_symmat[gene,]
	mat = t(apply(mat, 1, function(x) (x-mean(x))/sd(x)))
	mat[mat < -3] = -3
	mat[mat > 3] = 3
	tmpsvd = svd(mat)
	tmpmeta = rank(tmpsvd$v[,1])/ncol(mat)
	if (paths[i] %in% mgflip) {
		tmpmeta = -tmpmeta
	}
	tmpmeta = 1-tmpmeta
	ord = order(tmpmeta)
	ind = which(gene %in% checkgene[i])
	matheat= cr_symmat[gene,]
	matheat = matheat[ind:(ind+1),]
	matheat = t(apply(matheat, 1, function(x) (x-mean(x))/sd(x)))
	matheat[matheat < -3] = -3
	matheat[matheat > 3] = 3
	heatmap.2(matheat[,ord], trace='none',scale='none', col='bluered', ColSideColors = bluered(length(tmpmeta))[rank(tmpmeta)][ord], Colv=NA, main=paths[i], cexRow=1)
}
dev.off()

## make metagene with Gatza pathways in ICGC data samples:
gatzametalist = list()
for (i in 1:length(paths)) {
	pdfname = gsub('_probes', 'meta.pdf', paths[i])
	pdf(pdfname)
	gene = get(paths[i])
	mat = cr_symmat[gene,]
	mat = t(apply(mat, 1, function(x) (x-mean(x))/sd(x)))
	mat[mat < -3] = -3
	mat[mat > 3] = 3
	tmpsvd = svd(mat)
	trmat = diag(1/tmpsvd$d) %*% t(tmpsvd$u)
	tmpmeta = vector()
	for (j in 1:length(cancertypes)) {
		cancer = t(get(cancertypes[j]))
		mat = cancer[gene,]
		mat = standardise_data(mat)
		bmi = get(paste(cancertypes[j], 'bmi', sep=''))
		tmpsvd = t(trmat %*% mat)
		tmpsvd = tmpsvd[,1]
		tmpsvd = rank(tmpsvd)/length(tmpsvd)
		if (paths[i] %in% mgflip) {
			tmpsvd = 1-tmpsvd
		}
		main = gsub('_probes', '', paths[i])
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

## stick the BMI metagenes and Gatza metagenes together:
allmetagenes = cbind(bmimetalist, gatzametalist)
allmetagenes = t(data.matrix(allmetagenes))

## need to make a combined clinical data for all the cancer types
x = c('bcr_patient_barcode', 'tumor_tissue_site', 'weight', 'height')
allclin = BLCAclin[,x]
for(i in 2:length(cancertypes)) {
	txt = paste(cancertypes[i], 'clin', sep='')
	clin = get(txt)
	clin = clin[,x]
	allclin = rbind(allclin, clin)
}
rownames(allclin) = allclin$bcr_patient_barcode
allclin = allclin[colnames(allmetagenes),-1]
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

allcol = cbind(BMIStatus=bmicol, CancerType=cancertypecol)

pdf(file='gatzabmimeta.pdf',width=14, height=7)
heatmap.2x(allmetagenes, scale='none', trace='none', col='bluered', ColSideColors=cancertypecol)
heatmap.2x(allmetagenes[, 1:1872], scale='none', trace='none', col='bluered', ColSideColors=cancertypecol, Colv=F)
bmiord = order(bmicol)
heatmap3(allmetagenes[, bmiord], scale='none', ColSideColors=allcol, Colv=NA)
heatmap3(allmetagenes[,1:1872], scale='none', ColSideColors=bmicol, Colv=NA)
#heatmap.2x(allmetagenes, scale='none', trace='none', col="bluered", dendrogram="both")
dev.off()

pdf('allmetacor.pdf')
x = cor(t(allmetagenes), method='pearson')
heatmap.2(x, trace = 'none', scale='none', col='bluered', main="Pearson correlation")
x = cor(t(allmetagenes), method='spearman')
heatmap.2(x, trace = 'none', scale='none', col='bluered', main='Spearman correlation')
x = cor(gatzametalist, method='spearman')
heatmap.2(x, trace = 'none', scale='none', col='bluered', main='Spearman correlation')
dev.off()





















devtools::install_github("TomKellyGenetics/heatmap.2x", ref="supr")

# This example takes a subtype "#EF559E" subset and splits the heatmap be mutation status "CDH!_Mut"

bb3_Stat<-bb3[,bb3[8,]=="#EF559E" & cut == 1 & is.na(CDH1_Mt3)==F] ##bb3 is my column colour matrix of clinical and mutation data
dim(dataset)
## tree for genes
tree_exprSL_voom_corr_dist<-as.dendrogram(hclust(as.dist(1-cor(t(dataset)))))
## split columns by mutation status
data_low<-as.matrix(Data_Matrix_BRCA_RNASeq_Voom_GatzaYR[,bb3[8,]=="#EF559E" & CDH1_Mt3==1 & is.na(CDH1_Mt3)==F & cut == 1])
data_high<-as.matrix(Data_Matrix_BRCA_RNASeq_Voom_GatzaYR[,bb3[8,]=="#EF559E" & CDH1_Mt3==0 & is.na(CDH1_Mt3)==F & cut == 1])
## trees for each split dataset (correlation distance)
dist_low<-dist(as.dist(1-cor(data_low)))
hc_low<-hclust(dist_low,method='complete')
dist_high<-dist(as.dist(1-cor(data_high)))
hc_high<-hclust(dist_high,method='complete')
## join trees together
hc<-merge(as.dendrogram(hc_low), as.dendrogram(hc_high), height=max(hc_low$height, hc_high$height)+1)
#rr<-RowCols[match(rownames(dataset), rownames(Data_Matrix_BRCA_RNASeq_Voom_GatzaYR)),]
# plot clusters of genes
Cluster<-cutree(as.hclust(tree_exprSL_voom_corr_dist),4)
ColCluster<-c("red", "blue", "green", "orange")[Cluster]
#rr<-cbind(rr, ColCluster)  ## rr is my row colour data for genes
heatmap.mik.mod2(as.matrix(dataset[,match(labels(hc), colnames(dataset))]), scale='none', trace='none', col=bluered(50), ColSideColors=bb3_Stat[,match(labels(hc), colnames(dataset))], Colv=hc, Rowv=tree_exprSL_voom_corr_dist, margin=c(12, 12), dendrogram='both', main = "TCGA Breast Gene Expression Gatza 2011", xlab = "Sample", ylab = "Pathway", cexCol=1.15, cexRow=1.15, keysize=2.25)

## add legend to heatmap
legend("topleft",legend=c(rep("", 11), "Normal", "Tumour", "Metastasis", "", "Ductal", "Lobular", "", "Stage 1", "Stage 2", "Stage 3", "Stage 4", "", "Positive","Negative", "", "Basal","Her2","LumA","LumB","Normal","","Somatic Mutation","Somatic Mutation (Slient)","Wildtype","Somatic Mutation (Gene 1)","Somatic Mutation (Gene 2)", "Somatic Mutation (Gene 3)", "Somatic Mutation (Gene 1 Silent)","Somatic Mutation (Gene 2 Silent)", "Somatic Mutation (Gene 3 Silent)", "Somatic Mutation (Gene 1 and 2)", "", "CDH1 Low",  "CDH1 High"),

tree_exprSL_voom_eu_dist<-as.dendrogram(hclust(dist(dataset)))
data_low<-as.matrix(Data_Matrix_BRCA_RNASeq_Voom_GatzaYR[,bb3[8,]=="#EF559E" & CDH1_Mt3==1 & is.na(CDH1_Mt3)==F & cut == 1])
data_high<-as.matrix(Data_Matrix_BRCA_RNASeq_Voom_GatzaYR[,bb3[8,]=="#EF559E" & CDH1_Mt3==0 & is.na(CDH1_Mt3)==F & cut == 1])
dist_low<-dist(t(data_low), method = "euclidean")
hc_low<-hclust(dist_low,method='complete')
dist_high<-dist(t(data_high), method = "euclidean")
hc_high<-hclust(dist_high,method='complete')
hc<-merge(as.dendrogram(hc_low), as.dendrogram(hc_high))
hc<-merge(as.dendrogram(hc_low), as.dendrogram(hc_high), height=max(hc_low$height, hc_high$height)+1)#rr<-RowCols[match(rownames(dataset), rownames(Data_Matrix_BRCA_RNASeq_Voom_GatzaYR)),]
Cluster<-cutree(as.hclust(tree_exprSL_voom_eu_dist),4)
ColCluster<-c("red", "blue", "green", "orange")[Cluster]
#rr[, ncol(rr)]<-ColCluster
heatmap.mik.mod2(as.matrix(dataset[,match(labels(hc), colnames(dataset))]), scale='none', trace='none', col=bluered(50), ColSideColors=bb3_Stat[,match(labels(hc), colnames(dataset))], Colv=hc, Rowv=tree_exprSL_voom_eu_dist, margin=c(12, 12), dendrogram='both', main = "TCGA Breast Gene Expression Gatza 2011", xlab = "Sample", ylab = "Pathway", cexCol=1.15, cexRow=1.15, keysize=2.25)







###############################################################################
## Common genes across ICGC data stuff

## use codes from task11?



###############################################################################
## Pathway enrichment stuff

## TODO: do the enrichment analysis with reactome and kegg






###############################################################################
## Extra stuff

## TODO: Gatza signatures







