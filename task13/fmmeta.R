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

pdf('pdf/fmmeta1.pdf')

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

pdf('pdf/fmmeta2.pdf')
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

pdf('pdf/fmmeta3.pdf')

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

pdf('pdf/fmtcga.pdf')
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


