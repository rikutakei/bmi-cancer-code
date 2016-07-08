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

###############################################################################
## Gatza metagene in Creighton data:

# Check the correlation between the  metagenes created from raw and standardised
pdf('pdf/crgatzarawvsstdmeta.pdf')
plot_raw_vs_std(cr_symmat, paths, main = 'Creighton')
dev.off()

# need to check if the metagenes are going in the same direction

# list of genes related to/representing the pathway:
checkgene = c('AKT1', 'CTNNB1', 'E2F1', 'EGFR', 'ESR1', 'ERBB2', 'IFNA1', 'IFNG', 'MYC', 'TP53', 'TP63', 'PIK3CA', 'PGR', 'HRAS', 'SRC', 'STAT3', 'TGFB1', 'TNF')

# make data matrix for the heatmap:
matheat = cr_symmat[checkgene,]
matheat = t(apply(matheat, 1, function(x) (x-mean(x))/sd(x)))
matheat[matheat < -3] = -3
matheat[matheat > 3] = 3

# make row colours for heatmap
col = make_col(as.vector(crclin$ERstatus), continuous=F)
col = rbind(ER = col, PR = make_col(as.vector(crclin$PRstatus), continuous=F), HER2 = make_col(as.vector(crclin$HER2status), continuous=F), LN = make_col(as.vector(crclin$LNstatus), continuous=F))

pdf('pdf/gatzametadirection.pdf')
crgatzametalist = list()
for (i in 1:length(paths)) {
	gene = get(paths[i])
	mat = cr_symmat[gene,]
	mat = t(apply(mat, 1, function(x) (x-mean(x))/sd(x)))
	mat[mat < -3] = -3
	mat[mat > 3] = 3
	tmpsvd = svd(mat)
	tmpmeta = rank(tmpsvd$v[,1])/ncol(mat)
	ord = order(tmpmeta)
	#heatmap.2(matheat[,ord], trace='none',scale='none', col='bluered', ColSideColors = bluered(length(tmpmeta))[rank(tmpmeta)][ord], Colv=NA, main=paths[i], cexRow=1.0)
	col1 = bluered(length(tmpmeta))[rank(tmpmeta)]
	col2 = bluered(length(tmpmeta))[rank(cr_symmat[checkgene[i],])]
	tmpcol = rbind(col2, col, meta=col1)
	txt = checkgene[i]
	rownames(tmpcol)[1] = txt
	heatmap.2x(matheat[,ord], trace='none',scale='none', col='bluered', ColSideColors = tmpcol[,ord], Colv=NA, main=paths[i], cexRow=1.0)
	crgatzametalist[[i]] = tmpmeta
}
dev.off()
crgatzametalist = as.data.frame(crgatzametalist)
colnames(crgatzametalist) = gsub('_probes', '', paths)
crgatzametalist = t(as.matrix(crgatzametalist))

# calculate the correlation of the metagenes with their representing genes:
gatzacorval = c()
for (i in 1:nrow(crgatzametalist)) {
	tmpcor = cor(crgatzametalist[i,], cr_symmat[checkgene[i],])
	gatzacorval = c(gatzacorval, tmpcor)
}
names(gatzacorval) = rownames(crgatzametalist)

# BMI metagenes were all going in the same direction as the gatza metagenes,
# and these metagenes from Gatza paper needs to be flipped:
mgflip = c('akt_probes', 'e2f1_probes', 'egfr_probes', 'er_probes', 'ifna_probes', 'myc_probes', 'p53_probes', 'ras_probes', 'src_probes', 'stat3_probes', 'tnfa_probes')
rawmgflip = c('akt_probes', 'e2f1_probes', 'egfr_probes', 'er_probes', 'ifna_probes', 'myc_probes', 'p53_probes', 'ras_probes', 'src_probes', 'stat3_probes', 'tnfa_probes')

# Check if it works:
pdf('pdf/gatzametadirectioncheck.pdf')
crgatzametalist = list()
for (i in 1:length(paths)) {
	gene = get(paths[i])
	mat = cr_symmat[gene,]
	mat = t(apply(mat, 1, function(x) (x-mean(x))/sd(x)))
	tmpsvd = svd(mat)
	tmpmeta = rank(tmpsvd$v[,1])/ncol(mat)
	if (paths[i] %in% mgflip) {
		tmpmeta = 1-tmpmeta
	}
	ord = order(tmpmeta)
	#heatmap.2(matheat[,ord], trace='none',scale='none', col='bluered', ColSideColors = bluered(length(tmpmeta))[rank(tmpmeta)][ord], Colv=NA, main=paths[i], cexRow=1.0)
	col1 = bluered(length(tmpmeta))[rank(tmpmeta)]
	col2 = bluered(length(tmpmeta))[rank(cr_symmat[checkgene[i],])]
	tmpcol = rbind(col2, col, meta=col1)
	txt = checkgene[i]
	rownames(tmpcol)[1] = txt
	heatmap.2x(matheat[,ord], trace='none',scale='none', col='bluered', ColSideColors = tmpcol[,ord], Colv=NA, main=paths[i], cexRow=1.0)
	crgatzametalist[[i]] = tmpmeta
}
crgatzametalist = as.data.frame(crgatzametalist)
colnames(crgatzametalist) = gsub('_probes', '', paths)
crgatzametalist = t(as.matrix(crgatzametalist))
heatmap.2x(crgatzametalist, trace='none',scale='none', col='bluered', ColSideColors = col, main='All Gatza Pathway Metagenes', cexRow=1.0)
dev.off()

gatzacor = cor(t(crgatzametalist), method='pearson')
gatzacor2 = cor(t(crgatzametalist), method='spearman')

gatzaord = c('er', 'pr', 'p53', 'bcat', 'e2f1', 'pi3k', 'myc', 'ras', 'ifna', 'ifng', 'akt', 'p63', 'src', 'her2', 'egfr', 'tgfb', 'stat3', 'tnfa')

# check if I get similar clustering as Gatza paper:
pdf('pdf/gatzacheck.pdf')
#pdf('pdf/gatzachecknoflip.pdf')
heatmap.2(crgatzametalist, trace='none',scale='none', col=matlab.like, cexRow=1)
heatmap.2(crgatzametalist[gatzaord,], trace='none',scale='none', col=matlab.like, cexRow=1, dendrogram='none', Rowv=F)
ord = hclust(dist(crgatzametalist))
dend = as.dendrogram(ord)
dend = reorder(dend, rowMeans(crgatzametalist))
#ord = rev(ord)
heatmap.2(gatzacor, trace='none',scale='none', col=matlab.like, cexRow=1, main='pearson', Rowv=dend, Colv=dend)
heatmap.2(gatzacor2, trace='none',scale='none', col=matlab.like, cexRow=1, main='spearman', Rowv=dend, Colv=dend)
#heatmap.2(gatzacor[gatzaord,gatzaord], trace='none',scale='none', col=matlab.like, cexRow=1, dendrogram='none', Rowv=F, Colv=F)
dev.off()

###############################################################################
## Gatza metagene in Fuentes-Mattei data:

# Check the correlation between the  metagenes created from raw and standardised
pdf('pdf/fmgatzarawvsstdmeta.pdf')
plot_raw_vs_std(fm_symmat, paths, main='Fuentes-Mattei')
dev.off()

# make data matrix for the heatmap:
matheat = fm_symmat[checkgene,]
matheat = t(apply(matheat, 1, function(x) (x-mean(x))/sd(x)))
matheat[matheat < -3] = -3
matheat[matheat > 3] = 3

# make row colours for heatmap
tmpcol = c('blue', 'red', 'white')
ER = make_col(as.vector(fmclin$characteristics..ER_status), continuous=F, colours = tmpcol)
PR = make_col(as.vector(fmclin$characteristics..PR_status), continuous=F, colours = tmpcol)
HER2 = make_col(as.vector(fmclin$HER2.Status), continuous=F, colours = tmpcol)
tmpcol = brewer.pal(3, 'YlOrRd')
Grade = make_col(as.vector(fmclin$BMNgrd), continuous=F, colours = tmpcol)

col = rbind(ER = ER, PR = PR, HER2 = HER2, Grade = Grade)

mgflip = c('bcat_probes', 'egfr_probes', 'er_probes', 'her2_probes', 'myc_probes', 'pr_probes', 'ras_probes', 'src_probes', 'tnfa_probes')

pdf('pdf/gatzametadirectionfm.pdf')
fmgatzametalist = list()
fmgatzametacor = vector()
for (i in 1:length(paths)) {
	gene = get(paths[i])
	mat = fm_symmat[gene,]
	mat = t(apply(mat, 1, function(x) (x-mean(x))/sd(x)))
	mat[mat < -3] = -3
	mat[mat > 3] = 3
	tmpsvd = svd(mat)
	tmpmeta = rank(tmpsvd$v[,1])/ncol(mat)
	if (paths[i] %in% mgflip) {
		tmpmeta = 1-tmpmeta
	}
	ord = order(tmpmeta)
	#heatmap.2(matheat[,ord], trace='none',scale='none', col='bluered', ColSideColors = bluered(length(tmpmeta))[rank(tmpmeta)][ord], Colv=NA, main=paths[i], cexRow=1.0)
	col1 = bluered(length(tmpmeta))[rank(tmpmeta)]
	col2 = bluered(length(tmpmeta))[rank(fm_symmat[checkgene[i],])]
	tmpcol = rbind(col2, col, meta=col1)
	txt = checkgene[i]
	rownames(tmpcol)[1] = txt
	heatmap.2x(matheat[,ord], trace='none',scale='none', col='bluered', ColSideColors = tmpcol[,ord], Colv=NA, main=paths[i], cexRow=1.0)
	fmgatzametalist[[i]] = tmpmeta
	cor = cor(tmpmeta, fm_symmat[checkgene[i],])
	fmgatzametacor = c(fmgatzametacor, cor)
}
fmgatzametalist = as.data.frame(fmgatzametalist)
colnames(fmgatzametalist) = gsub('_probes', '', paths)
fmgatzametalist = t(as.matrix(fmgatzametalist))
heatmap.2x(fmgatzametalist, trace='none',scale='none', col='bluered', ColSideColors = col, main='All Gatza Pathway Metagenes in FM data', cexRow=1.0)
names(fmgatzametacor) = paths
dev.off()

gatzacor = cor(t(fmgatzametalist), method='pearson')
gatzacor2 = cor(t(fmgatzametalist), method='spearman')

# check if I get similar clustering as Gatza paper:
pdf('pdf/gatzacheckfm.pdf')
heatmap.2x(fmgatzametalist, trace='none',scale='none', col=matlab.like, ColSideColors = col, cexRow=1)
ord = hclust(dist(fmgatzametalist))
dend = as.dendrogram(ord)
dend = reorder(dend, rowMeans(fmgatzametalist))
heatmap.2(gatzacor, trace='none',scale='none', col=matlab.like, cexRow=1, main='pearson', Rowv=dend, Colv=dend)
heatmap.2(gatzacor2, trace='none',scale='none', col=matlab.like, cexRow=1, main='spearman', Rowv=dend, Colv=dend)
dev.off()

###############################################################################
## Gatza metagene in Cris Print's breast cancer data:

# Check the correlation between the  metagenes created from raw and standardised
pdf('pdf/crisgatzarawvsstdmeta.pdf')
plot_raw_vs_std(cris_symmat, paths, main='Cris BC data')
dev.off()

# make data matrix for the heatmap:
matheat = cris_symmat[checkgene,]
matheat = t(apply(matheat, 1, function(x) (x-mean(x))/sd(x)))
matheat[matheat < -3] = -3
matheat[matheat > 3] = 3

# make row colours for heatmap
tmpcol = c('blue', 'red', 'white')
ER = make_col(as.vector(crisclin$ER.status), continuous=F, colours = tmpcol)
PR = make_col(as.vector(crisclin$PgR.status), continuous=F, colours = tmpcol)
LN = make_col(as.vector(crisclin$LN.status), continuous=F, colours = tmpcol)
tmpcol = brewer.pal(3, 'YlOrRd')
Grade = make_col(as.vector(crisclin$Grade), continuous=F, colours = tmpcol)
tmpcol = brewer.pal(6, 'Dark2')
SubType = make_col(as.vector(crisclin$subtype), continuous=F, colours = tmpcol)

col = rbind(ER = ER, PR = PR, LN = LN, Grade = Grade, SubType = SubType)

mgflip = c('bcat_probes', 'e2f1_probes', 'egfr_probes', 'ifng_probes', 'p63_probes', 'pr_probes', 'ras_probes', 'stat3_probes', 'tgfb_probes')

pdf('pdf/gatzametadirectioncris.pdf')
crisgatzametalist = list()
crisgatzametacor = vector()
for (i in 1:length(paths)) {
	gene = get(paths[i])
	mat = cris_symmat[gene,]
	mat = t(apply(mat, 1, function(x) (x-mean(x))/sd(x)))
	mat[mat < -3] = -3
	mat[mat > 3] = 3
	tmpsvd = svd(mat)
	tmpmeta = rank(tmpsvd$v[,1])/ncol(mat)
	if (paths[i] %in% mgflip) {
		tmpmeta = 1-tmpmeta
	}
	ord = order(tmpmeta)
	#heatmap.2(matheat[,ord], trace='none',scale='none', col='bluered', ColSideColors = bluered(length(tmpmeta))[rank(tmpmeta)][ord], Colv=NA, main=paths[i], cexRow=1.0)
	col1 = bluered(length(tmpmeta))[rank(tmpmeta)]
	col2 = bluered(length(tmpmeta))[rank(cris_symmat[checkgene[i],])]
	tmpcol = rbind(col2, col, meta=col1)
	txt = checkgene[i]
	rownames(tmpcol)[1] = txt
	heatmap.2x(matheat[,ord], trace='none',scale='none', col='bluered', ColSideColors = tmpcol[,ord], Colv=NA, main=paths[i], cexRow=1.0)
	crisgatzametalist[[i]] = tmpmeta
	cor = cor(tmpmeta, cris_symmat[checkgene[i],])
	crisgatzametacor = c(crisgatzametacor, cor)
}
crisgatzametalist = as.data.frame(crisgatzametalist)
colnames(crisgatzametalist) = gsub('_probes', '', paths)
crisgatzametalist = t(as.matrix(crisgatzametalist))
heatmap.2x(crisgatzametalist, trace='none',scale='none', col='bluered', ColSideColors = col, main='All Gatza Pathway Metagenes in Cris BC data', cexRow=1.0)
names(crisgatzametacor) = paths
dev.off()

gatzacor = cor(t(crisgatzametalist), method='pearson')
gatzacor2 = cor(t(crisgatzametalist), method='spearman')

gatzaord = c('er', 'pr', 'p53', 'bcat', 'e2f1', 'pi3k', 'myc', 'ras', 'ifna', 'ifng', 'akt', 'p63', 'src', 'her2', 'egfr', 'tgfb', 'stat3', 'tnfa')

# check if I get similar clustering as Gatza paper:
pdf('pdf/gatzacheckcris.pdf')
#pdf('pdf/gatzachecknoflip.pdf')
heatmap.2x(crisgatzametalist, trace='none',scale='none', col=matlab.like, ColSideColors = col, cexRow=1)
#heatmap.2x(crisgatzametalist[gatzaord,], trace='none',scale='none', col=matlab.like, ColSideColors = col, cexRow=1)
ord = hclust(dist(crisgatzametalist))
dend = as.dendrogram(ord)
dend = reorder(dend, rowMeans(crisgatzametalist))
#ord = rev(ord)
heatmap.2(gatzacor, trace='none',scale='none', col=matlab.like, cexRow=1, main='pearson', Rowv=dend, Colv=dend)
heatmap.2(gatzacor2, trace='none',scale='none', col=matlab.like, cexRow=1, main='spearman', Rowv=dend, Colv=dend)
#heatmap.2(gatzacor[gatzaord,gatzaord], trace='none',scale='none', col=matlab.like, cexRow=1, dendrogram='none', Rowv=F, Colv=F)
dev.off()

