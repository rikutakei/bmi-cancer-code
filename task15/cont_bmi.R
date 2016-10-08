## Get obesity-associated genes from the correlation of the genes with the BMI
## values, rather than the sample BMI status

# Use RMA normalised Creighton data:
mat = crrma
dim(mat) # 22283 103

# Correlate it with the BMI value of the samples:
bmicor = cor(t(mat), crclin$bmi, method = "spearman")
colnames(bmicor) = "correlation"

# Try it with RMA normalised, standardised data:
stdmat = crstdrma
dim(mat) # 22283 103

# Correlate it with the BMI value of the samples:
bmicorstd = cor(t(stdmat), crclin$bmi, method = "spearman")
colnames(bmicorstd) = "correlation"

# The values didn't have any difference in the results, so use either data:
summary(bmicor)
summary(bmicorstd)

# Do some data scrambling
set.seed(1) # set seed for reproducibility

# Randomly pick samples:
tmpbmi = sample(crclin$bmi, length(crclin$bmi))

# Use the random BMI to find correlation you would get by chance:
sismres1 = cor(t(stdmat), tmpbmi, method = "spearman")

# Repeat this 1000 times to get a distribution of correlation values:
for (i in 2:1000) {
	tmpbmi = sample(crclin$bmi, length(crclin$bmi))
	tmpres = cor(t(stdmat), tmpbmi, method = "spearman")
	sismres1 = cbind(sismres1, tmpres)
	if ((i %% 50) == 0)
	print(paste( 'finished ', i))
}

## Get gene that are significant:
## First get the p-value of the genes:
v = vector()
for (i in 1:nrow(sismres1)) {
	## probability of observing the true value by chance:
	tmp1 = sum(bmicorstd[i,] < sismres1[i,])/ncol(sismres1)
	tmp2 = sum(bmicorstd[i,] > sismres1[i,])/ncol(sismres1)
	val = min(tmp1, tmp2)
	v = c(v, val)
}
names(v) = rownames(bmicorstd)

## Pull out the genes with p <= 0.05:
ind = which(v <= 0.05)
length(ind) ## 4804 genes significant
contbmideg = bmicorstd[ind,]

## Order the genes based on the p-value:
ord = order(v[ind])
contbmideg = contbmideg[ord]

## Take the top 799 genes:
contbmigene = names(contbmideg)[1:799]

## Check if you get any significant genes with adjusted p-value:
adjv = p.adjust(v, method='BH', n = length(v))
ind = which(adjv <= 0.05)
length(ind) ## 162 genes significant

contbmiadjgene = rownames(bmicorstd)[ind]

###############################################################################
## Repeat this in other versions of Creighton data

## TODO: re-check the "get significant genes" sections for each of them
## Start with Caucasian-only data:
mat = camat

# Correlate it with the BMI value of the samples:
bmicor = cor(t(mat), caclin$bmi, method = "spearman")
colnames(bmicor) = "correlation"
summary(bmicor)

## Simulate with random bmi (n = 1000):
set.seed(1)
sismres2 = sim_bmi_cor(mat, caclin, n = 1000)
pval = get_pval(bmicor, sismres2)

## Get significant genes:
ind = which(pval <= 0.05)
length(ind) ## 2969 genes significant
contbmideg = bmicor[ind,]
ord = order(pval[ind])
contbmideg = contbmideg[ord]
cacontbmigene = names(contbmideg)[1:799] ## Take the top 799 genes:

## Check if you get any significant genes with adjusted p-value:
adjv = p.adjust(pval, method='BH', n = length(pval))
ind = which(adjv <= 0.05)
length(ind) ## 46 genes significant

cacontbmiadjgene = rownames(bmicor)[ind]

## Repeat with normal residual data:
mat = residuals

# Correlate it with the BMI value of the samples:
bmicor = cor(t(mat), crclin$bmi, method = "spearman")
colnames(bmicor) = "correlation"
summary(bmicor)

## Simulate with random bmi (n = 1000):
set.seed(1)
sismres3 = sim_bmi_cor(mat, crclin, n = 1000)
pval = get_pval(bmicor, sismres3)

## Get significant genes:
ind = which(pval <= 0.05)
length(ind) ## 2976 genes significant
contbmideg = bmicor[ind,]
ord = order(pval[ind])
contbmideg = contbmideg[ord]
rescontbmigene = names(contbmideg)[1:799] ## Take the top 799 genes:

## Check if you get any significant genes with adjusted p-pvalalue:
adjv = p.adjust(pval, method='BH', n = length(pval))
ind = which(adjv <= 0.05)
length(ind) ## 45 genes significant

rescontbmiadjgene = rownames(bmicor)[ind]

## Repeat it with Caucasian-only residual data:
mat = caresiduals

# Correlate it with the BMI value of the samples:
bmicor = cor(t(mat), caclin$bmi, method = "spearman")
colnames(bmicor) = "correlation"
summary(bmicor)

## Simulate with random bmi (n = 1000):
set.seed(1)
sismres4 = sim_bmi_cor(mat, caclin, n = 1000)
pval = get_pval(bmicor, sismres4)

## Get significant genes:
ind = which(pval <= 0.05)
length(ind) ## 2164 genes significant
contbmideg = bmicor[ind,]
ord = order(pval[ind])
contbmideg = contbmideg[ord]
carescontbmigene = names(contbmideg)[1:799] ## Take the top 799 genes:

## Check if you get any significant genes with adjusted p-value:
adjv = p.adjust(pval, method='BH', n = length(pval))
ind = which(adjv <= 0.05)
length(ind) ## 16 genes significant

carescontbmiadjgene = rownames(bmicor)[ind]

###############################################################################
## Get gene symbol version of the continuous BMI genes:

contbmigenes = c('contbmigene', 'cacontbmigene', 'rescontbmigene', 'carescontbmigene', 'contbmiadjgene', 'cacontbmiadjgene', 'rescontbmiadjgene', 'carescontbmiadjgene')

for (i in 1:length(contbmigenes)) {
	tmpgenes = mapIds(hgu133a.db, keys = get(contbmigenes[i]), column = 'SYMBOL', keytype = "PROBEID", multiVals = 'first')
	tmpgenes = tmpgenes[which(!is.na(tmpgenes))]
	tmpgenes = unique(tmpgenes)
	tmpgenes = tmpgenes[which(tmpgenes %in% colnames(BLCAraw))]
	print(length(tmpgenes))
	assign(contbmigenes[i], tmpgenes)
}
# 670 contbmigene
# 663 cacontbmigene
# 666 rescontbmigene
# 662 carescontbmigene
# 144 contbmiadjgene
# 40 cacontbmiadjgene
# 41 rescontbmiadjgene
# 15 carescontbmiadjgene

###############################################################################
## Check if the continuous BMI genes correlate with gene expression and/or BMI
## in Creighton's data

namelist = c("Continuous BMI metagene (original)", "Continuous BMI metagene (Caucasian-only)", "Continuous BMI metagene (Residual)", "Continuous BMI metagene (Caucasian residual)", "Continuous BMI metagene (Adjusted p-value)", "Continuous BMI metagene (Adj. Caucasian-only)", "Continuous BMI metagene (Adj. Residual)", "Continuous BMI metagene (Adj. Caucasian residual)")

## Use RMA normalised and standardised Creighton data:
pdf('pdf/creighton_mg_results/contbmimeta_cr.pdf')
dendlist = list()
conttransmatlist = list()
for (i in 1:length(contbmigenes)) {
	genes = get(contbmigenes[i])
	mat = crsymrma[genes,]
	## Make metagene and plot it:
	contbmimeta = make_metagene(mat, flip = F, raw = F)
	main = namelist[i]
	mgheatmap(mat, contbmimeta, main = main)
	bmiplot(crclin, contbmimeta, main = main)
	## Get dendrogram:
	dendlist[[i]] = get_dend(mat)
	## Make transformation matrix:
	conttransmatlist[[i]] = make_trans_mat(mat)
}
dev.off()

###############################################################################
## Check if the continuous BMI genes correlate with gene expression and/or BMI
## in ICGC data

cancertypes = c('BLCA', 'CESC', 'COAD', 'KIRP', 'LIHC', 'READ', 'SKCM', 'UCEC')

for (j in 1:length(contbmigenes)) {
	pdfname = "pdf/creighton_mg_results/"
	pdfname = paste(pdfname, contbmigenes[j], sep = '')
	pdfname = paste(pdfname, 'ICGC.pdf', sep = '_')
	pdf(pdfname)
	gene = get(contbmigenes[j])
	transmat = conttransmatlist[[j]]
	dend = dendlist[[j]]
	for (i in 1:length(cancertypes)) {
		txt = paste(cancertypes[i], 'raw', sep='')
		cancer = t(get(txt))
		mat = cancer[gene,]
		mat = standardise_data(mat)
		clin = get(paste(cancertypes[i], 'clin', sep=''))
		cancermeta = make_transmeta(mat, transmat, flip = F)
		main = paste(namelist[j], '\n(', sep = '')
		main = paste(main, cancertypes[i], sep = '')
		main = paste(main, ')', sep = '')
		mgheatmap(mat, cancermeta, dend = dend, main = main)
		bmiplot(clin, cancermeta, main = main)
	}
	dev.off()
}

###############################################################################
## Check if the continuous BMI genes correlate with gene expression and/or BMI
## in Cris' data

pdf('pdf/creighton_mg_results/contbmimeta_cris.pdf')
for (i in 1:length(contbmigenes)) {
	genes = get(contbmigenes[i])
	mat = crissymrma[genes,]
	meta = make_transmeta(mat, conttransmatlist[[i]], flip = F)
	dend = dendlist[[i]]
	main = "Continuous BMI metagene (Cris' data)"
	main = paste(namelist[i], "\n(Cris' data)", sep = '')
	mgheatmap(mat, meta, dend = dend, main = main)
	bmiplot(crisclin, meta, main = main)
}
dev.off()

