###############################################################################
## Gene expression analysis using Cris Print's BC data

# Get patients with BMI information (see crismeta.R to see how I made crisobclin):
crisobclin = crisobclin[which(crisobclin$bmiStatus != "NA"), ]
samplenames = paste(crisobclin$GeoFile, crisobclin$CelFile, sep='_')
crisobmat = cris_raw[,samplenames]

# Standardise data
crisobmatstd = standardise_data(crisobmat, log=F)

# Make groups based on BMI Status
group = crisobclin$bmiStatus
group = ifelse(group=='obese', 'obese', 'non-obese')
design = model.matrix(~group)

# find DEGs between the obese vs non-obese groups (raw data)
ttraw = make_tt(crisobmat, design)
sum(ttraw$P.Value < 0.05)   # 2147 genes
sum(ttraw$P.Value < 0.01)   # 366 genes
sum(ttraw$adj.P.Val < 0.05) # 0 genes
sum(ttraw$adj.P.Val < 0.01) # 0 genes
length(which(2^ttraw$logFC > 1.2)) # 1013 genes have FC > 1.2

# find DEGs between the obese vs non-obese groups (standardised data)
ttstd = make_tt(crisobmatstd, design)
sum(ttstd$P.Value < 0.05)   # 2227 genes
sum(ttstd$P.Value < 0.01)   # 377 genes
sum(ttstd$adj.P.Val < 0.05) # 0 genes
sum(ttstd$adj.P.Val < 0.01) # 0 genes
length(which(2^ttstd$logFC > 1.2)) # 3481 genes have FC > 1.2

# use p < 0.01 genes:
crisobgenesraw = rownames(ttraw[ttraw$P.Value < 0.01,])
crisobgenesstd = rownames(ttstd[ttstd$P.Value < 0.01,])

# check in Cris' BC data:
# obesity genes from the raw data first:
pdf(file='pdf/crisdegmeta.pdf', width=7, height=7)
checkmat = crisobmat[crisobgenesraw, ]
tmpsvd = svd(checkmat)
tmpmeta = tmpsvd$v[,1]
tmpmeta = rank(tmpmeta)/length(tmpmeta)
checkmat = t(apply(checkmat, 1, function(x) (x-mean(x))/sd(x)))
stdsvd = svd(checkmat)
stdmeta = stdsvd$v[,1]
stdmeta = rank(stdmeta)/length(stdmeta)
checkmat[checkmat < -3] = -3
checkmat[checkmat > 3] = 3
metaplot3(checkmat, tmpmeta, crisobclin, name = "Raw metagene from Cris' BC data in \nCris' BC data (raw)")
metaplot3(checkmat, stdmeta, crisobclin, name = "Raw metagene from Cris' BC data in \nCris' BC data (std)")
# obesity genes from the standardised data:
checkmat = crisobmat[crisobgenesstd, ]
tmpsvd = svd(checkmat)
tmpmeta = tmpsvd$v[,1]
tmpmeta = rank(tmpmeta)/length(tmpmeta)
checkmat = t(apply(checkmat, 1, function(x) (x-mean(x))/sd(x)))
stdsvd = svd(checkmat)
stdmeta = stdsvd$v[,1]
stdmeta = rank(stdmeta)/length(stdmeta)
checkmat[checkmat < -3] = -3
checkmat[checkmat > 3] = 3
metaplot3(checkmat, tmpmeta, crisobclin, name = "Std. metagene from Cris' BC data in \nCris' BC data (raw)")
metaplot3(checkmat, stdmeta, crisobclin, name = "Std. metagene from Cris' BC data in \nCris' BC data (std)")
dev.off()

###############################################################################
# Need to convert gene probe ID into gene symbols, and check which genes are in ICGC data

crisgenesymraw = mapIds(hgu133plus2.db, keys = crisobgenesraw, column = 'SYMBOL', keytype = "PROBEID", multiVals = 'first')
crisgenesymraw = crisgenesymraw[which(!is.na(crisgenesymraw))]
crisgenesymraw = unique(as.vector(crisgenesymraw))
crisgenesymstd = mapIds(hgu133plus2.db, keys = crisobgenesstd, column = 'SYMBOL', keytype = "PROBEID", multiVals = 'first')
crisgenesymstd = crisgenesymstd[which(!is.na(crisgenesymstd))]
crisgenesymstd = unique(as.vector(crisgenesymstd))

# Just a quick check in Cris' data to see that it's still related to BMI:
# First convert the gene probes into gene symbols
crisobmatsym = crisobmat
genes = mapIds(hgu133plus2.db, keys = rownames(crisobmatsym), column = 'SYMBOL', keytype = "PROBEID", multiVals = 'first')
rownames(crisobmatsym)  = as.vector(genes)
crisobmatsym = crisobmatsym[-(which(is.na(rownames(crisobmatsym)))),]
crisobmatsym = collapseRows(crisobmatsym, rowGroup = unique(rownames(crisobmatsym)), rowID = unique(rownames(crisobmatsym)))
crisobmatsym = crisobmatsym$datETcollapsed

# all gene symbols were in the converted Cris' data:
length(crisgenesymraw) #271 genes
which(crisgenesymraw %in% rownames(crisobmatsym)) %>% length #271 genes
length(crisgenesymstd) #277 genes
which(crisgenesymstd %in% rownames(crisobmatsym)) %>% length #277 genes

# need to remove some cancer type specific
cancergenes = vector()
cancertypes = c('BLCA', 'CESC', 'COAD', 'KIRP', 'LIHC', 'READ', 'SKCM', 'UCEC')
for (i in 1:length(cancertypes)) {
	tmp = get(cancertypes[i])
	d = colnames(tmp)[which(!(colnames(tmp) %in% cancergenes))]
	cancergenes = d
}
# few genes were lost, when looking for genes common in both Cris' and ICGC data:
commongenesraw = crisgenesymraw[which(crisgenesymraw %in% colnames(BLCA))]
length(commongenesraw) # 221 obesity-associated genes common in ICGC data
commongenesstd = crisgenesymstd[which(crisgenesymstd %in% colnames(BLCA))]
length(commongenesstd) # 228 obesity-associated genes common in ICGC data

# check in Cris' BC data:
pdf(file='pdf/crisicgccommonmeta.pdf', width=7, height=7)
checkmat = crisobmatsym[commongenesraw, ]
checkmat = t(apply(checkmat, 1, function(x) (x-mean(x))/sd(x)))
tmpsvd = svd(checkmat)
tmpmeta = tmpsvd$v[,1]
tmpmeta = rank(tmpmeta)/length(tmpmeta)
checkmat[checkmat < -3] = -3
checkmat[checkmat > 3] = 3
metaplot3(checkmat, tmpmeta, crisobclin, name = "Common raw metagene from Cris' BC data in \nCris' BC data (raw)")
checkmat = crisobmatsym[commongenesstd, ]
checkmat = t(apply(checkmat, 1, function(x) (x-mean(x))/sd(x)))
tmpsvd = svd(checkmat)
tmpmeta = tmpsvd$v[,1]
tmpmeta = rank(tmpmeta)/length(tmpmeta)
checkmat[checkmat < -3] = -3
checkmat[checkmat > 3] = 3
metaplot3(checkmat, tmpmeta, crisobclin, name = "Common std. metagene from Cris' BC data in \nCris' BC data (raw)")
dev.off()

# Make transformation matrix in standardised Cris' data:
tmpmat = crisobmatsym[commongenesraw,]
tmpmat = t(apply(tmpmat, 1, function(x) (x-mean(x))/sd(x)))
tmpsvd = svd(tmpmat)
cristransmatraw = diag(1/tmpsvd$d) %*% t(tmpsvd$u)
tmpmat = crisobmatsym[commongenesstd,]
tmpmat = t(apply(tmpmat, 1, function(x) (x-mean(x))/sd(x)))
tmpsvd = svd(tmpmat)
cristransmatstd = diag(1/tmpsvd$d) %*% t(tmpsvd$u)

# Transform ICGC data:
pdf(file='pdf/crismetarawicgc.pdf', width=7, height=7)
for (i in 1:length(cancertypes)) {
	genes = commongenesraw
	bmi = get(paste(cancertypes[i], 'bmi', sep=''))
	mat = t(get(cancertypes[i]))
	mat = mat[genes,]
	mat = standardise_data(mat)
	tmpmeta = t(cristransmatraw %*% mat)
	tmpmeta = tmpmeta[,1]
	tmpmeta = rank(tmpmeta)/length(tmpmeta)
	mat[mat < -3] = -3
	mat[mat > 3] = 3
	txt = paste("Raw metagene from Cris' BC data in \n", cancertypes[i], sep='')
	metaplot3(mat, tmpmeta, bmi, name = txt)
}
dev.off()



# compare genes with Creighton's and my obesity-associated gene probes





# repeat, but with other variables taken into account (e.g. PR/ER/HER2)





# GE in different Auckland cohorts (?):































