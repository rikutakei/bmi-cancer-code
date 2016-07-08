###############################################################################
## Creighton metagene in Cris' data

# extract Creighton's obesity-associated genes from Cris' data:
crisobmat = cris_raw[cr_obsgene,]

# need to find which samples have BMI information:
crisobclin = crisclin

# add BMIStatus variable to the data:

BMIStatus = crisobclin$BMI
for (i in 1:length(BMIStatus)) {
	if (is.na(BMIStatus[i])){
		BMIStatus[i] = 'NA'
	} else {
		if (BMIStatus[i] >= 30) {
			BMIStatus[i] = 'obese'
		} else if (BMIStatus[i] < 25) {
			BMIStatus[i] = 'normal'
		} else {
			BMIStatus[i] = 'overweight'
		}
	}
}
crisobclin = cbind(crisobclin, bmiStatus = BMIStatus)
names(crisobclin)[26] = 'bmi'

# adjust the colnames to match the clinical data:
n = colnames(cris_raw)
n = gsub('GSM[[:digit:]_]*', '', n)
colnames(crisobmat) = n

crisobmat = crisobmat[,crisobclin$CelFile]

###############################################################################
## make transformation matrix in Creighton data, then use it to make metagene in Cris' data

cr_obsmat2 = cr_raw[cr_obsgene,]
cr_obsmat2_std = t(apply(cr_obsmat2, 1, function(x) (x-mean(x))/sd(x)))

crsvdraw = svd(cr_obsmat2)
crtransmatraw = diag(1/crsvdraw$d) %*% t(crsvdraw$u)

crsvdstd = svd(cr_obsmat2_std)
crtransmatstd = diag(1/crsvdstd$d) %*% t(crsvdstd$u)

pdf(file='pdf/criscrmeta.pdf', width=7, height=7)
trans = c('crtransmatraw', 'crtransmatstd')
for (i in 1:length(trans)) {
	trmat = get(trans[i])
	tmp = t(trmat %*% crisobmat)
	meta = tmp[,1]
	meta = rank(meta)/length(meta)
	if (i == 1) {
		meta = 1-meta
	}
	mat = crisobmat
	mat = t(apply(crisobmat, 1, function(x) (x-mean(x))/sd(x)))
	mat[mat > 3] = 3
	mat[mat < -3] = -3
	metaplot3(mat, meta, crisobclin, name="Cris Print's BC Data")
}
dev.off()

