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

pdf('pdf/crmeta1.pdf')
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

pdf('pdf/crtcga.pdf')

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


