## functions_final.R file
##
## Add any (useful/worth keeping) functions I have used in my project in this file

## Function to change the ICGC data into TCGA data:
icgc_to_tcga = function(x) {
    tcgaid = gsub(pattern = '-[[:alnum:]]{3}-[[:alnum:]]{3}-[[:alnum:]]{4}-[[:alnum:]]{2}', replacement = '', rownames(x))
    rownames(x) = tcgaid
    return(x)
}

## Function that converts a gene probe ID matrix into gene symbol matrix
make_sym_mat = function (mat) {
	tmpgenes = mapIds(hgu133a.db, keys = rownames(mat), column = 'SYMBOL', keytype = "PROBEID", multiVals = 'first')
	ind = which(!is.na(tmpgenes))
	tmpgenes = tmpgenes[ind]
	mat = mat[ind,]
	mat = collapseRows(mat, rowID = rownames(mat), rowGroup = tmpgenes)
	mat = mat$datETcollapsed
	return(mat)
}

## Function that makes a metagene from a data
## mat is the matrix of interest, flip is to flip the metagene or not, and raw
## is to make a raw metagene rather than the ranked metagene
make_metagene = function(mat, flip = T, raw = T) {
	metagene = svd(mat)
	metagene = metagene$v[,1] # first principle component (i.e. the metagene)
	if (!raw) {
		metagene = rank(metagene)/length(metagene)
	}
	if (flip) {
		metagene= 1-metagene
	}
	return(metagene)
}

## Function to get the row dendrogram of a heatmap, given a matrix
get_dend = function(mat) {
	tmpmat = mat
	tmpmat[tmpmat > 3] = 3
	tmpmat[tmpmat < -3] = -3
	dend = heatmap.2(tmpmat, trace = 'none', scale = 'none')$rowDendrogram
	return(dend)
}

## Function to plot a heatmap, given a metagene.
mgheatmap = function(mat, metagene, dend = NULL, main = '') {
	tmpmat = mat
	tmpmat[tmpmat > 3] = 3
	tmpmat[tmpmat < -3] = -3
	metaord = order(rank(metagene))
	if (!is.null(dend )) {
		heatmap.2(tmpmat, trace='none', scale='none', col='bluered', Rowv = dend, main=main)
		heatmap.2(tmpmat, trace='none', scale='none', col='bluered', Rowv = dend, ColSideColors = bluered(length(metagene))[rank(metagene)], main=main)
		heatmap.2(tmpmat[,metaord], trace='none', scale='none', col='bluered', Rowv = dend, ColSideColors = bluered(length(metagene))[rank(metagene)][metaord], Colv=F, main=main)
	} else {
		heatmap.2(tmpmat, trace='none', scale='none', col='bluered', main=main)
		heatmap.2(tmpmat, trace='none', scale='none', col='bluered', ColSideColors = bluered(length(metagene))[rank(metagene)], main=main)
		heatmap.2(tmpmat[,metaord], trace='none', scale='none', col='bluered', ColSideColors = bluered(length(metagene))[rank(metagene)][metaord], Colv=F, main=main)
	}
}

## Function to plot boxplot and scatterplot of metagene vs. BMI/BMI status.
bmiplot = function(clin, metagene, main = '') {
	# boxplot
	bmifactor = factor(clin$bmiStatus, levels=c("normal","overweight", "obese"))
	boxplot(metagene~bmifactor, main=paste(main, "vs. BMI Status"), ylab = 'Metagene value', xlab = 'BMI Status', ylim = c(-0.05, 1.1))

	# p-value/legend for boxplot
	normind = which('normal' == bmifactor)
	ovind   = which('overweight' == bmifactor)
	obind   = which('obese' == bmifactor)
	txt = t.test(metagene[c(normind,ovind)]~bmifactor[c(normind,ovind)], alternative = 'two.sided')$p.value
	txt2 = t.test(metagene[c(normind,obind)]~bmifactor[c(normind,obind)], alternative = 'two.sided')$p.value
	txt3 = summary(aov(metagene~bmifactor))[[1]]$Pr[1]
	legend(x = 1.6, y = 1.08, bty='n', legend = as.expression(bquote(p == .(format(txt, digits=4)))))
	legend(x = 2.6, y = 1.08, bty='n', legend = as.expression(bquote(p == .(format(txt2, digits=4)))))
	legend('bottomright', bty='n', legend = as.expression(bquote("ANOVA p" == .(format(txt3, digits=4)))))

	# scatter plot
	plot(clin$bmi,metagene, pch = 20, main=paste(main, "vs. BMI"), ylab = 'Metagene value', xlab = 'BMI')

	# p-value/legend for scatter plot
	abline(lm(metagene~clin$bmi))
	txt = summary(lm(metagene~clin$bmi))$adj.r.squared
	txt2 = summary(lm(metagene~clin$bmi))$coef[2,4]
	legend('bottomright', bty='n', legend = c(as.expression(bquote(R^2 == .(format(txt, digits = 4)))), as.expression(bquote(p == .(format(txt2, digits=4))))))
}

## Function to make a transformation matrix from the given data
make_trans_mat = function(mat) {
	svd = svd(mat)
	transmat = diag(1/svd$d) %*% t(svd$u)
	return(transmat)
}

# Function to log and standardise RNA-seq data for heatmap visualisation.
# mat is the matrix to be standardised, and log is used if it is set to TRUE
standardise_data = function(mat, log = T) {
	if (log) {
		mat = log10(mat + 1) # add 1 before logging to avoid logging 0 values
	}
	mat = t(apply(mat, 1, function(mat) (mat-mean(mat))/sd(mat)))
	mat[is.nan(mat)] = 0
	return(mat)
}

# Function to make a metagene by applying a transformation matrix
# mat is the matrix to be standardised, transmat is the transformation matrix
# and the metagene is flipped if flip is set to TRUE
make_transmeta = function(mat, transmat, flip = T) {
	t = transmat %*% mat
	t = t[1,]
	t = rank(t)/length(t)
	if (flip) {
		t = 1-t
	}
	return(t)
}

# Function to get the top DEGs.
# mat is the gene expression matrix, group is a vector of the groups you want
# to do the DEG analysis against, and p is the p-value threshold for the top
# DEGs found.
# The analysis is FDR (Benjamini & Hochberg) corrected
get_deg = function(mat, group, p = 0.05) {
	design = model.matrix(~group)
	fit    = lmFit(mat, design)
	fit    = eBayes(fit)
	tt     = topTable(fit, coef = 2, adjust = 'BH', n=nrow(mat))
	tt     = tt[tt$P.Value <= p,]
	return(tt)
}

## Function to normalise RNA-seq data for gene expression analysis:
## Takes in a raw RNA-seq data (genes by samples) with a design matrix
voom_norm = function(x, group) {
	design = model.matrix(~group)
	a = x[rowSums(cpm(x)) > 9,]
	b = calcNormFactors(a)
	a = voom(a, design, lib.size = colSums(a)*b)
	return(a)
}

## Function to make all the Gatza metagenes
## mat is the data set you're making the metagenes in, metalist is the vector
## of the variable names of the gene signature, flip is a T/F vector that specifies
## whether to flip a certain metagene in the metalist, raw is to keep the
## metagene in the raw or rank-based form
make_multi_meta = function (mat, metalist, flip = flip, raw = T) {
	res = matrix(0, 0, ncol(mat))
	for (i in 1:length(metalist)) {
		genes = get(metalist[i])
		tmp = mat[genes,]
		tmpflip = flip[i]
		tmpmeta = make_metagene(tmp, flip = tmpflip, raw = raw)
		res = rbind(res, tmpmeta)
	}
	rownames(res) = metalist
	return(res)
}

## Function to simulate the genes that correlate with BMI
## Take in a gene by sample matrix and corresponding clinical data
sim_bmi_cor = function(mat, clin, n = 1000) {
	tmpbmi = sample(clin$bmi, length(clin$bmi))
	result = cor(t(mat), tmpbmi, method = "spearman")
	for (i in 2:n) {
		tmpbmi = sample(clin$bmi, length(clin$bmi))
		tmpres = cor(t(mat), tmpbmi, method = "spearman")
		result = cbind(result, tmpres)
		if ((i %% 50) == 0)
			print(paste( 'finished ', i))
	}
	return(result)
}

## Function to calculate p-values from a simulation.
## true_val is a vector of the true correlation value you got, and the mat
## contains the correlation values for each simulation where rows are the
## genes, and the columns are the simulation
get_pval = function(true_val, mat) {
	v = vector()
	for (i in 1:nrow(result)) {
		## probability of observing the true value by chance:
		tmp1 = sum(true_val[i] < mat[i,])/ncol(mat)
		tmp2 = sum(true_val[i] > mat[i,])/ncol(mat)
		val = min(tmp1, tmp2)
		v = c(v, val)
	}
	names(v) = rownames(mat)
	return(v)
}


## Function to return the pathways associated with the gene
## genename is the name of the gene you want to look up, and dbtype is either
## KEGG, GO, or reactome
lookup = function(genename, dbtype) {
	## check which database you're looking up in
	if (dbtype == "KEGG") {
		genetoentrez = KEGG.list
		entreztopath = keggpath
	} else if (dbtype == "GO") {
		genetoentrez = GO.list
		entreztopath = gopath
	} else if (dbtype == "reactome"){
		genetoentrez = reactome.list
		entreztopath = reactomepath
	} else {
		print('Unknown database type')
		return()
	}

	entrezID = 0
	entrezID = genetoentrez[[genename]] #get entrezID from the gene name

	#pull out the paths from the entrezID
	if(is.null(entrezID)) {
		return(entrezID)
	} else if (length(entrezID) > 1) {
		path = vector()
		for (i in 1:length(entrezID)) {
			path = c(path, entreztopath[[entrezID[i]]])
		}
	} else {
		path = entreztopath[[entrezID]]
	}

	return(path)
}
