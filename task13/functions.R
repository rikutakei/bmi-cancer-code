## All the functions for task13 goes in here:

## Function to change the ICGC data into TCGA data:
icgc_to_tcga = function(x) {
    tcgaid = gsub(pattern = '-[[:alnum:]]{3}-[[:alnum:]]{3}-[[:alnum:]]{4}-[[:alnum:]]{2}', replacement = '', rownames(x))
    rownames(x) = tcgaid
    return(x)
}

#Make a function to log and standardise cancer data
# make sure that the cancer data have genes as their rows, not columns
# if the matrix to be standardised doesn't have to be logged, then set the log variable to F
standardise_data = function(x, log = T) {
    if (log) {
        x = log10(x + 1)
    }
    x = t(apply(x, 1, function(x) (x-mean(x))/sd(x)))
    x[x < -3] = -3
    x[x > 3] = 3
	x[is.nan(x)] = 0
	#x = x[complete.cases(x),]
    return(x)
}

## Function to Voom normalise ICGC data:
## The RNA-seq data has to be genes by samples
voomICGC = function(x, bmi) {
	group = ifelse(bmi[,4] == 'obese', 'obese', 'non-obese')
	design = model.matrix(~group)
	a = x[rowSums(cpm(x)) > 9,]
	b = calcNormFactors(a)
	a = voom(x, design, lib.size = colSums(a)*b)
	return(a)
}

#make_tt function from task09:
# Function to make topTable from RNA-seq and model
make_tt = function(x, model) {
    fit = lmFit(x, model)
    fit = eBayes(fit)
    tt = topTable(fit, coef = 2, adjust = 'BH', n = nrow(x))
    return(tt)
}

# Function to pull out the DEGs from the top table
pull_deg = function(x, adj = F, y = 0.05) {
    if(adj) {
        x = x[x$adj.P.Val <= y,]
    } else {
        x = x[x$P.Value <= y,]
    }
    return(x)
}

## Function to print heatmaps and other plots:
metaplot = function(x, meta, bmi, dendrogram, name = '') {
	ord = order(meta)
	heatmap.2(x, trace='none',scale='none', col='bluered', main=name)
	heatmap.2(x, trace='none',scale='none', col='bluered', ColSideColors = bluered(length(meta))[rank(meta)], Colv=NA, Rowv = as.dendrogram(dendrogram), dendrogram='row', main=name)
	heatmap.2(x[,ord], trace='none',scale='none', col='bluered', ColSideColors = bluered(length(meta))[rank(meta)][ord], Colv=NA, Rowv = as.dendrogram(dendrogram), dendrogram='row', main=name)
	bmifactor = factor(bmi$bmiStatus, levels=c("normal","overweight", "obese"))
	boxplot(meta~bmifactor, main=paste(name, "vs. BMI Status"), ylab = 'Metagene Score', xlab = 'BMI Status')
	txt = summary(aov(meta~bmifactor))[[1]]$Pr[1]
	legend('top', bty='n', legend = as.expression(bquote(p == .(format(txt, digits=4)))))
	plot(bmi$bmi,meta, pch = 20, main=paste(name, "vs. BMI"), ylab = 'Metagene Score', xlab = 'BMI')
	abline(lm(meta~bmi$bmi))
	txt = summary(lm(meta~bmi$bmi))$adj.r.squared
	txt2 = summary(lm(meta~bmi$bmi))$coef[2,4]
	legend('bottomright', bty='n', legend = c(as.expression(bquote(R^2 == .(format(txt, digits = 4)))), as.expression(bquote(p == .(format(txt2, digits=4))))))
}

metaplot2 = function(x, meta, bmi, name = '') {
	ord = order(meta)
	heatmap.2(x, trace='none',scale='none', col='bluered', main=name)
	heatmap.2(x, trace='none',scale='none', col='bluered', ColSideColors = bluered(length(meta))[rank(meta)], Colv=NA, main=name)
	heatmap.2(x[,ord], trace='none',scale='none', col='bluered', ColSideColors = bluered(length(meta))[rank(meta)][ord], Colv=NA, main=name)
	bmifactor = factor(bmi$bmiStatus, levels=c("normal","overweight", "obese"))
	boxplot(meta~bmifactor, main=paste(name, "vs. BMI Status"), ylab = 'Metagene Score', xlab = 'BMI Status')
	txt = summary(aov(meta~bmifactor))[[1]]$Pr[1]
	legend('top', bty='n', legend = as.expression(bquote(p == .(format(txt, digits=4)))))
	plot(bmi$bmi,meta, pch = 20, main=paste(name, "vs. BMI"), ylab = 'Metagene Score', xlab = 'BMI')
	abline(lm(meta~bmi$bmi))
	txt = summary(lm(meta~bmi$bmi))$adj.r.squared
	txt2 = summary(lm(meta~bmi$bmi))$coef[2,4]
	legend('bottomright', bty='n', legend = c(as.expression(bquote(R^2 == .(format(txt, digits = 4)))), as.expression(bquote(p == .(format(txt2, digits=4))))))
}

## Function to print heatmaps and other plots with p-values:
metaplot3 = function(x, meta, bmi, name = '') {
	ord = order(meta)
	#heatmap
	heatmap.2(x, trace='none',scale='none', col='bluered', main=name)
	heatmap.2(x, trace='none',scale='none', col='bluered', ColSideColors = bluered(length(meta))[rank(meta)], Colv=NA, main=name)
	heatmap.2(x[,ord], trace='none',scale='none', col='bluered', ColSideColors = bluered(length(meta))[rank(meta)][ord], Colv=NA, main=name)

	# boxplot
	bmifactor = factor(bmi$bmiStatus, levels=c("normal","overweight", "obese"))
	boxplot(meta~bmifactor, main=paste(name, "vs. BMI Status"), ylab = 'Metagene Score', xlab = 'BMI Status', ylim = c(-0.05, 1.1))

	#p-value/legend for boxplot
	normind = which('normal' == bmifactor)
	ovind = which('overweight' == bmifactor)
	obind = which('obese' == bmifactor)
	txt = t.test(meta[c(normind,ovind)]~bmifactor[c(normind,ovind)], alternative= 'two.sided')$p.value
	txt2 = t.test(meta[c(normind,obind)]~bmifactor[c(normind,obind)], alternative= 'two.sided')$p.value
	txt3 = summary(aov(meta~bmifactor))[[1]]$Pr[1]
	legend(x = 1.6, y = 1.08, bty='n', legend = as.expression(bquote(p == .(format(txt, digits=4)))))
	legend(x = 2.6, y = 1.08, bty='n', legend = as.expression(bquote(p == .(format(txt2, digits=4)))))
	legend('bottomright', bty='n', legend = as.expression(bquote("ANOVA p" == .(format(txt3, digits=4)))))

	# scatter plot
	plot(bmi$bmi,meta, pch = 20, main=paste(name, "vs. BMI"), ylab = 'Metagene Score', xlab = 'BMI')
	abline(lm(meta~bmi$bmi))
	txt = summary(lm(meta~bmi$bmi))$adj.r.squared
	txt2 = summary(lm(meta~bmi$bmi))$coef[2,4]
	legend('bottomright', bty='n', legend = c(as.expression(bquote(R^2 == .(format(txt, digits = 4)))), as.expression(bquote(p == .(format(txt2, digits=4))))))
}

# Function to check the given metagene in a list of data.
# You need to give a list of the variable names of the data you want to check
# your metagene in.
# genelist should contain the list of genes used to create the transformation
# matrix, and the transmat is the transformation matrix.
# If you want to flip the metagene, then set the flip to TRUE
# This function is designed for a list of ICGC data, but I guess it can work
# for other cancer data as well
check_data = function(datlist, bmilist, genelist, transmat, log=T, flip=F, name='') {
	for (i in 1:length(datlist)) {
		testdat = t(get(datlist[i]))
		mat = testdat[genelist,]
		mat = standardise_data(mat, log=log)
		bmi = get(bmilist[i])
		tmpsvd = t(transmat %*% mat)
		tmpsvd = tmpsvd[,1]
		tmpsvd = rank(tmpsvd)/length(tmpsvd)
		if(flip) {
			tmpsvd = 1-tmpsvd
		}
		main = paste(name, '(')
		main = paste(main, datlist[i], sep = '')
		main = paste(main, ')', sep = '')
		metaplot3(mat, tmpsvd, bmi, name = main)
	}
}

# Function to make a set of colour for heatmap.2x
make_col = function(values, continuous = T, colours = c()) {
	if (continuous) {
		col = bluered(length(values))[rank(values)]
		return(col)
	} else {
		u = sort(unique(values))
		if (length(u) > 8) {
			print('There were too many values')
			return()
		}
		else if (length(colours) > 0){
			col = values
			brewercol = colours
			for (i in 1:length(col)) {
				for(j in 1:length(u)) {
					if(is.na(col[i])){
						col[i] = "#FFFFFF"
						break()
					}
					if (col[i] == u[j]) {
						col[i] = brewercol[j]
						break()
					}
				}
			}
		}
		else {
			col = values
			brewercol = brewer.pal(length(u), "Set1")
			for (i in 1:length(col)) {
				for(j in 1:length(u)) {
					if(is.na(col[i])){
						col[i] = "#FFFFFF"
						break()
					}
					if (col[i] == u[j]) {
						col[i] = brewercol[j]
						break()
					}
				}
			}
		}
		return(col)
	}
}

# Function to plot metagenes made from raw and standardised data:
plot_raw_vs_std <- function (dat, pathways, main='') {
	for (i in 1:length(pathways)) {
		gene = get(pathways[i])
		mat = dat[gene,]
		rawsvd = svd(mat)
		rawmeta = rank(rawsvd$v[,1])/ncol(mat)

		stdmat = t(apply(mat, 1, function(x) (x-mean(x))/sd(x)))
		stdsvd = svd(stdmat)
		stdmeta = rank(stdsvd$v[,1])/ncol(stdmat)

		txt = paste(main, '(')
		txt = paste(txt, pathways[i], sep='')
		txt = paste(txt, ')', sep='')
		plot(stdmeta, rawmeta, main=txt)
	}
}

# Function to make a heatmap of metagenes
metamap <- function (metalist, col, main='') {
	
	heatmap.2x(metalist, trace='none',scale='none', col='bluered', ColSideColors = col, Colv=NA, main=main, cexRow=1.0)
	heatmap.2x(metalist, trace='none',scale='none', col=matlab.like(ncol(metalist)), ColSideColors = col, Colv=NA, main=main, cexRow=1.0)

	c = cor(t(metalist), method = 'spearman')
	c2 = cor(t(metalist), method = 'pearson')

	txt = paste("Spearman correlation of\n", main, sep='')
	heatmap.2x(c, trace='none',scale='none', col='bluered', main=txt, cexRow=1.0)
	heatmap.2x(c, trace='none',scale='none', col=matlab.like(ncol(metalist)),  main=txt, cexRow=1.0)

	txt = paste("Pearson correlation of\n", main, sep='')
	heatmap.2x(c2, trace='none',scale='none', col='bluered', main=txt, cexRow=1.0)
	heatmap.2x(c2, trace='none',scale='none', col=matlab.like(ncol(metalist)),  main=txt, cexRow=1.0)

}

# Function to do gatza pathway stuff:
# mat = gene expression matrix
# matheat = gene expression matrix that is used ONLY for the heatmap
# pathlist = vector of the variable names of the gatza pathways
# flip = vector of the pathways to be flipped
# metalist = name of the variable you want to store the metagene values
# corlist = name of the variable you want to store the correlation values
# checkgene = vector of genes that represent the expression of the specific gatza pathway
# rank = boolean; whether to use probit or rank-based metagene scaling
gatzaPath <- function (mat, matheat, pathlist, flip, metalist, corlist, rank = T, checkgene) {
	tmpmetalist = list()
	tmpcorlist = list()
	for (i in 1:length(pathlist)) {
		gene = get(pathlist[i])
		tmpmat = mat[gene,]

		tmpsvd = svd(tmpmat)
		tmpmeta = tmpsvd$v[,1]

		# flip metagene:
		if (pathlist[i] %in% flip) {
			tmpmeta = 1-tmpmeta
		}
		
		# rank-based, or probit scaling:
		if (rank) {
			tmpmeta = rank(tmpmeta)/length(tmpmeta)
		} else {
			tmpmeta = pnorm(scale(tmpmeta))
		}

		ord = order(tmpmeta)
		col = bluered(length(tmpmeta))[rank(tmpmeta)]
		col2 = bluered(length(tmpmeta))[rank(matheat[checkgene[i],])]
		col = rbind(col2, Metagene=col)
		rownames(col)[1] = checkgene[i]
		heatmap.2x(matheat[,ord], trace='none',scale='none', col='bluered', ColSideColors = col[,ord], Colv=NA, main=pathlist[i], cexRow=1.0)
		tmpmetalist[[i]] = tmpmeta
		cor = cor(tmpmeta, matheat[i,])
		tmpcorlist[[i]] = cor
	}

	tmpmetalist = as.data.frame(tmpmetalist)
	colnames(tmpmetalist) = gsub('_probes', '', pathlist)
	tmpmetalist = t(as.matrix(tmpmetalist))

	names(tmpcorlist) = pathlist
	tmpcorlist = as.data.frame(tmpcorlist)

	tmpord = c('er', 'pr', 'p53', 'bcat', 'e2f1', 'pi3k', 'myc', 'ras', 'ifna', 'ifng', 'akt', 'p63', 'src', 'her2', 'egfr', 'tgfb', 'stat3', 'tnfa')

	heatmap.2(tmpmetalist, trace='none',scale='none', col=matlab.like, cexRow=1, main='Gatza metagenes')
	heatmap.2(tmpmetalist[tmpord,], trace='none',scale='none', col=matlab.like, cexRow=1, Rowv=F, main='Gatza metagenes ordered as in their paper')

	cor = cor(t(tmpmetalist), method='pearson')
	heatmap.2(cor, trace='none',scale='none', col=matlab.like, cexRow=1)
	heatmap.2(cor[tmpord, tmpord], trace='none',scale='none', col=matlab.like, cexRow=1, Rowv=F, Colv=F)
	return(list(metagene = tmpmetalist, correlation = tmpcorlist))
}

# Function that applies the Gatza transformation matrices to a dataset and plots a heatmap of the metagenes.
# mat = dataset for the transformation matrix to be applied
# metalist = variable names of the pathways used to create the transformation matrix (length should be the same as the length of translist)
# translist = list containing all of the transformation matrix, in the order of the metalist names
# flip = a vector of pathway names that describes which metagene should be flipped
# main = name of the heatmap
gttransfun <- function (mat, metalist, translist, flip, main='') {
	tmpmetalist = list()
	for (i in 1:length(metalist)) {
		# prepare the dataset for transformation:
		gene = get(metalist[i])
		tmpmat = mat[gene,]

		# transform the data:
		tmpmat = translist[[i]] %*% tmpmat
		tmpmeta = tmpmat[1,]

		# rank the metagenes:
		tmpmeta = rank(tmpmeta)/length(tmpmeta)

		# flip metagenes that are in the list:
		if (metalist[i] %in% flip) {
			tmpmeta = 1-tmpmeta
		}

		# save the metagene:
		tmpmetalist[[i]] = tmpmeta
		names(tmpmetalist)[i] = metalist[i]
	}
	# return(tmpmetalist)
	tmp = as.data.frame(tmpmetalist)
	tmp= t(as.matrix(tmp))
	tmpcor = cor(t(tmp), method = 'pearson')

	# make some heatmaps:
	heatmap.2(x=tmp, scale='none', trace='none', col='bluered', main=main)
	heatmap.2(x=tmp, scale='none', trace='none', col=matlab.like, main=main)
	heatmap.2(x=tmpcor, scale='none', trace='none', col='bluered', main=main)
	heatmap.2(x=tmpcor, scale='none', trace='none', col=matlab.like, main=main)

	reslist = list(metagene=tmp, correlation=tmpcor)
	return(reslist)
}

# Function to get the metagene value from both transformation matrix and svd by itself
# mat = data matrix to get the metagene from
# metalist = variable names of the pathways
# translist = list that contains the transformation matrix
getmeta <- function (mat, metalist, translist) {
	corval = list()
	svdmeta = list()
	transmeta = list()
	for (i in 1:length(metalist)) {
		gene = get(metalist[i])
		tmp = mat[gene,]

		# make metagene from svd:
		tmpsvd = svd(tmp)
		tmpmeta = tmpsvd$v[,1]
		tmpmeta = rank(tmpmeta)/length(tmpmeta)

		# make metagene from transformation matrix:
		trans = translist[[i]]
		metatrans = trans %*% tmp
		metatrans = metatrans[1,]
		metatrans = rank(metatrans)/length(metatrans)

		svdmeta[[i]] = tmpmeta
		names(svdmeta)[i] = metalist[i]
		transmeta[[i]] = metatrans
		names(svdmeta)[i] = metalist[i]
		tmpcor = cor(tmpmeta, metatrans, method = 'spearman')
		corval[[i]] = tmpcor
	}
	corval = unlist(corval)
	names(corval) = metalist
	# return(corval)
	svdmeta = as.data.frame(svdmeta)
	transmeta = as.data.frame(transmeta)
	return(list(svd = svdmeta, trans = transmeta, correlation = corval))
}


















