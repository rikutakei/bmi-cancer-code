## All the functions for task13 goes in here:

## Function to change the ICGC data into TCGA data:
icgc_to_tcga = function(x) {
    tcgaid = gsub(pattern = '-[[:alnum:]]{3}-[[:alnum:]]{3}-[[:alnum:]]{4}-[[:alnum:]]{2}', replacement = '', rownames(x))
    rownames(x) = tcgaid
    return(x)
}

#Make a function to log and standardise cancer data
# make sure that the cancer data have genes as their columns, not rows
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
