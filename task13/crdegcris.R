###############################################################################
## BMI metagene from Creighton data in Cris' data

## Remember to rerun the DEG analysis to get the raw (gene probe ID, not gene symbol) BMI gene probe sets
## (See creightonGE.R)

# see if there are any probe ID not in Cris' data
num = c()
num2 = c()
for (i in 1:length(allobsname)) {
	ind = which(get(allobsname[i]) %in% rownames(cris_raw))
	ind2 = length(get(allobsname[i]))
	num = c(num, length(ind))
	num2 = c(num2, ind2)
}

## all probes are in Cris' data

## start making transformation matrix:

for (i in 1:length(allobsname)) {
	gene = get(allobsname[i])
	mat = cr_raw[gene,]
	tmpsvd = svd(mat)
	mat = t(apply(mat, 1, function(x) (x-mean(x))/sd(x)))
	tmpsvd2 = svd(mat)

	txt = gsub("genes", "rawtransmat", allobsname[i])
	trmat = diag(1/tmpsvd$d) %*% t(tmpsvd$u)
	assign(txt, trmat)

	txt = gsub("genes", "stdtransmat", allobsname[i])
	trmat = diag(1/tmpsvd2$d) %*% t(tmpsvd2$u)
	assign(txt, trmat)
}

## Transform Cris' data using transformation matrix and get metagene:
rawtransmat = gsub("genes", "rawtransmat", allobsname)
stdtransmat = gsub("genes", "stdtransmat", allobsname)
bmimetaraw = list()
bmimetastd = list()
for (i in 1:length(allobsname)) {
	gene = get(allobsname[i])
	mat = cris_raw[gene,]
	trmat = get(rawtransmat[i])
	tmpmeta = t(trmat %*% mat)
	tmpmeta = tmpmeta[,1]
	tmpmeta = rank(tmpmeta)/length(tmpmeta)
	bmimetaraw[[i]] = tmpmeta

	mat = t(apply(mat, 1, function(x) (x-mean(x))/sd(x)))
	trmat = get(stdtransmat[i])
	tmpmeta = t(trmat %*% mat)
	tmpmeta = tmpmeta[,1]
	tmpmeta = rank(tmpmeta)/length(tmpmeta)
	bmimetastd[[i]] = tmpmeta
}
bmimetaraw = t(as.data.frame(bmimetaraw))
rownames(bmimetaraw) = allobsname
bmimetastd = t(as.data.frame(bmimetastd))
rownames(bmimetastd) = allobsname

## Check for the direction of the metagenes:

# first, find the genes that are common across all of the metagenes:
checkgenes= c(rawobsgenes,crolgenes,resobsgenes,rescrolgenes, caobsgenes,cacrolgenes,caresobsgenes,carescrolgenes)
checkgenes = table(checkgenes)[which(table(checkgenes) == 8)]
checkgenes = names(checkgenes)

# make row colours for heatmap
tmpcol = c('blue', 'red', 'white')
ER = make_col(as.vector(crisclin$ER.status), continuous=F, colours = tmpcol)
PR = make_col(as.vector(crisclin$PgR.status), continuous=F, colours = tmpcol)
LN = make_col(as.vector(crisclin$LN.status), continuous=F, colours = tmpcol)
tmpcol = brewer.pal(6, 'Dark2')
SubType = make_col(as.vector(crisclin$subtype), continuous=F, colours = tmpcol)
tmpcol = brewer.pal(3, 'YlOrRd')
Grade = make_col(as.vector(crisclin$Grade), continuous=F, colours = tmpcol)
BMIStatus = as.vector(crisclin$bmiStatus)
BMIStatus[which(BMIStatus == "NA")] = "#FFFFFF"
BMIStatus[which(BMIStatus == "normal")] = tmpcol[1]
BMIStatus[which(BMIStatus == "overweight")] = tmpcol[2]
BMIStatus[which(BMIStatus == "obese")] = tmpcol[3]


col = rbind(ER = ER, PR = PR, LN = LN, BMIStatus = BMIStatus, Grade = Grade, SubType = SubType)

namelist = c("Raw metagene","CR overlap metagene (raw)","Residual metagene","CR overlap metagene (residual)","Raw metagene (Caucasian/raw)","CR overlap metagene (Caucasian/raw)","Residual metagene (Caucasian/residual)","CR overlap metagene (Caucasian/residual)")

# make a matrix for the heatmap:
matheat = cris_raw[checkgenes,]
matheat = t(apply(matheat, 1, function(x) (x-mean(x))/sd(x)))
matheat[matheat < -3] = -3
matheat[matheat > 3] = 3

# transform raw data, using raw transformation matrix:

pdf(file='pdf/bmimetacrisraw.pdf', width=7, height=7)
for (i in 1:nrow(bmimetaraw)) {
	tmpmeta = bmimetaraw[i,]
	if (i == 7) {
		tmpmeta = 1-bmimetaraw[i,]
	}
	ord = order(tmpmeta)
	tmpcol = rbind(col, meta = bluered(length(tmpmeta))[rank(tmpmeta)])
	txt = namelist[i]
	txt = paste(txt, "\nin Cris' raw data (transformed)")
	heatmap.2x(matheat[,ord], trace='none',scale='none', col='bluered', ColSideColors = tmpcol[,ord], Colv=NA, main=txt, cexRow=1.0)
}
dev.off()

# transform standardised data, using standardised transformation matrix:
pdf(file='pdf/bmimetacrisstd.pdf', width=7, height=7)
for (i in 1:nrow(bmimetastd)) {
	tmpmeta = bmimetastd[i,]
	if (i == 7) {
		tmpmeta = 1-bmimetastd[i,]
	}
	ord = order(tmpmeta)
	tmpcol = rbind(col, meta = bluered(length(tmpmeta))[rank(tmpmeta)])
	txt = namelist[i]
	txt = paste(txt, "\nin Cris' standardised data (transformed)")
	heatmap.2x(matheat[,ord], trace='none',scale='none', col='bluered', ColSideColors = tmpcol[,ord], Colv=NA, main=txt, cexRow=1.0)
}
dev.off()

# In both cases, need to flip the 7th metagene:
bmimetaraw[7,] = 1-bmimetaraw[7,]
bmimetastd[7,] = 1-bmimetastd[7,]

## Just a quick look at raw vs std metagene values (transformation matrix not used)
pdf(file='pdf/crrawvsstdcris.pdf', width=7, height=7)
plot_raw_vs_std(cris_raw, allobsname, main='Cris BC data')
dev.off()

###############################################################################
## Just a quick check before looking for clustering with gatza metagenes, to see if any of the BMI metagenes correlate with sample BMI in Cris' data

## get data with samples that contain BMI data:
crisobmat = cris_raw[,crisobclin$CelFile]

pdf(file='pdf/bmimetacrisrawclin.pdf', width=7, height=7)
for (i in 1:length(allobsname)) {
	gene = get(allobsname[i])
	mat = crisobmat[gene,]
	trmat = get(rawtransmat[i])
	tmpmeta = t(trmat %*% mat)
	tmpmeta = tmpmeta[,1]
	tmpmeta = rank(tmpmeta)/length(tmpmeta)
	if(i == 7) {
		tmpmeta = 1-tmpmeta
	}
	txt = namelist[i]
	txt = paste(txt, "\nin Cris' standardised data (transformed)")
	matheat = t(apply(mat, 1, function(x) (x-mean(x))/sd(x)))
	matheat[matheat < -3] = -3
	matheat[matheat > 3] = 3
	metaplot3(matheat, tmpmeta, crisobclin, name = txt)
}
dev.off()

pdf(file='pdf/bmimetacrisstdclin.pdf', width=7, height=7)
for (i in 1:length(allobsname)) {
	gene = get(allobsname[i])
	mat = crisobmat[gene,]
	mat = t(apply(mat, 1, function(x) (x-mean(x))/sd(x)))
	trmat = get(stdtransmat[i])
	tmpmeta = t(trmat %*% mat)
	tmpmeta = tmpmeta[,1]
	tmpmeta = rank(tmpmeta)/length(tmpmeta)
	if(i == 7) {
		tmpmeta = 1-tmpmeta
	}
	txt = namelist[i]
	txt = paste(txt, "\nin Cris' standardised data (transformed)")
	matheat = t(apply(mat, 1, function(x) (x-mean(x))/sd(x)))
	matheat[matheat < -3] = -3
	matheat[matheat > 3] = 3
	metaplot3(matheat, tmpmeta, crisobclin, name = txt)
}
dev.off()

###############################################################################
## See which Gatza pathway metagene the BMI metagenes cluster together with (in Cris' data)

crisallmetalist = rbind(crisgatzametalist, bmimetastd)
rownames(crisallmetalist) = gsub('genes', '', rownames(crisallmetalist))

pdf(file='pdf/allmetacris.pdf', width=7, height=7)
heatmap.2x(crisallmetalist, trace='none',scale='none', col='bluered', ColSideColors = col, Colv=NA, main="All metagenes in Cris' BC data", cexRow=1.0)
heatmap.2x(crisallmetalist, trace='none',scale='none', col=matlab.like(ncol(crisallmetalist)), ColSideColors = col, Colv=NA, main="All metagenes in Cris' BC data", cexRow=1.0)
c = cor(t(crisallmetalist), method = 'spearman')
c2 = cor(t(crisallmetalist), method = 'pearson')
heatmap.2x(c, trace='none',scale='none', col='bluered', main="Spearman correlation of all\n metagenes in Cris' BC data", cexRow=1.0)
heatmap.2x(c, trace='none',scale='none', col=matlab.like(ncol(crisallmetalist)),  main="Spearman correlation of all \nmetagenes in Cris' BC data", cexRow=1.0)
heatmap.2x(c2, trace='none',scale='none', col='bluered', main="Pearson correlation of all\n metagenes in Cris' BC data", cexRow=1.0)
heatmap.2x(c2, trace='none',scale='none', col=matlab.like(ncol(crisallmetalist)),  main="Pearson correlation of all \nmetagenes in Cris' BC data", cexRow=1.0)

## BMI metagenes by themselves:
heatmap.2x(bmimetastd, trace='none',scale='none', col='bluered', ColSideColors = col, Colv=NA, main="BMI metagenes in Cris' BC data", cexRow=1.0)
heatmap.2x(bmimetastd, trace='none',scale='none', col=matlab.like(ncol(bmimetastd)), ColSideColors = col, Colv=NA, main="BMI metagenes in Cris' BC data", cexRow=1.0)
c = cor(t(bmimetastd), method = 'spearman')
c2 = cor(t(bmimetastd), method = 'pearson')
heatmap.2x(c, trace='none',scale='none', col='bluered', main="Spearman correlation of BMI\n metagenes in Cris' BC data", cexRow=1.0)
heatmap.2x(c, trace='none',scale='none', col=matlab.like(ncol(crisallmetalist)),  main="Spearman correlation of BMI \nmetagenes in Cris' BC data", cexRow=1.0)
heatmap.2x(c2, trace='none',scale='none', col='bluered', main="Pearson correlation of BMI\n metagenes in Cris' BC data", cexRow=1.0)
heatmap.2x(c2, trace='none',scale='none', col=matlab.like(ncol(crisallmetalist)),  main="Pearson correlation of BMI \nmetagenes in Cris' BC data", cexRow=1.0)

## Gatza metagenes by themselves:
heatmap.2x(crisgatzametalist, trace='none',scale='none', col='bluered', ColSideColors = col, Colv=NA, main="Gatza metagenes in Cris' BC data", cexRow=1.0)
heatmap.2x(crisgatzametalist, trace='none',scale='none', col=matlab.like(ncol(crisgatzametalist)), ColSideColors = col, Colv=NA, main="Gatza metagenes in Cris' BC data", cexRow=1.0)
c = cor(t(crisgatzametalist), method = 'spearman')
c2 = cor(t(crisgatzametalist), method = 'pearson')
heatmap.2x(c, trace='none',scale='none', col='bluered', main="Spearman correlation of Gatza\n metagenes in Cris' BC data", cexRow=1.0)
heatmap.2x(c, trace='none',scale='none', col=matlab.like(ncol(crisallmetalist)),  main="Spearman correlation of Gatza \nmetagenes in Cris' BC data", cexRow=1.0)
heatmap.2x(c2, trace='none',scale='none', col='bluered', main="Pearson correlation of Gatza\n metagenes in Cris' BC data", cexRow=1.0)
heatmap.2x(c2, trace='none',scale='none', col=matlab.like(ncol(crisallmetalist)),  main="Pearson correlation of Gatza \nmetagenes in Cris' BC data", cexRow=1.0)
dev.off()

###############################################################################
# see if get any difference when I split the samples into Auckland A and Auckland B cohort groups

aukA = crisobclin[which(crisobclin$Auck.cohort == "Auckland A"),]
aukA = paste(aukA$GeoFile, aukA$CelFile, sep='_') # 30 patients
aukB = crisobclin[which(crisobclin$Auck.cohort == "Auckland B"),]
aukB = paste(aukB$GeoFile, aukB$CelFile, sep='_') # 77 patients

pdf(file='pdf/auckAmetamap.pdf', width=7, height=7)

aukAmetalist = crisallmetalist[,aukA]
aukAcol = col[,78:107] # Auckland A cohort is the last 30 patients in the data

metamap(metalist = aukAmetalist, col = aukAcol, main = "All metagenes in Cris' BC data (AuckA)")

## BMI metagenes by themselves:
metamap(metalist = aukAmetalist[19:26,], col = aukAcol, main = "BMI metagenes in Cris' BC data (AuckA)")

## Gatza metagenes by themselves:
metamap(metalist = aukAmetalist[1:18,], col = aukAcol, main = "Gatza metagenes in Cris' BC data (AuckA)")

dev.off()

pdf(file='pdf/auckBmetamap.pdf', width=7, height=7)

aukBmetalist = crisallmetalist[,aukB]
aukBcol = col[,1:77] # Auckland A cohort is the last 30 patients in the data

metamap(metalist = aukBmetalist, col = aukBcol, main = "All metagenes in Cris' BC data (AuckB)")

## BMI metagenes by themselves:
metamap(metalist = aukBmetalist[19:26,], col = aukBcol, main = "BMI metagenes in Cris' BC data (AuckB)")

## Gatza metagenes by themselves:
metamap(metalist = aukBmetalist[1:18,], col = aukBcol, main = "Gatza metagenes in Cris' BC data (AuckB)")

dev.off()



















