###############################################################################
## miscallaneous stuff

###############################################################################
## Nick's potential BMI-related gene stuff

## Genes from TWAS/GWAS and/or eQTL studies that were related to BMI/obesity

## import genes:
nickbmigenes = read.csv(file = 'obsgenes/nickbmigene.csv')
nicknovelgenes = nickbmigenes[which(nickbmigenes$DataSource == 'novel eqtl'),]$GeneName
nicknovelgenes = as.vector(nicknovelgenes)
nickgenes = nickbmigenes$GeneName
nickgenes = as.vector(nickgenes)

## see how many genes match the cancer matrix:
which(nickgenes %in% colnames(BLCA)) %>% length ## 148 out of 176 genes in the cancer data
which(nicknovelgenes %in% colnames(BLCA)) %>% length ## 17 out of 20 novel genes in the cancer data

## extract the genes from the cancer data and normalise it:
cancertypes = c("BLCA", "CESC", "COAD", "KIRP", "LIHC", "READ", "SKCM", "UCEC")
for (i in 1:length(cancertypes)) {
	mat = get(cancertypes[i])
	genes = nickgenes[which(nickgenes %in% colnames(BLCA))]
	mat = t(mat[,genes])
	mat = standardise_data(mat, log=T)
	txt = paste('nick', cancertypes[i], sep='')
	assign(txt, mat)
}

## Make metagene and transformation matrix in BLCA data (training data set):
nicksvd = svd(nickBLCA)
nickmeta = nicksvd$v[,1]
nickmeta = rank(nickmeta)/length(nickmeta)
nick_transmat = diag(1/nicksvd$d) %*% t(nicksvd$u)

cancertypes = c("CESC", "COAD", "KIRP", "LIHC", "READ", "SKCM", "UCEC")
pdf(file='pdf/nickmeta.pdf', width=7, height=7)
metaplot3(nickBLCA, nickmeta, BLCAbmi, name = "Nick's metagene in BLCA")
for (i in 1:length(cancertypes)) {
	mat = get( paste('nick', cancertypes[i], sep=''))
	tmpmeta = t(nick_transmat %*% mat)
	tmpmeta = tmpmeta[,1]
	tmpmeta = rank(tmpmeta)/length(tmpmeta)
	bmi = get( paste(cancertypes[i], 'bmi', sep=''))
	txt = paste( "Nick's metagene in", cancertypes[i], sep=' ')
	txt = paste(txt, '(transformed)', sep=' ')
	metaplot3(mat, tmpmeta, bmi, name = txt)
}
dev.off()















