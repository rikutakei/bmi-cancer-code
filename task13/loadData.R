###############################################################################
## Creighton et al data:

files = readLines('./raw/creighton/files.txt')
files = paste('./raw/creighton/', files, sep='')
cr_raw = ReadAffy(filenames = files)
cr_raw = rma(cr_raw) ## RMA normalise the data
cr_raw = exprs(cr_raw) ## change the format into matrix

cr_obsgene = read.csv('./obsgenes/crobsgenes.txt', header=F)
cr_obsgene = as.vector(cr_obsgene[,1])

crclin = read.csv('./clindata/crclin.csv', sep=',', header=T)

cr_symmat = cr_raw
tmpgenes = mapIds(hgu133a.db, keys = rownames(cr_symmat), column = 'SYMBOL', keytype = "PROBEID", multiVals = 'first')
rownames(cr_symmat) = tmpgenes
cr_symmat = cr_symmat[-(which(is.na(rownames(cr_symmat)))),]
cr_symmat = collapseRows(cr_symmat, unique(rownames(cr_symmat)), unique(rownames(cr_symmat)))
cr_symmat = cr_symmat$datETcollapsed
dim(cr_symmat) #13031 genes

###############################################################################
## Fuentes-Mattei et al data:

files = readLines('./raw/fuentes-mattei/files.txt')
files = paste('./raw/fuentes-mattei/', files, sep='')
fm_raw = ReadAffy(filenames = files)
fm_raw = rma(fm_raw) ## RMA normalise the data
fm_raw = exprs(fm_raw) ## change the format into matrix

fm_obsgene = read.csv('./obsgenes/fmobsgenes.txt', header=T)

fmclin = read.csv('./clindata/fmclin.csv', sep=',', header=T)

fm_symmat = fm_raw
tmpgenes = mapIds(hgu133a.db, keys = rownames(fm_symmat), column = 'SYMBOL', keytype = "PROBEID", multiVals = 'first')
rownames(fm_symmat) = tmpgenes
fm_symmat = fm_symmat[-(which(is.na(rownames(fm_symmat)))),]
fm_symmat = collapseRows(fm_symmat, unique(rownames(fm_symmat)), unique(rownames(fm_symmat)))
fm_symmat = fm_symmat$datETcollapsed
dim(fm_symmat) #13031 genes

###############################################################################
## Cris Print's Breast cancer data:

files = readLines('./raw/cris/files.txt')
files = paste('./raw/cris/', files, sep='')
cris_raw = ReadAffy(filenames = files)
cris_raw = rma(cris_raw) ## RMA normalise the data
cris_raw = exprs(cris_raw) ## change the format into matrix

crisclin = read.csv('./clindata/crisclin.csv', sep=',', header=T)

cris_symmat = cris_raw
tmpgenes = mapIds(hgu133plus2.db, keys = rownames(cris_symmat), column = 'SYMBOL', keytype = "PROBEID", multiVals = 'first')
rownames(cris_symmat) = tmpgenes
cris_symmat = cris_symmat[-(which(is.na(rownames(cris_symmat)))),]
cris_symmat = collapseRows(cris_symmat, unique(rownames(cris_symmat)), unique(rownames(cris_symmat)))
cris_symmat = cris_symmat$datETcollapsed
dim(cris_symmat) #21359 genes

###############################################################################
## ICGC data:

## (may have to use the server for extra RAM (ulimit -s 20480))

files = c('BLCA', 'CESC', 'COAD', 'KIRP', 'LIHC', 'READ', 'SKCM', 'UCEC')
cancerfiles = paste('./raw/raw-cancer-data/exp_seq.', files, sep='')
cancerfiles = paste(cancerfiles, '-US.tsv', sep='')
cancerbmi = paste('./clindata/TCGA/', tolower(files), sep='')
cancerbmi = paste(cancerbmi, '_clinical_patient.txt', sep='')
for (i in 1:length(files)) {
	## process sequence data:
	seq = read.table(cancerfiles[i], sep = '\t', header=T)
	seq = tbl_df(seq)
	dup = duplicated(paste(seq$submitted_sample_id,seq$gene_id,sep=''))
	seq = seq[!dup,c('submitted_sample_id','gene_id','raw_read_count')]
	genes = unique(seq$gene_id)
	seq  = spread(seq,submitted_sample_id,raw_read_count)
	rownames(seq) = genes
	seq  = data.matrix(seq)
	seq = seq[-1,-1]
	seq = icgc_to_tcga(t(seq))
    if (length(rownames(seq)) > length(unique(rownames(seq)))) {
        seq = collapseRows(seq, unique(rownames(seq)), unique(rownames(seq)))
        seq = seq$datETcollapsed
    }

	## process clinical data:
	clin = read.csv(cancerbmi[i], sep = '\t', skip=1, header=T)
	clin = clin[-1,]

	## process bmi data:
	bmi = clin[,c('height', 'weight')]
	rownames(bmi) = clin$bcr_patient_barcode
	if (length(which(bmi[,1] == "[Not Available]")) > 0){
		bmi = bmi[-which(bmi[,1] == "[Not Available]"),]
	}
	if (length(which(bmi[,2] == "[Not Available]")) > 0){
		bmi = bmi[-which(bmi[,2] == "[Not Available]"),]
	}

	seq = seq[which(rownames(seq) %in% rownames(bmi)),]
	bmi = bmi[which(rownames(bmi) %in% rownames(seq)),]

	## assign variable names:
	assign(files[i], seq)
	assign(paste(files[i], 'clin', sep=''), clin)
	assign(paste(files[i], 'bmi', sep=''), bmi)
}

###############################################################################
## Pathway data base:

#Import Human Gene Symbols
SYMBOL.list<-as.list(org.Hs.egSYMBOL)

KEGG.list<-as.list(org.Hs.egPATH) ##Import KEGG pathways:
names(KEGG.list)<-unlist(SYMBOL.list) #name KEGG pathway lists with corresponding gene symbols
keggpath<-as.list(KEGGPATHID2NAME) #mapping KEGG path IDs to human read pathway name

GO.list<-as.list(org.Hs.egGO) ##Import GO pathways:
tmp = GO.list #... and reformat so matching KEGG.list

# pull out the GOID
for (i in 1:length(tmp)) {
    ind = names(tmp)[i]
    if(class(tmp[[ind]]) == 'list') {
        tmp[[ind]] = names(tmp[[ind]])
    }
}
tmp = tmp[which(!is.na(tmp))]##filter out the NA values
GO.list = tmp
names(GO.list)<-SYMBOL.list[names(GO.list)] #name GO pathway lists with corresponding gene symbols
gopath<-as.list(GOTERM) #mapping GO IDs to human read pathway name
#... and reformat so matching KEGG.list
tmp = list()
for(i in 1:length(gopath)) tmp[[i]]<-gopath[[i]]@Term
names(tmp) = names(gopath)
gopath = tmp

#Import Human Reactome pathways
reactome.list <-as.list(reactomeEXTID2PATHID)
reactome.list = reactome.list[order(as.numeric(names(reactome.list)))] #Sort the names in the list:
reactome.list = reactome.list[which(names(reactome.list) %in% names(SYMBOL.list))] #get paths that have gene symbols:
names(reactome.list) = unlist(SYMBOL.list[names(reactome.list)]) #rename the entrez gene ID into gene symbol
reactomepath <- as.list(reactomePATHID2NAME) #mapping readtome path IDs to human read pathway name
reactomepath = reactomepath[grep('Homo sapiens',reactomepath)] #pull out all human-related pathways
reactomepath = lapply(reactomepath, function(x) gsub('Homo sapiens: ', '', x)) #cut out the 'Homo sapiens: ' bit so it's only the pathway names.


