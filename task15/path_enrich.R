library(DESeq)
library(affy)
library(limma)
library(GO.db)
library(KEGG.db)
library(org.Hs.eg.db)
library(reactome.db)
library(WGCNA)
library(dplyr)
library(hgu133a.db)
library(tidyr)

# KEGG matrix
genenames = colnames(BLCAraw) #Get all the gene names in the cancer data
pathnames = unique(unlist(keggpath)) #Get all the pathways in the keggpath list
tmp = matrix(nrow = length(genenames), ncol = length(pathnames)) #make an empty matrix to fill
colnames(tmp) = pathnames
rownames(tmp) = genenames

for(i in 1:length(rownames(tmp))) {
    v = colnames(tmp) %in% lookup(rownames(tmp)[i], "KEGG")
    if(length(v) > 0) {
        tmp[rownames(tmp)[i],] = v
    }
}

keggTFmat = ifelse(tmp == F, 0, 1)

# Reactome matrix
pathnames = unique(unlist(reactomepath))
tmp = matrix(nrow = length(genenames), ncol = length(pathnames))
colnames(tmp) = pathnames
rownames(tmp) = genenames

for(i in 1:length(rownames(tmp))) {
    v = colnames(tmp) %in% lookup(rownames(tmp)[i], "reactome")
    if(length(v) > 0) {
        tmp[rownames(tmp)[i],] = v
    }
}

reactomeTFmat = ifelse(tmp == F, 0, 1)

# GO matrix
# For GO, there are too many pathways in total, so will have to make a matrix
# using the pathways relevant to the genes in the cancer data.
goentrez = GO.list[genenames] ## pull out pathways that are relevant
goentrez = goentrez[which(!is.na(names(goentrez)))] #remove NAs
goentrez = unique(unlist(goentrez)) #get unique entrezID

pathnames = gopath[goentrez] ## pull out necessary paths using the entrezID
pathnames = unlist(pathnames)

tmp = matrix(nrow = length(genenames), ncol = length(pathnames))
rownames(tmp) = genenames
colnames(tmp) = pathnames

for(i in 1:length(rownames(tmp))) {
    v = colnames(tmp) %in% lookup(rownames(tmp)[i], "GO")
    if(length(v) > 0) {
        tmp[rownames(tmp)[i],] = v
    }
}

goTFmat = ifelse(tmp == F, 0, 1)

##remove pathways that have no genes involved (i.e. 0s)
keggTFmat = keggTFmat[,!apply(keggTFmat==0, 2, all)]
goTFmat = goTFmat[,!apply(goTFmat ==0, 2, all)]
reactomeTFmat = reactomeTFmat[,!apply(reactomeTFmat==0, 2, all)]

## save output
dput(keggTFmat, file = 'kegggenebypath.txt')
dput(reactomeTFmat, file = 'reactomegenebypath.txt')
dput(goTFmat, file = 'gogenebypath.txt')
