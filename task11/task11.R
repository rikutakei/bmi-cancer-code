######################################################################

                        ##R code for task 11

######################################################################

## Do similar stuff as task 10, but with pathway enrichment analysis.

######################################################################
## pre-analysis stuff
setwd('~/Documents/masters/data/task11/')

library(data.table)
library(gplots)

#libraries for data wrangling
library(WGCNA)
library(hgu133a.db)
library(tidyr)
library(dplyr)

#libraries for DEG analysis
library(limma)
library(edgeR)
library(DESeq)
library(affy)

#libraries for pathway analysis
library(org.Hs.eg.db)
library(KEGG.db)
library(GO.db)
library(reactome.db)

#pathway analysis packages from bioconductor
source('http://bioconductor.org/biocLite.R')
biocLite(c('topGO','GOstats','ReactomePA','KEGGprofile','PathNet','safe'))
library(topGO)
library(GOstats)
library(ReactomePA)
library(KEGGprofile)
library(PathNet)
library(safe)

######################################################################

## R code from Tom:

##Need to search for the pathways using the gene symbols
## obtained from DEG (or somewhere else).

## Search/match the gene symbols with the entrez ID by pulling out the
## ID using the ****.list variable.

# Once you know the entrez ID, you can search the pathways that match
# the entrez ID, using the path* variables.

######################################################################

#Import Human Gene Symbols
SYMBOL.list<-as.list(org.Hs.egSYMBOL)

##Import KEGG pathways:
KEGG.list<-as.list(org.Hs.egPATH)
#name KEGG pathway lists with corresponding gene symbols
names(KEGG.list)<-unlist(SYMBOL.list)
#mapping KEGG path IDs to human read pathway name
keggpath<-as.list(KEGGPATHID2NAME)

##Import GO pathways:
GO.list<-as.list(org.Hs.egGO)
#... and reformat so matching KEGG.list
tmp = GO.list

## the for-loop below gives you duplicate values -> DON'T USE IT (for now at least)
#for(i in 1:length(GO.list))  tmp[[i]]<-names(GO.list[[i]] ) ## pull out GOID from GO.list

#Use this to pull out the GOID
for (i in 1:length(tmp)) {
    ind = names(tmp)[i]
    if(class(tmp[[ind]]) == 'list') {
        tmp[[ind]] = names(tmp[[ind]])
    }
}
tmp = tmp[which(!is.na(tmp))]##filter out the NA values
GO.list = tmp
#name GO pathway lists with corresponding gene symbols
names(GO.list)<-SYMBOL.list[names(GO.list)]
#mapping GO IDs to human read pathway name
gopath<-as.list(GOTERM)
#... and reformat so matching KEGG.list
tmp = list()
for(i in 1:length(gopath)) tmp[[i]]<-gopath[[i]]@Term
names(tmp) = names(gopath)
gopath = tmp

#Import Human Reactome pathways
reactome.list <-as.list(reactomeEXTID2PATHID)
#Sort the names in the list:
reactome.list = reactome.list[order(as.numeric(names(reactome.list)))]
#get paths that have gene symbols:
reactome.list = reactome.list[which(names(reactome.list) %in% names(SYMBOL.list))]
#rename the entrez gene ID into gene symbol
names(reactome.list) = unlist(SYMBOL.list[names(reactome.list)])
#mapping readtome path IDs to human read pathway name
reactomepath <- as.list(reactomePATHID2NAME)
#pull out all human-related pathways
reactomepath = reactomepath[grep('Homo sapiens',reactomepath)]
#cut out the 'Homo sapiens: ' bit so it's only the pathway names.
reactomepath = lapply(reactomepath, function(x) gsub('Homo sapiens: ', '', x))

######################################################################

#some codes from task 10

#load all the data:
files = readLines('../task09/files.txt')
files = paste('../task09/',files,sep='')
bmifiles = gsub('mat','bmi',files)

#create variables for cancer data matrix
for (i in 1:length(files)) {
    txt = gsub('t9.txt','',files[i])
    txt = gsub('../task09/','',txt)
    assign(txt, t(dget(files[i]))) #transpose the matrix so the columns are samples
    files[i] = txt #rename files to its variable names

    txt = gsub('t9.txt','',bmifiles[i])
    txt = gsub('../task09/','',txt)
    assign(txt, dget(bmifiles[i]))
    bmifiles[i] = txt #rename files to its variable names
}

#make a matrix of the sample sizes of each cancer type
samples = matrix(nrow = 8, ncol = 3)
colnames(samples) = unique(BLCAbmi[,4])[c(1,3,2)]
rownames(samples) = gsub('mat','',files)
for (i in 1: length(files)) {
    samples[i,] = as.vector(table(get(bmifiles[i])[,4]))[c(1,3,2)]
}


#source the file to get all the functions from task10.R
source('~/Documents/codes/bmi-cancer-code/task10/func10.R')

set.seed(1) ##set the seed for reproducibility

#make a gene-by-pathway matrix to use it in the pathway enrichment analysis

#first, make a function to pull out pathways from the database, using the genes in the cancer data
#genename is the name of the gene (character), databasetype is one of "KEGG", "GO" or "reactome" (character)
lookup = function(genename, databasetype) {
    ## check which database you're looking up in
    if (databasetype == "KEGG") {
        genetoentrez = KEGG.list
        entreztopath = keggpath
    } else if (databasetype == "GO") {
        genetoentrez = GO.list
        entreztopath = gopath
    } else {
        genetoentrez = reactome.list
        entreztopath = reactomepath
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

#make a function to make a conditional vector describing which pathway(s) is present for a given gene
#makeTFvector = function(path, allpath) {
    #v = allpath %in% path
    #return(v)
#}


#KEGG matrix
pathnames = unique(unlist(keggpath)) #Get all the pathways in the keggpath list
genenames = rownames(BLCAmat) #Get all the gene names in the cancer data
tmp = matrix(nrow = length(genenames), ncol = length(pathnames)) #make an empty matrix to fill
colnames(tmp) = pathnames
rownames(tmp) = genenames

for(i in 1:length(rownames(tmp))) {
    v = colnames(tmp) %in% lookup(rownames(tmp)[i], "KEGG")
    if(length(v) > 0) {
        tmp[rownames(tmp)[i],] = v
    }
}

#GO matrix


#reactome matrix





#Pathway enrichment analysis function
#It should take in DEGs and pull out relevant data from a pre-made gene-by-pathway matrix.















































sample_test2 = function(files = files, samples = samples, mat, n = 100, p = 0.05) {
    originalList = as.list(gsub('mat','',files))
    names(originalList) = gsub('mat','',files)
    for (i in 1:length(files)) {
        originalList[[i]] = get(files[i])
    }

    for (j in 1:n) {
        l = originalList
        l = mclapply(seq_along(l), function(x) { ##seq_along gives you the length of the list (i.e. the indices of the cancer subtypes in the list)
            group = sample(1:sum(samples[x,]), samples[x,3], replace=F) #pick samples randomly
            group = c(1:sum(samples[x,])) %in% group #make a vector out of the randomly chosen samples
            design = model.matrix(~group) #make model matrix using the randomly chosen samples.
            dat = normVoom(l[[x]], design) #normalise the data
            top = make_tt(dat, design) #make top table from the data
            l[[x]] = rownames(pull_deg(top, y = p))
        })

        t = unlist(l)
        t = table(table(t))
        v = c(rep(0,8))

        for (i in 1:length(t)) {
            v[i] = t[i]
        }

        mat[j,] = v
    }
    return(mat)
}





































