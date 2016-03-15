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

keggTFmat = ifelse(tmp == F, 0, 1)

#reactome matrix

pathnames = unique(unlist(reactomepath)) #Get all the pathways in the reactomepath list
tmp = matrix(nrow = length(genenames), ncol = length(pathnames)) #make an empty matrix to fill
colnames(tmp) = pathnames
rownames(tmp) = genenames

for(i in 1:length(rownames(tmp))) {
    v = colnames(tmp) %in% lookup(rownames(tmp)[i], "reactome")
    if(length(v) > 0) {
        tmp[rownames(tmp)[i],] = v
    }
}

reactomeTFmat = ifelse(tmp == F, 0, 1)

#GO matrix
#for GO, there are too many pathways in total, so will have to make a matrix
# using the pathways relevant to the genes in the cancer data.
goentrez = GO.list[genenames] ## pull out pathways that are relevant
goentrez = goentrez[which(!is.na(names(goentrez)))] #remove NAs
genenames = names(goentrez) #update genenames
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

##save the files as txt
dput(keggTFmat, file = 'kegggenebypath.txt')
dput(reactomeTFmat, file = 'reactomegenebypath.txt')
dput(goTFmat, file = 'gogenebypath.txt')

#Pathway enrichment analysis function
#It should take in DEGs and pull out relevant data from a pre-made gene-by-pathway matrix.
# deg should be in the top table format (containing only significant genes)
# db is the type of database you want to use in string format (i.e. 'KEGG', 'reactome', or 'GO')

pathenrich = function(deg, db, gn = genenames) {
    #check which database you're using
    if (db == "KEGG") {
        db = keggTFmat
    } else if (db == "GO") {
        db = goTFmat
    } else {
        db = reactomeTFmat
    }

    ##pull out relevant genes from the TF mat of the database
    #mat = db[which(rownames(db) %in% rownames(deg)),]
    ##remove columns with 0s (i.e. pathways with no DEGs)
    #mat = mat[,!apply(mat==0, 2, all)]

    # matrix to dump the p-values of the pathway enrichment
    v = matrix(nrow = ncol(db), ncol = 1, dimnames = list(colnames(db),'p.value'))
    x = gn %in% rownames(deg)
    for (i in 1:ncol(db)){
        y = gn %in% rownames(db)[db[,i]==1]
        if(length(table(y)) > 1 ) { ##do Fisher's test only if there is a '1' inthere
            f = fisher.test(table(x,y))
            v[i,1] = f$p.value
        } else {
            v[i,1] = NA
        }
    }

    return(v)
}

#set genenames to all 20501 genes in the  cancer data
genenames = rownames(BLCAmat)

##This function does path enrichment analysis n times on a single cancer type
# dat is the cancer data, sample_dat is the sample distribution in original data (normal/overweight/obese), gn is a vector containing all the genes in the cancer data (i.e. rownames(BLCAmat)), n is the number of analyses to be carried out.
npathenrich = function(dat, sample_dat, gn = genenames, n = 100, db = "KEGG"){
    for(i in 1:n) {
        #pick samples randomly
        group = sample(seq(sum(sample_dat)),sample_dat[3], replace = F)
        #make a vector out of the randomly chosen samples
        group = seq(sum(sample_dat)) %in% group
        #make model matrix using the randomly chosen samples.
        design = model.matrix(~group)
        #normalise the data
        x = normVoom(dat, design)
        #make top table from the data
        tt = make_tt(x, design)
        # pull out DEGs
        deg = pull_deg(tt)

        tmppath = pathenrich(deg, db, gn = genenames)

        if(i == 1) {
            mat = tmppath
        } else if (i > 1){
            mat = cbind(mat, tmppath)
        }
    }
    return(mat)
}


#############################################################################
# start pathway enrichment analysis:
#############################################################################

#begin by doing the pathway enrichment analysis on the ob vs lean/ov on each
# cancer type

files = gsub('mat','bmi',files)

originalpathlist = originalList

## do pathway enrichment analysis on each cancer type
originalpathlist = lapply(seq_along(originalpathlist), function(x) {
                       group = get(files[x])[,4]
                       group = ifelse(group == 'obese', 'obese', 'normal/overweight')
                       model = model.matrix(~group)
                       dat = originalpathlist[[x]]
                       dat = normVoom(dat, model)
                       tt = make_tt(dat, model)
                       deg = pull_deg(tt)
                       path = pathenrich(deg, db = "GO", gn = rownames(dat))
})

## adjust the p-values using FDR:
originalpathlist = lapply(originalpathlist, function(x) {
                        p.adj = p.adjust(x, method = 'BH', n = nrow(x))
                        x = cbind(x, p.adj)
})

## pull out the names of enriched paths:
enrpaths = lapply(originalpathlist, function(x) {
                  ind = which(x[,2] <= 0.05)
                  rownames(x)[ind]
})

names(enrpaths) = gsub('bmi', '', files)

## ready to do the pathway enrichment analysis:

set.seed(1) ##set the seed for reproducibility

#make copy just in case
ntestrun = originalList

# repeat the test run, but with npathenrich function instead of pathenrich:
ntestrun = lapply(seq_along(ntestrun ), function(x) {
                          dat = ntestrun [[x]]
                          path = npathenrich( dat, sample_dat = samples[x,] ,n = 10, db = "GO", gn = rownames(dat))
})

ntestrun2 = ntestrun

ntestrun2 = lapply(ntestrun2 , function(x) {
                   tmp = apply(x, 2, function(y) p.adjust(y, method = 'BH', n = length(y)))
                   colnames(tmp) = paste('adj.', colnames(tmp), sep = '')
                   return(tmp)
})


























