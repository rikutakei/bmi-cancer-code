#load bioconductor annotation packages
if(require("org.Hs.eg.db")){
  SYMBOL.list<-as.list(org.Hs.egSYMBOL)
} else{
  source("http://bioconductor.org/biocLite.R")
  biocLite("org.Hs.eg.db", suppressUpdates=T)
  library("org.Hs.eg.db")
  SYMBOL.list<-as.list(org.Hs.egSYMBOL)
}
if(!require("KEGG.db")){
  source("http://bioconductor.org/biocLite.R")
  biocLite("KEGG.db", suppressUpdates=T)
  library("KEGG.db")
}
if(!require("GO.db")){
  source("http://bioconductor.org/biocLite.R")
  biocLite("GO.db", suppressUpdates=T)
  library("GO.db")
}
if(!require("reactome.db")){
  source("http://bioconductor.org/biocLite.R")
  biocLite("reactome.db", suppressUpdates=T)
  library("reactome.db")
}
#Import Human KEGG pathways
KEGG.list<-as.list(org.Hs.egPATH)
#Import Human GO pathways
GO.list<-as.list(org.Hs.egGO)
#... and reformat so matching KEGG.list
GO.list2<-list()
for(i in 1:length(GO.list))  GO.list2[[i]]<-names(GO.list[[i]])
rm(GO.list)
#Import Human Gene Symbols
SYMBOL.list<-as.list(org.Hs.egSYMBOL)
#name KEGG pathway lists with corresponding gene symbols
names(KEGG.list)<-unlist(SYMBOL.list)
#name GO pathway lists with corresponding gene symbols
names(GO.list2)<-SYMBOL.list[names(GO.list2)]
#mapping KEGG path IDs to human read pathway name
path<-as.list(KEGGPATHID2NAME)
#mapping GO IDs to human read pathway name
path2<-as.list(GOTERM)
#... and reformat so matching KEGG.list
pathx<-list()
for(i in 1:length(path2)) pathx[[i]]<-path2[[i]]@Term
names(pathx)<-names(path2)
path2<-pathx
rm(pathx)

#download WikiPathways gpml format and extract Entrez IDS for each pathway (shell script)
#import WikiPathways into R, extract human read path names (path3), path IDs, and WikiPathways for each Entrez ID / Gene Symbol (R script)
system("bash download_WikiPathways.sh")
  #source("download_WikiPathways.R")
load("WikiPathways.db.RData")

#Import Human Reactome pathways
reactome.list<-as.list(reactomeEXTID2PATHID)
#name reactome pathway lists with corresponding gene symbols
reactome.list<-reactome.list[match(names(SYMBOL.list), names(reactome.list))]
names(reactome.list)<-unlist(SYMBOL.list)
#mapping readtome path IDs to human read pathway name
path4<-as.list(reactomePATHID2NAME)
#extract human pathways
path4<-lapply(path4[grep("Homo sapiens", path4)], function(x) strsplit(x, split=": ")[[1]][2:length(strsplit(x, split=": ")[[1]])])
path4<-lapply(path4, function(x) if(length(x)>1) paste0(x[1], ": ", x[2]) else x)
#save files
save.image("Gene_Set_Data.RData")
