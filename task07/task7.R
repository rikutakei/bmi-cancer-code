######################################################################

                        ##R code for task6

######################################################################

## Now looking at other cancer data.

######################################################################
## pre-analysis stuff
setwd('~/Documents/masters/data/task07/')

library(data.table)
library(tidyr)
######################################################################

#eq<-fread(paste(path,'exp_seq.SKCM-US.tsv',sep=''))

#dup<-duplicated(paste(expSeq$gene_id,expSeq$submitted_sample_id,sep=''))

#mik<-expSeq[!dup,c("submitted_sample_id","gene_id","raw_read_count"),with=FALSE]
#rm('expSeq')
#aa<-spread(mik,submitted_sample_id,raw_read_count)
#expSeqMat<-round(as.matrix(aa[,-1,with=FALSE]),4)
#rownames(expSeqMat)<-as.vector(data.frame(aa[,1,with=F])[,1])
#save(list="expSeqMat",file='expSeqMat-skcm.RData')

cancertypes = c('BLCA','CESC','COAD','KIRP','LIHC','READ','SKCM','UCEC')
files = readLines('files.txt')

path = getwd()
path = gsub('task07','raw/raw-cancer-data/',path)

for(i in 1:length(files)){
     files[i] = paste(path,files[i],sep='')
}

for(i in 1:length(files)) {
     assign(cancertypes[i],fread(files[i]))
}

#Convert the TCGA data into matrix form using tidyr:spread function

for(i in 1:length(cancertypes)) {
    t = get(cancertypes[i])
    t = t[,c('submitted_sample_id','gene_id','raw_read_count'),with=F] # pull out the required information
    dup = duplicated(paste(t$submitted_sample_id,t$gene_id,sep='')) #find duplicated genes
    t = t[!dup,] # remove duplicated genes
    t = spread(t,submitted_sample_id,raw_read_count) # make matrix
    rownames(t) = t[,gene_id] # change rownames into gene names
    t = data.matrix(t) #change it into matrix from data frame/table
    t = t[-1,-1] # remove the gene column and ? row
    txt = paste(cancertypes[i],'mat',sep='')
    assign(txt,t)
}

for(i in 1:length(cancertypes)) {
    txt = paste(cancertypes[i],'mat',sep='')
    t = get(txt)
    dput(t,file=paste(txt,'.txt',sep=''))
}
