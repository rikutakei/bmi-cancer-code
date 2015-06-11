######################################################################

                        ##R code for task2

######################################################################

## set working directory:
setwd('~/Documents/masters/data/task02')

## read the RNA seq data:
seq = read.table('exp_seq.UCEC-US.tsv', sep = '\t', header = T)

## the file format is useless, so use dplyr to manipulate it:
library(dplyr)

seq = tbl_df(seq) ## change data frame into tbl class

## extract all the raw count data and put it in a matrix:
initsamplenames = filter(seq, gene_id == 'PAN2')$submitted_sample_id
initgenenames = unique(seq$gene_id)
raw = matrix(nrow = 541, ncol = 20502)
rownames(raw) = initsamplenames
colnames(raw) = initgenenames
for (i in 1:nrow(raw)) {
    gene_name = raw[1,1]
    tmp = filter(seq, gene_id == gene_name) >%>
    select(raw)
}

