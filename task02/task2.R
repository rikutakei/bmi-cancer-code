######################################################################

                        ##R code for task2

######################################################################

## set working directory:
setwd('~/Documents/masters/data/task02')

## read the RNA seq data:
seq = read.table('exp_seq.UCEC-US.tsv', sep = '\t', header = T)

## the file format is useless, so use dplyr to manipulate it:
library(dplyr)
library(gplots)

seq = tbl_df(seq) ## change data frame into tbl class

## extract all the raw count data and put it in a matrix:
initsamplenames = filter(seq, gene_id == 'PAN2')$submitted_sample_id
initgenenames = unique(seq$gene_id)
raw = matrix(nrow = 541, ncol = 20502)
rownames(raw) = initsamplenames
colnames(raw) = initgenenames

for (i in 1:ncol(raw)) {
    gene = as.vector(initgenenames[i])
    tmpreadcount = filter(seq, gene_id == gene) %>%
    select(submitted_sample_id, raw_read_count)

    for (j in 1:nrow(tmpreadcount)) {
        sample = as.matrix(tmpreadcount[j,1])
        raw[sample, gene] = as.numeric(tmpreadcount[j, 2])
    }
    print(paste("Finished processing ", gene))
}

## save the matrix:
## dput(raw, file = 'raw_read_count.txt')

## need to SVD the raw read count and validate the Creighton et al results

## load obesity gene names and series matrix data
library(hgu133a.db) ## load library for microarray annotation

## read series matrix (microarray) data
## gse24185exp = read.table('./GSE24185_series_matrix.txt', header=T,
##                          sep='',comment.char='!', row.names=1)

## read in the obesity gene file:
obsgenenames = readLines('obsGenes.txt')

## check the valid columns to get from the hgu133a.db
columns(hgu133a.db)

## select the gene symbols that are assigned to the probes:
genesym = select(hgu133a.db, keys = obsgenenames, columns = 'SYMBOL')

## OR map the probe IDs to the symbols:
genesym2 = mapIds(hgu133a.db, keys = obsgenenames, column = 'SYMBOL', 
                  keytype = 'PROBEID', multiVals = 'first')

## NOTE: this only takes the first gene symbol that the probe ID is
## matched to, if the ID has multiple gene values.

## going to use the genesym2 data from here:

## get unique gene symbols and remove NA:
genesym2 = unique(genesym2)[2:length(genesym2)]
genesym2 = sort(genesym2) ## sort the symbols

## find which columns of the raw data I need:
colind = which(colnames(raw) %in% genesym2)

## extract the columns:
obsgeneraw = raw[,colind]

## there is a row of NAs in the data, so delete that:
which(is.na(obsgeneraw[,1] == T)) ## to find which row has the NAs
obsgeneraw = obsgeneraw[-221,] ## delete the row of NAs

## need to find which samples are actually from the primary tumour site
## load specimen data:
specdata = read.table(file = 'specimen.UCEC-US.tsv', sep = '\t',
                      header = T)
primtum = specdata[specdata$specimen_type == 'Primary tumour - solid tissue',]
primtumid = as.vector(primtum$submitted_specimen_id)

## scale the data so it has mean of 0, standard deviation of 1:
scaledobsgene = scale(obsgeneraw)
scaledobsgene = scaledobsgene[order(rownames(scaledobsgene)),] ## order by rownames

## make the max and min values to be 3 or -3:
scaledobsgene[scaledobsgene >= 3] = 3
scaledobsgene[scaledobsgene <= -3] = -3

## change the rownames so that it matches the TCGA sample names
specimenid = gsub(pattern = '-[[:alnum:]]{3}-[[:alnum:]]{4}-[[:alnum:]]{2}',
                               replacement = '', rownames(scaledobsgene))

## find which of these samples are from the primary tumour site:
ind = which(!(primtumid %in% specimenid))
scaledobsgene = scaledobsgene[-ind,]

## Import the clinical data to match the scaledobsgene.
clin = read.table('ucec_clinical_patient.txt', sep = '\t', skip = 1,
                  header = T)
clin = clin[-1,] ## remove random header

## get height and weight data:
hwdata = clin[,c('weight','height')]
rownames(hwdata) = clin$bcr_patient_barcode ## get rownames

## remove any samples with no height or weight data:
hwdata = hwdata[-(which(hwdata == '[Not Available]')),]
hwdata = hwdata[-(which(hwdata[,2] == '[Not Available]')),]

## convert to numerical matrix
rows = rownames(hwdata)
hwdata = as.matrix((hwdata))
hwdata = apply(hwdata, 2, as.numeric)
rownames(hwdata) = rows

## found a sample with a height of 66 cm (probs should be 166cm)
hwdata[hwdata[,2] < 100, 2] = 100 + hwdata[hwdata[,2] < 100, 2]

## calculate BMI:
bmi = hwdata[,1]/((hwdata[,2] / 100) ^ 2) ## need to convert height into m
bmi = round(bmi, 2)
hwdata = cbind(hwdata, bmi)

tcgaid = gsub(pattern = '-[[:alnum:]]{3}-[[:alnum:]]{3}-[[:alnum:]]{4}-[[:alnum:]]{2}',
                               replacement = '', rownames(scaledobsgene))

 #check which samples have BMI data:
#bmidataind = which(!(rownames(scaledobsgene) %in% rownames(hwdata)))
#scaledobsgene = scaledobsgene[-bmidataind,]

 #select the unique BMI samples:
#hwdata = hwdata[-(which(!(rownames(hwdata) %in% rownames(scaledobsgene)))),]


## need to identify which sample names are repeated:
tmprow = rownames(scaledobsgene) ## get rownames
repind = table(tmprow) > 1 ## find any names that are repeated
repnames = unique(tmprow)[repind] ## get the sample names that are repeated

## for loop to graph the replicates to decide which one should be
## discarded:
for (i in 1:length(repnames)) {
    ind = which(rownames(scaledobsgene) == repnames[i])
    heatmap.2(scaledobsgene[ind,], trace = 'none', col = 'bluered',
              scale = 'none', labRow = ind)
}


