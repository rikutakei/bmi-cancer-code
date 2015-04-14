######################################################################

                        ##R code for task1

######################################################################

## set working directory, etc:
setwd('~/Documents/masters/data/task01/patient_data/')
file = readLines('names.txt')

## read in the patient data from the downloaded files:
raw = list()
for (i in 1:length(file)){
    txt = gsub("*.txt","",file[i])
    assign(txt,read.csv(file[i],sep = '\t', header = T))
}

## make a list of 33 and put all the patient data in this list:
dat_list = vector("list",33) ## make an empty list with 33 items
txt = gsub("*.txt","",file)
for (i in 1:length(txt)) {
    dat_list[[i]] = get(txt[i])
}
names(dat_list) = txt ## rename the names of the list

## need to find which patient data has both weight and height data:
## First, need to see which cancer subtypes have weight and height columns
subt_ind = vector()
for (i in 1:length(dat_list)) {
    tmp_col = colnames(dat_list[[i]])
    c = c("weight_kg_at_diagnosis","height_cm_at_diagnosis")
    if (c[1] %in% tmp_col & c[2] %in% tmp_col) {
        subt_ind = append(subt_ind,i)
    }
}

## this gave me indices for 14 subtypes. when i vectorised all of the colnames
## from the list, there were no other weight-like variables so, i think i've 
## got all of the subtypes. (see the code below)
v = vector()
for (i in 1:33){
    v = append(v,colnames(dat_list[[i]]))
}
t = table(v)
t = as.matrix(t)
## rownames(t) ## This will give all of the names of the variables.

## Now I'll have to find out which patient has both the weight and height
## data, within these cancer subtypes.

## subt_ind contains all of the indices of the subtypes we need
bmi_pat_ind = vector("list", 14)
for (i in 1:length(subt_ind)) {
    tmp = as.vector(dat_list[[subt_ind[i]]]$height_cm_at_diagnosis)
    tmp = as.numeric(tmp)
    tmp2 = as.vector(dat_list[[subt_ind[i]]]$weight_kg_at_diagnosis)
    tmp2 = as.numeric(tmp2)
    
    bmi_pat_ind[[i]] = which((tmp > 0) & (tmp2 > 0))
}

## Make a table/list of cancer subtypes only with patients with height
## and weight data. 

bmi_pat_list = bmi_pat_ind ## make a copy (will be overwritten in the for loop)
for (i in 1:length(subt_ind)){
    ind = bmi_pat_ind[[i]]
    bmi_pat_list[[i]] = dat_list[[subt_ind[i]]][ind,] ## copy the data of patients with height and weight
}

## dput(x = bmi_pat_list, file = 'BMI_data.txt') ## to send data to Mik
