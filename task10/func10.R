# Functions from task 10:

#normVoom function from task09:
normVoom = function(x, design) {
    a = x[rowSums(cpm(x)) > 9,]
    b = calcNormFactors(a)
    a = voom(a, design, lib.size = colSums(a)*b)
    return(a)

}

#make_tt function from task09:
make_tt = function(x, model) {
    fit = lmFit(x, model)
    fit = eBayes(fit)
    tt = topTable(fit, coef = 2, adjust = 'BH', n = nrow(x))
    return(tt)

}

#pull_deg function from task09:
pull_deg = function(x, adj = F, y = 0.05) {
    if(adj) {
        x = x[x$adj.P.Val <= y,]

    } else {
        x = x[x$P.Value <= y,]

    }
    return(x)

}

#function to do DEG analysis n times, using randomly chosen samples
###################################################################
# files is a vector of the names of the variable
# samples is a matrix containing how many samples are in normal/overweight/obese group for each variable in files
# mat is an empty matrix of n by length(files)
# n is the total number of testing you want to carry out
# p is the p-value to be used for the pull_deg function
###################################################################
sample_test = function(files = files, samples = samples, mat, n = 100, p = 0.05) {
    for (j in 1:n) {
        l = as.list(gsub('mat','',files))
        names(l) = gsub('mat','',files)

        for (i in 1:length(files)) {
            group = sample(1:sum(samples[i,]), samples[i,3], replace=F) #pick samples randomly
            group = c(1:sum(samples[i,])) %in% group #make a vector out of the randomly chosen samples
            design = model.matrix(~group) #make model matrix using the randomly chosen samples.
            dat = get(files[i])
            dat = normVoom(dat, design) #normalise the data
            top = make_tt(dat, design) #make top table from the data
            l[[i]] = rownames(pull_deg(top, y = p))

        }

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

#attempt at optimising the above function
sample_test2 = function(files = files, samples = samples, mat, n = 100, p = 0.05) {
    originalList = as.list(gsub('mat','',files))
    names(originalList) = gsub('mat','',files)
    for (i in 1:length(files)) {
        originalList[[i]] = get(files[i])

    }

    for (j in 1:n) {
        l = originalList
        l = mclapply(seq_along(l), function(x) {
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

