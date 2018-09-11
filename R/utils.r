# filter pairwise order to help speed up clonal ordering
# 

#' get list of acceptable pairwise ordering (parent -> child)
#' using the vaf cutoff threshold
#' 
getOrderPassingThreshold <- function(v, vaf.col.names, vaf.diff=0.05, method='mean', cluster.col.name='cluster'){
    cat('Buidling parent child order using diff in VAF cutoff\n')
    clusters = unique(sort(as.integer(v$cluster[!is.na(v$cluster)])))
    n = length(clusters)
    co = c()
    for (i in 1:(n-1)){
        for (j in (i+1):n){
            co = rbind(co, c(clusters[i],clusters[j]))
            co = rbind(co, c(clusters[j],clusters[i]))
        }
    }
    cat(nrow(co), 'pairwise order(s) generated.\n')

    vafs = estimate.clone.vaf(v, cluster.col.name, vaf.col.names, method=method)
    rownames(vafs) = as.character(vafs$cluster)

    violated = apply(vafs[as.character(co[,2]),vaf.col.names] -
        vafs[as.character(co[,1]),vaf.col.names] > vaf.diff, 1, any, na.rm=T)
    cat(sum(violated), ' pairwise order(s) removed due to higher VAF in child.\n')
    cat(sum(!violated), ' pairwise order(s) retained.\n')
    table(violated)
    co = co[!violated,,drop=F]
    colnames(co) = c('parent', 'child')
    return(as.data.frame.matrix(co))
}

# print scores of models
printScores <- function(x){
    # print ppos for individual samples and the product
    cat('*** Product of ppos:\n')
    print(x$matched$probs)

    # min ppos of all clones per samples
    for (i in 1:x$num.matched.models){
        # pneg matrix
        p = x$matched$clone.ccf.pvalues[[i]][, vaf.col.names];
        # print ppos
        cat('*** Min ppos:\n')
        cat('model ', i, '\tmin ppos = ', min(1-p, na.rm=T), '\t' , min(p[p > 0], na.rm=T), '\n')
    }
}

