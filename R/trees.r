# Utilities for converting and plotting clonal evolution tree
# Motivated by examples from: http://stackoverflow.com/questions/28163979/how-can-i-plot-a-tree-and-squirrels-in-r

#' Create a tree from the merged.tree data.frame in clonevol output
#' @description
#' @param t: merged.tree data frame from clonevol
#' @param branch.lens: named vector of branch lengths, named according
#' to the cluster/clone ID/name. If not provided, all branches will receive
#' a same length
#' @output Returns the same data.frame with additional columns, including:
#' $branches --> symbols of the branches corresponding to the clone
#' $blengths --> branch lengths
#' 
convert.clone.to.branch <- function(t, branch.lens = NULL){
    # symbols used for each branch
    syms = c(seq(1,9), unlist(strsplit('abcdefghijklmnopqrstuvwxyz', '')))
    t = t[!is.na(t$parent) & !is.na(t$excluded) & !t$excluded,]
    t$branches = NA
    t$blengths = NA
    rownames(t) = t$lab
    grow.tree <- function(t, lab, parent.symbol=''){
        #print(paste0('---', lab))
        if (t[lab,'parent'] == '-1'){
            t[lab,'branches'] = 'Y'
        }
        children = !is.na(t$parent) & t$parent == lab
        if (any(children)){
            children.labs = t$lab[children]
            num.children = length(children.labs)
            if (num.children == 1 && parent.symbol != ''){
                children.syms = '0'
            }else{
                children.syms = syms[1:num.children]
            }
            children.syms = paste0(parent.symbol, children.syms)
            t$branches[children] = children.syms
            #print(children.labs)
            #print(children.syms)
            for (i in 1:length(children.labs)){
                t = grow.tree(t, children.labs[i], children.syms[i])
            }
        }
        return(t)
    }
    tg = grow.tree(t, t$lab[!is.na(t$parent) & t$parent == '-1'])
    #tg$samples.with.nonzero.cell.frac = gsub(',+$', '', gsub('\\s*:\\s*[^:]+(,|$)', ',', tg$sample.with.nonzero.cell.frac.ci))
    tg$samples.with.nonzero.cell.frac = gsub(',+$', '', gsub('°[^,]+(,|$)', '', gsub('\\s*:\\s*[^:]+(,|$)', ',', tg$sample.with.cell.frac.ci)))
    if (is.null(branch.lens)){
        tg$blengths = 5
    }else{
        tg$blengths = branch.lens[tg$lab]
    }
    return(tg)
    
}

#' Create trees for all merged.trees in clonevol output
#' @description
#' @param x: output of infer.clonal.models
#' @param cluster.col: cluster column name (used to count variants
#' in each cluster/clone, and scale tree branch length
#' @param branch.scale: values=c('none', 'sqrt', 'log2'), scale branch
#' length by corresponding transformations. Note, branch length will
#' be estimated by the number of variants in each cluster/clone
#
convert.merged.tree.clone.to.branch <- function(x, cluster.col='cluster', branch.scale='none'){
    num.trees = length(x$matched$merged.trees)
    if (num.trees > 0){
        res = list()
        blens = table(x$variants[[cluster.col]])
        if (branch.scale=='log2'){
            blens = log2(blens)
        }else if (branch.scale=='sqrt'){
            blens = sqrt(blens)
        }
        for (i in 1:num.trees){
            res[[i]] = convert.clone.to.branch(x$matched$merged.trees[[i]], blens)
        }        
        #x$matched$merged.trees.clone.as.branch = res
        x$matched$merged.trees = res
    }
    return(x)
}

#' Plot tree
#' 
plot.tree.clone.as.branch <- function(mt, angle=15, event.sep.char=','){
    mt$events = gsub(',', '\n', mt$events)
    g <- germinate(list(trunk.height=32,
                       branches=mt$branches,
                       lengths=mt$blengths,
                       branch.colors=mt$color,
                       node.colors=mt$color,
                       node.labels=mt$lab,
                       node.texts=mt$samples.with.nonzero.cell.frac,
                       branch.texts=mt$events),
                       angle=angle
    )
}




