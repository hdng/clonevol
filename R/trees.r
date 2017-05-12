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
convert.clone.to.branch <- function(t, branch.lens = NULL,
    merged.tree.node.annotation='sample.with.nonzero.cell.frac.ci'){
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
            if (num.children == 1){# && parent.symbol != ''){
                # one child, grow straight branch (no left right)
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
    if (merged.tree.node.annotation=='sample.with.nonzero.cell.frac.ci'){
        #tg$samples.with.nonzero.cell.frac = gsub(',+$', '',
        #    gsub('\\s*:\\s*[^:]+(,|$)', ',', tg$sample.with.nonzero.cell.frac.ci))
        tg$samples.with.nonzero.cell.frac = gsub(',+$', '', gsub('°[^,]+(,|$)', '',
             gsub('\\s*:\\s*[^:]+(,|$)', ',', tg$sample.with.cell.frac.ci)))
    }else{
        cat(paste0('WARN: merged.tree.node.annotation = ',
            merged.tree.node.annotation, ' not supported! No node annotation made.\n'))
    }
    if (is.null(branch.lens)){
        tg$blengths = 5
    }else{
        tg$blengths = branch.lens[tg$lab]
    }
    # color founding clone of met with diff. border
    #tg$node.border.color = ifelse(
    #    grepl('*', gsub('*P', '', tg[[merged.tree.node.annotation]], fixed=T), fixed=T),
    #    'red', 'black')
    tg$node.border.color = 'black'
    tg$node.border.width = 1
    tg$branch.border.color = 'white'
    tg$branch.border.linetype = 'solid'
    tg$branch.border.width = 0.5

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
convert.merged.tree.clone.to.branch <- function(x, cluster.col='cluster',
                                                branch.scale='none'){
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

        res = list()
        for (i in 1:length(x$matched$trimmed.merged.trees)){
            res[[i]] = convert.clone.to.branch(x$matched$trimmed.merged.trees[[i]], blens)
        }        
        #x$matched$merged.trees.clone.as.branch = res
        x$matched$trimmed.merged.trees = res

    }
    return(x)
}

#' Plot tree
#' 
plot.tree.clone.as.branch <- function(mt, angle=15, branch.width=1, branch.text.size=0.3,
    node.size=3, node.label.size=0.75, node.text.size=0.5, event.sep.char=',', show.event=TRUE,
    tree.rotation=0, text.angle=NULL,
    tree.label=NULL, branch.border.width=NULL,...){
    if ('events' %in% colnames(mt)){
        mt$events = gsub(event.sep.char, '\n', mt$events)
    }else{
        mt$events = ''
    }

    if (!is.null(branch.border.width)){
        mt$branch.border.width = branch.border.width
    }
    if (!show.event){mt$events=''}
    g <- germinate(list(trunk.height=32,#not used
                       branches=mt$branches,
                       lengths=mt$blengths,
                       branch.colors=mt$color,
                       branch.border.colors=mt$branch.border.color,
                       branch.border.linetypes=mt$branch.border.linetype,
                       branch.border.widths=mt$branch.border.width,
                       node.colors=mt$color,
                       node.border.colors=mt$node.border.color,
                       node.border.widths=mt$node.border.width,
                       node.labels=mt$lab,
                       node.texts=mt$samples.with.nonzero.cell.frac,
                       branch.texts=mt$events),
                       angle=angle
    )
    plot(g, branch.width=branch.width, branch.text.size=branch.text.size,
        node.size=node.size, node.label.size=node.label.size,
        node.text.size=node.text.size, tree.rotation=tree.rotation, text.angle=text.angle, tree.label=tree.label)
}




