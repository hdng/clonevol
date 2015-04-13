# Clonevol version 0.1.1
# Created by: Ha Dang <hdang@genome.wustl.edu>
# Date: Dec. 25, 2014
# Last modified:
#   Dec. 27, 2014  -- more than 2 samples model inference and plotting
#   Feb. 02, 2015  -- polyclonal model supported, bugs fixed.
#   Mar. 07, 2015  -- subclonal bootstrap test implementation, bugs fixed.
#
# Dependencies: igraph, ggplot2, grid, reshape2
#
# Purposes:
#  Infer and visualize clonal evolution in multi cancer samples
#  using somatic mutation clusters and their variant allele frequencies
#
# How-to-run example (see more examples at the end of this file):
#   c = read.table('clusters.tsv', header=T)
#   x = infer.clonal.models(c)
#   plot.clonal.models(x$models, out.dir='out', matched=x$matched,
#                     out.format='png')
#   plot.clonal.models(x$models, out.dir='out', matched=x$matched,
#                     out.format='pdf', overwrite.output=TRUE)
#
#

#' Create a data frame to hold clonal structure of a single sample
#'
#' @description Create a data frame ready for clonal structure enumeration. This
#' data frame will hold VAFs of variant clusters, descending order, together
#' with other additional columns convenient for clonal structure enumeration.
#'
#' @usage clones.df = make.clonal.data.frame(vafs, labels, add.normal=FALSE)
#'
#' @param vafs: VAFs of the cluster (values from 0 to 0.5)
#' @param labels: labels of the cluster (ie. cluster numbers)
#' @param add.normal: if TRUE, normal cell clone will be added as the superclone
#' (ie. the founding clones will be originated from the normal clone). This is
#' used only in polyclonal models where all clones can be separate founding
#' clones as long as their total VAFs <= 0.5. Default = FALSE
#' @param founding.label: label of the founding cluster/clone
#'
#' @details Output will be a data frame consisting of the following columns
#' vaf: orginial VAF
#' lab: labels of clusters
#' occupied: how much VAF already be occupied by the child clones (all zeros)
#' free: how much VAF is free for other clones to fill in (all equal original
#' VAF)
#' color: colors to plot (string)
#' parent: label of the parent clone (all NA, will be determined in
#' enumerate.clones)
#'
#' @examples
#' clones <- data.frame(cluster=c(1,2,3), sample.vaf=c(0.5, 0.3, 0.1))
#' clones.df <- make.clonal.data.frame(clones$sample.vaf, clones$cluster)
#'
#'
# TODO: Historically, the evolution tree is stored in data.frame for
# convenience view/debug/in/out, etc. This can be improved by using some
# tree/graph data structure
make.clonal.data.frame <- function (vafs, labels, add.normal=FALSE,
                                    founding.label=NULL, colors=NULL){
    v = data.frame(lab=as.character(labels), vaf=vafs, stringsAsFactors=F)
    if (is.null(colors)){
        #colors = c('#a6cee3', '#b2df8a', '#cab2d6', '#fdbf6f', '#fb9a99',
        #           '#d9d9d9','#999999', '#33a02c', '#ff7f00', '#1f78b4',
        #           '#fca27e', '#ffffb3', '#fccde5', '#fb8072', '#b3de69',
        #           'f0ecd7', rep('#e5f5f9',1))
        colors=get.clonevol.colors(nrow(v))
    }
    clone.colors = colors[seq(1,nrow(v))]
    v$color = clone.colors
    v = v[order(v$vaf, decreasing=T),]
    if (!is.null(founding.label)){
        #make founding clone first in data frame
        v1 = v[v$lab == founding.label,]
        v2 = v[v$lab != founding.label,]
        v = rbind(v1, v2)
    }
    # add dummy normal cluster to cover
    if (add.normal){
        v = rbind(data.frame(vaf=0.5, lab='0', color=colors[length(colors)],
                             stringsAsFactors=F), v)
    }

    v$parent = NA
    v$ancestors = '-'
    v$occupied = 0
    v$free = v$vaf
    v$free.mean = NA
    v$free.lower = NA
    v$free.upper = NA
    v$free.confident.level = NA
    v$free.confident.level.non.negative = NA
    v$p.value = NA
    v$num.subclones = 0
    v$excluded = NA
    #rownames(v) = seq(1,nrow(v))
    rownames(v) = v$lab
    #print(str(v))
    #print(v)
    return(v)
}

#' Check if clone a is ancestor of clone b in the clonal evolution tree
#' @usage
#' x = is.ancestor(v, a, b)
#' @param v: clonal structure data frame
#' @param a: label of the cluster/clone a
#' @param b: label of the cluster/clone a
#'
is.ancestor <- function(v, a, b){
    #cat('Checking if', a, '---ancestor-->', b, '\n')
    if (is.na(b) || b == '-1'){
        return(FALSE)
    }else{
        par = v$parent[v$lab == b]
        if(is.na(par)){
            return(FALSE)
        }else if (par == a){
            return(TRUE)
        }else{
            return(is.ancestor(v, a, par))
        }
    }
}


#' Enumerate all possible clonal structures for a single sample, employing the
#' subclonal test
#'
#' @description Enumerate all possible clonal structures for a single sample
#' using monoclonal (ie. the primary tumor is originated from a single
#' cancer cell) or polyclonal model (ie. the primary tumor can originate from
#' multi cancer cells)
#'
#' @usage
#'
#' enumerate.clones (v, sample, variants=NULL, subclonal.test.method='bootstrap'
#' , boot=NULL, p.value.cutoff=0.1)
#'
#' @param v: a data frame output of make.clonal.data.frame function
#'
#' @details This function return a list of data frames. Each data frame
#' represents a clonal structure, similar to the output format of
#' make.clonal.data.frame output but now have 'parent' column identified
#' indicating the parent clone, and other columns (free, occupied) which
#' can be used to calculate cellular fraction and plotting. The root clone
#' has parent = -1, clones that have VAF=0 will have parent = NA
#'
#' @examples --
#'
#'
enumerate.clones <- function(v, sample=NULL, variants=NULL,
                             founding.cluster = NULL,
                             ignore.clusters=NULL,
                             subclonal.test.method='bootstrap',
                             boot=NULL,
                             p.value.cutoff=0.05,
                             alpha=0.05,
                             min.cluster.vaf=0){
    vv = list() # to hold list of output clonal models
    findParent <- function(v, i){
        #print(i)
        if (i > nrow(v)){
            #debug
            #print(v)
            vv <<- c(vv, list(v))
        }else{
            #print(head(v))
            vaf = v[i,]$vaf
            if (!is.na(v[i,]$parent) && v[i,]$parent == '-1'){# root
                vx = v
                findParent(vx, i+1)
            }else if (v[i,]$excluded){
                vx = v
                vx$parent[i] = NA
                findParent(vx, i+1)
            }else{
                #for (j in 1:(i-1)){
                for (j in 1:nrow(v)){
                    parent.cluster = as.character(v[j,]$lab)
                    current.cluster = as.character(v[i,]$lab)
                    is.ancestor = is.ancestor(v, current.cluster,
                                              parent.cluster)
                    #print(v)
                    #cat('i=', i, 'j=', j, '\n')
                    if (i != j && !v$excluded[j] && !is.ancestor){
                        # assign cluster in row j as parent of cluster in row i
                        sub.clusters = as.character(c(v$lab[!is.na(v$parent) &
                                                    v$parent == parent.cluster],
                                                    current.cluster))
                        #print(str(v))
                        # debug
                        #cat('Testing...', sample, '-', j, parent.cluster,
                        #  'sub clusters:', sub.clusters, '\n')
                        t = subclonal.test(sample, parent.cluster, sub.clusters,
                                           boot=boot,
                                           cdf=v,
                                           min.cluster.vaf=min.cluster.vaf,
                                           alpha=alpha)

                        if(t$p.value >= p.value.cutoff){
                            vx = v
                            # debug
                            #cat(i, '<-', j, 'vaf=', vaf, '\n')
                            #print(head(v))
                            #print(head(vx))
                            vx$p.value[j] = t$p.value
                            vx$free.mean[j] = t$free.vaf.mean
                            vx$free.lower[j] = t$free.vaf.lower
                            vx$free.upper[j] = t$free.vaf.upper
                            vx$free.confident.level[j] =
                                t$free.vaf.confident.level
                            vx$free.confident.level.non.negative[j] =
                                t$free.vaf.confident.level.non.negative
                            vx$free[j] = vx$free[j] - vaf
                            vx$occupied[j] = vx$occupied[j] + vaf
                            vx$num.subclones[j] = length(sub.clusters)
                            #vx$parent[i] = vx[j,]$lab
                            vx$parent[i] = parent.cluster
                            vx$ancestors[i] = paste0(vx$ancestors[j],
                                paste0('#',parent.cluster,'#'))

                            # calculate confidence interval for vaf estimate of
                            # the subclone if it does not contain other
                            # subclones (will be overwrite later
                            # if subclones are added to this subclone)
                            #if (is.na(vx$free.lower[i])){
                            if (vx$num.subclones[i] == 0){
                                t = subclonal.test(sample,
                                       as.character(vx[i,]$lab),
                                       sub.clusters=NULL, boot=boot,
                                       cdf=vx,
                                       min.cluster.vaf=min.cluster.vaf,
                                       alpha=alpha)
                                vx$free.mean[i] = t$free.vaf.mean
                                vx$free.lower[i] = t$free.vaf.lower
                                vx$free.upper[i] = t$free.vaf.upper
                                vx$p.value[i] = t$p.value
                                vx$free.confident.level[i] =
                                    t$free.vaf.confident.level
                                vx$free.confident.level.non.negative[i] =
                                    t$free.vaf.confident.level.non.negative
                            }
                            findParent(vx, i+1)
                        }
                    }
                }
            }
        }
    }

    # exclude some cluster with VAF not significantly diff. from zero
    # print(v)
    for (i in 1:nrow(v)){
        cl = as.character(v[i,]$lab)

        # TODO: Find a better way
        # Currently this test won't work well to determine if a VAF is zero
        #t = subclonal.test(sample, parent.cluster=cl, sub.clusters=NULL,
        #                   boot=boot, min.cluster.vaf=min.cluster.vaf,
        #                   alpha=alpha)
        #v[i,]$excluded = ifelse(t$p.value < p.value.cutoff, TRUE, FALSE)

        # if the median/mean (estimated earlier) VAF < e, do
        # not consider this cluster in this sample
        v[i,]$excluded = ifelse(v[i,]$vaf <= min.cluster.vaf, TRUE, FALSE)
        #cat(sample, 'cluster ', cl, 'exclude p = ', t$p.value,
        #     'excluded=', v[i,]$excluded, '\n')
    }

    # also exlude clusters in the ignore.clusters list
    if (!is.null(ignore.clusters)){
        ignore.idx = v$lab %in% as.character(ignore.clusters)
        v$excluded[ignore.idx] = TRUE
        message('WARN: The following clusters are ignored:',
            paste(v$lab[ignore.idx], collapse=','), '\n')
    }
    #print(v)

    # if normal sample (0) is included, the normal sample
    # will be root (polyclonal model), otherwise find the
    # founding clone and place it first
    if (v[1,]$lab == 0){
        v[1,]$parent = -1
        findParent(v, 2)
    }else{
        #print(founding.cluster)
        if (is.null(founding.cluster)){
            max.vaf = max(v$vaf)
            roots = rownames(v)[v$vaf == max.vaf]
        }else{
            roots = rownames(v)[v$lab == founding.cluster]
        }
        # debug
        #cat('roots:', paste(roots, collapse=','), '\n')
        for (r in roots){
            #print(roots)
            vr = v
            vr[r,]$parent = -1
            #print(vr)
            findParent(vr,1)
        }
    }
    return(vv)
}



#' Check if two clonal structures are compatible (one evolve to the other)
#'
#' @description Check if two clonal structures are compatible (one evolve to
#' the other); ie. if structure v1 evolves to v2, all nodes in v2 must have
#' the same parents as in v1. This function returns TRUE if the two clonal
#' structures are compatible, otherwise, return FALSE
#'
#' @param v1: first clonal structure data frame
#' @param v2: first clonal structure data frame
#'
#' @details --
#' @examples --
#'
match.sample.clones <- function(v1, v2){
    compatible = TRUE
    for (i in 1:nrow(v2)){
        vi = v2[i,]
        parent2 = vi$parent
        if (is.na(parent2)){next}
        parent1 = v1[v1$lab == vi$lab,]$parent
        #debug
        #cat(vi$lab, ' par1: ', parent1, 'par2: ', parent2, '\n')
        if (!is.na(parent1) && parent1 != parent2){
            compatible = FALSE
            break
        }
    }
    return(compatible)
}

#' Draw a polygon representing a clone evolution, annotated with cluster label
#' and cellular fraction
#'
#' @description Draw a polygon representing a clone, annotated with cluster
#' label and cellular fraction
#'
#' @usage draw.clone(x, y, wid=1, len=1, col='gray', label=NA, cell.frac=NA)
#'
#' @param x: x coordinate
#' @param y: y coordinate
#' @param wid: width of the polygon (representing cellular fraction)
#' @param len: length of the polygon
#' @param col: fill color of the polygon
#' @param label: name of the clone
#' @param cell.frac: cellular fraction of the clone
#' @param cell.frac.position: position for cell.frac =
#' c('top.left', 'top.right', 'top.mid',
#' 'right.mid', 'right.top', 'right.bottom',
#' 'side', 'top.out')
#'
draw.clone <- function(x, y, wid=1, len=1, col='gray', label=NA, cell.frac=NA,
                       cell.frac.position='top.left', cell.frac.angle=NULL,
                       cell.frac.side.arrow=TRUE,
                       cell.frac.side.arrow.col='black',
                       text.size=1){
    beta = min(wid/5, (wid+len)/20)
    gamma = wid/2
    xx = c(x, x+beta, x+len, x+len, x+beta)
    yy = c(y, y+gamma, y+gamma, y-gamma, y-gamma)

    polygon(xx, yy, border='black', col=col, lwd=0.2)
    if (!is.na(label)){
        text(x+0.2*text.size, y, label, cex=text.size, adj=c(0,0.5))
    }
    if (!is.na(cell.frac)){
        cell.frac.x = 0
        cell.frac.y = 0
        angle = 0
        adj = c(0,0)
        if (cell.frac.position == 'top.left'){
            cell.frac.x = max(x+beta, x + 0.4)
            cell.frac.y = y+gamma#-0.3*text.size
            adj = c(0, 1)
        }else if (cell.frac.position == 'top.right'){
            cell.frac.x = x+len
            cell.frac.y = y+gamma
            adj = c(1, 1)
        }else if (cell.frac.position == 'top.mid'){
            cell.frac.x = x+beta+(len-beta)/2
            cell.frac.y = y+gamma
            adj = c(0.5, 1)
        }else if (cell.frac.position == 'right.mid'){
            cell.frac.x = x+len
            cell.frac.y = y
            adj = c(0, 0.5)
            angle = 45
        }else if (cell.frac.position == 'right.top'){
            cell.frac.x = x+len
            cell.frac.y = y+gamma
            adj = c(0, 1)
            angle = 45
        }else if (cell.frac.position == 'side'){
            angle = atan(gamma/beta)*(180/pi)# - 5
            cell.frac.x = x+beta/3+0.3*text.size
            cell.frac.y = y+gamma/2-0.3*text.size
            adj = c(0.5, 0.5)
        }else if (cell.frac.position == 'top.out'){
            cell.frac.x = x+len
            cell.frac.y = y.out
            adj = c(1, 0.5)
            # increase y.out so next time, text will be plotted a little higher
            # to prevent overwritten! Also, x.out.shift is distance from arrow
            # to polygon
            y.out <<- y.out + 0.75
            x.out.shift <<- x.out.shift + 0.1
        }
        #debug
        #cat('x=', cell.frac.x, 'y=', cell.frac.y, '\n')
        if (is.null(cell.frac.angle)){
            cell.frac.angle = angle
        }
        text(cell.frac.x, cell.frac.y, cell.frac,
             cex=text.size*0.7, srt=cell.frac.angle, adj=adj)
        if(cell.frac.side.arrow && cell.frac.position=='top.out'){
            # draw arrow
            x0 = cell.frac.x
            y0 = cell.frac.y
            x1 = cell.frac.x+x.out.shift
            y1 = y0
            x2 = x2 = x1
            y2 = y+gamma
            x3 = x0
            y3 = y2
            segments(x0, y0, x1, y1, col=cell.frac.side.arrow.col, lwd=1.5)
            segments(x1, y1, x2, y2, col=cell.frac.side.arrow.col, lwd=1.5)
            arrows(x0=x2, y0=y2, x1=x3, y1=y3,col=cell.frac.side.arrow.col,
                   length=0.025, lwd=1.5)
        }
    }
}

#' Rescale VAF of subclones s.t. total VAF must not exceed parent clone VAF
#'
#' @description Rescale VAF of subclones s.t. total VAF must not exceed parent
#' clone VAF. When infered using statistical test, sometime the estimated
#' total mean/median VAFs of subclones > VAF of parent clones which makes
#' drawing difficult (ie. subclone receive wider polygon than parent clone.
#' This function rescale the VAF of the subclone for drawing purpose only,
#' not for the VAF estimate.
#'
#' @param v: clonal structure data frame as output of enumerate.clones
#'
#'
rescale.vaf <- function(v, down.scale=0.99){
    #v = vx
    #print(v)
    rescale <- function(i){
        #print(i)
        parent = v[i,]$lab
        parent.vaf = v[i, ]$vaf
        subclones.idx = which(v$parent == parent)
        sum.sub.vaf = sum(v[subclones.idx,]$vaf)
        scale = ifelse(sum.sub.vaf > 0, parent.vaf/sum.sub.vaf, 1)
        #debug
        #cat('parent.vaf=', parent.vaf, ';sum.sub.vaf=', sum.sub.vaf,
        #    ';scale=', scale, '\n')
        for (idx in subclones.idx){
            if(scale < 1){
                #cat('Scaling...\n')
                #print(v)
                v[idx,]$vaf <<- down.scale*scale*v[idx,]$vaf
                #print(v)
            }
        }
        for (idx in subclones.idx){
            rescale(idx)
        }

    }
    root.idx = which(v$parent == -1)
    rescale(root.idx)
    return(v)
}


#' Set vertical position of clones within a sample clonal polygon visualization
#'
#' @description All subclones of a clone will be positioned next to each
#' other from the bottom of the polygon of the parent clone.
#' this function set the y.shift position depending on VAF
#'
#' @param v: clonal structure data frame as output of enumerate.clones
#'
#'
set.position <- function(v){
    v$y.shift = 0
    max.vaf = max(v$vaf)
    scale = 0.5/max.vaf
    #debug
    #print(v)
    for (i in 1:nrow(v)){
        vi = v[i,]
        subs = v[!is.na(v$parent) & v$parent == vi$lab,]
        if (nrow(subs) == 0){next}
        vafs = subs$vaf
        margin = (vi$vaf - sum(vafs))/length(vafs)*scale
        spaces = rep(0, length(vafs))
        if (length(spaces) >= 2){
            for (j in 2:length(spaces)){
                spaces[j] = sum(vafs[1:j-1]+margin)
            }
        }else{
            # re-centering if only 1 subclone inside another
            spaces = (vi$vaf-vafs)/2
        }
        #debug
        #print(subs)
        v[!is.na(v$parent) & v$parent == vi$lab,]$y.shift = spaces
    }
    #print(v)
    return(v)
}

get.cell.frac.ci <- function(vi, include.p.value=T, sep=' - '){
    cell.frac.lower = ifelse(vi$free.lower <= 0, '0',
                             gsub('\\.[0]+$|0+$', '',
                                  sprintf('%0.2f', 200*vi$free.lower)))
    cell.frac.upper = ifelse(vi$free.upper >= 0.5, '100%',
                             gsub('\\.[0]+$|0+$', '',
                                  sprintf('%0.2f%%', 200*vi$free.upper)))
    cell.frac = paste0(cell.frac.lower, sep , cell.frac.upper)
    if(include.p.value){
        cell.frac = paste0(cell.frac, '(',sprintf('%0.2f',
                            vi$free.confident.level.non.negative),
                            #',p=', 1-vi$p.value,
                           ')')
    }

    #debug
    #print(vi$free.confident.level.non.negative)
    #if (vi$free.confident.level.non.negative == 0.687){
    #    print(vi$free.confident.level.non.negative)
    #    print(cell.frac)
    #    print(vi)
    #    vii <<- vi
    #}

    return(cell.frac)
}

#' Draw clonal structures/evolution of a single sample
#'
#' @description Draw clonal structure for a sample with single or multiple
#' clones/subclones using polygon and tree plots
#'
#' @param v: clonal structure data frame (output of enumerate.clones)
#'
draw.sample.clones <- function(v, x=2, y=0, wid=30, len=8,
                               label=NULL, text.size=1,
                               cell.frac.ci=F,
                               top.title=NULL){
    #print(v)
    v = rescale.vaf(v)
    # scale VAF so that set.position works properly, and drawing works properly
    max.vaf = max(v$vaf)
    scale = 0.5/max.vaf
    v$vaf = v$vaf*scale
    max.vaf = max(v$vaf)
    high.vaf = max.vaf - 0.02
    low.vaf = 0.2
    y.out <<- wid*max.vaf/2+0.5
    x.out.shift <<- 0.1

    #print(v)

    draw.sample.clone <- function(i){
        vi = v[i,]
        #debug
        #cat('drawing', vi$lab, '\n')
        if (vi$vaf > 0){
            #if (vi$parent == 0){# root
            #if (is.na(vi$parent)){
            if (!is.na(vi$parent) && vi$parent == -1){
                xi = x
                yi = y
                leni = len
            }else{
                x.shift = 1
                if (vi$y.shift + vi$vaf >= high.vaf && vi$vaf < low.vaf){
                    x.shift = 2
                }
                par = v[v$lab == vi$parent,]
                xi = par$x + x.shift
                yi = par$y - wid*par$vaf/2 + wid*vi$vaf/2 + vi$y.shift*wid
                leni = par$len - x.shift
            }
            #cell.frac.position = ifelse(vi$free.lower < 0.05 & vi$vaf > 0.25, 'side', 'top.right')
            #cell.frac.position = ifelse(vi$free.lower < 0.05, 'top.out', 'top.right')
            cell.frac.position = ifelse(vi$free < 0.05, 'top.out', 'top.right')
            #cell.frac.position = ifelse(vi$num.subclones > 0 , 'right.top', 'right.mid')
            #cell.frac.position = 'top.mid'
            cell.frac = paste0(gsub('\\.[0]+$|0+$', '',
                                    sprintf('%0.2f', vi$free.mean*2*100)), '%')
            if(cell.frac.ci){
                cell.frac = get.cell.frac.ci(vi, include.p.value=T)
            }
            draw.clone(xi, yi, wid=wid*vi$vaf, len=leni, col=vi$color,
                       label=vi$lab,
                       cell.frac=cell.frac,
                       cell.frac.position=cell.frac.position,
                       cell.frac.side.arrow.col=vi$color,
                       text.size=text.size)
            v[i,]$x <<- xi
            v[i,]$y <<- yi
            v[i,]$len <<- leni
            for (j in 1:nrow(v)){
                #cat('---', v[j,]$parent,'\n')
                if (!is.na(v[j,]$parent) && v[j,]$parent != -1 &&
                        v[j,]$parent == vi$lab){
                    draw.sample.clone(j)
                }
            }
        }
    }
    plot(c(1, 10),c(-10,10), type = "n", xlab='', ylab='', xaxt='n',
         yaxt='n', axes=F)
    if (!is.null(label)){
        text(x-1, y, label=label, srt=90, cex=text.size, adj=c(0.5,1))
    }
    if (!is.null(top.title)){
        text(x, y+10, label=top.title, cex=(text.size*1.25), adj=c(0,0.5))
    }

    # move root to the first row and plot
    root = v[!is.na(v$parent) & v$parent == -1,]
    v = v[is.na(v$parent) | v$parent != -1,]
    v = rbind(root, v)
    v = set.position(v)
    v$x = 0
    v$y = 0
    v$len = 0

    #debug
    #print(v)

    draw.sample.clone(1)
}

#' Construct igraph object from clonal structures of a sample
#'
make.graph <- function(v, cell.frac.ci=T){
    library(igraph)
    #v = v[!is.na(v$parent),]
    #v = v[!is.na(v$parent) | v$vaf != 0,]
    v = v[!is.na(v$parent),]
    #rownames(v) = seq(1,nrow(v))
    rownames(v) = v$lab
    g = matrix(0, nrow=nrow(v), ncol=nrow(v))
    rownames(g) = rownames(v)
    colnames(g) = rownames(v)
    for (i in 1:nrow(v)){
        par.lab = v[i,]$parent
        if (!is.na(par.lab) && par.lab != -1){
            #if(par.lab != 0){
            par.idx = rownames(v[v$lab == par.lab,])
            #debug
            #cat(par.lab, '--', par.idx, '\n')
            g[par.idx, i] = 1
            #}
        }
    }
    #print(g)
    g <- graph.adjacency(g)
    cell.frac = gsub('\\.[0]+$|0+$', '', sprintf('%0.2f%%', v$free.mean*2*100))
    if(cell.frac.ci){
        cell.frac = get.cell.frac.ci(v, include.p.value=F, sep=' -\n')
    }
    labels = paste0(v$lab,'\n', cell.frac)
    V(g)$name = labels
    V(g)$color = v$color
    return(list(graph=g, v=v))
}

#' Draw all enumerated clonal models for a single sample
#' @param x: output from enumerate.clones()
draw.sample.clones.all <- function(x, outPrefix, object.to.plot='polygon',
                                   ignore.clusters=NULL){
    pdf(paste0(outPrefix, '.pdf'), width=6, height=6)
    for(i in 1:length(x)){
        xi = x[[i]]
        #xi = scale.cell.frac(xi, ignore.clusters=ignore.clusters)
        if (object.to.plot == 'polygon'){
            draw.sample.clones(xi, cell.frac.ci=T)
        }else{
            plot.tree(xi, node.shape='circle', node.size=35, cell.frac.ci=T)
        }
    }
    dev.off()
    cat(outPrefix, '\n')
}


#' Plot clonal evolution tree
#'
#' @description Plot a tree representing the clonal structure of a sample
#' Return the graph object (with node name prefixed with node.prefix.to.add)
#'
#' @param v: clonal structure data frame as the output of enumerate.clones
#' @param display: tree or graph
#'
plot.tree <- function(v, node.shape='square', display='tree',
                      node.size=50,
                      tree.node.text.size=1,
                      cell.frac.ci=T,
                      node.prefix.to.add=NULL,
                      out.prefix=NULL, out.format='graphml'){
    library(igraph)
    x = make.graph(v, cell.frac.ci=cell.frac.ci)
    #print(v)
    g = x$graph
    v = x$v
    root.idx = which(!is.na(v$parent) & v$parent == '-1')
    #cell.frac = gsub('\\.[0]+$|0+$', '', sprintf('%0.2f%%', v$free*2*100))


    #V(g)$color = v$color
    #display = 'graph'
    if(display == 'tree'){
        layout = layout.reingold.tilford(g, root=root.idx)
    }else{
        layout = NULL
    }

    plot(g, edge.color='black', layout=layout,
         edge.arrow.size=0.75, edge.arrow.width=0.75,
         vertex.shape=node.shape, vertex.size=node.size,
         vertex.label.cex=tree.node.text.size)#, vertex.color=v$color,
    #vertex.label=labels)

    # remove newline char because Cytoscape does not support multi-line label
    V(g)$name = gsub('\n', ' ', V(g)$name, fixed=T)
    if (!is.null(node.prefix.to.add)){
        V(g)$name = paste0(node.prefix.to.add, V(g)$name)
    }

    if (!is.null(out.prefix)){
        out.file = paste0(out.prefix, '.', out.format)
        #cat('Writing tree to ', out.file, '\n')
        write.graph(g, file=out.file, format=out.format)
    }

    return(g)

}

#' Write clonal evolution tree to file
#' TODO: do!
write.tree <- function(v, out.file, out.format='tabular'){
    v = v[, c('parent', 'lab', '')]
}

get.model.score <- function(v){
    return(prod(v$p.value[!is.na(v$p.value)]))
}

#' Find matched models between samples
#' infer clonal evolution models, given all evolve from the 1st sample
find.matched.models <- function(vv, samples){
    nSamples = length(samples)
    matched = NULL
    scores = NULL
    find.next.match <- function(prim.idx, prim.model.idx,
                                met.idx, met.model.idx,
                                matched.models, model.scores){
        if (met.idx > nSamples){
            if (all(matched.models > 0)){
                matched <<- rbind(matched, matched.models)
                scores <<- rbind(scores, model.scores)
            }
        }else{
            for (j in 1:length(matched.models)){
                match.with.all.models = TRUE
                if (matched.models[j] > 0){
                    if(!match.sample.clones(vv[[j]][[matched.models[j]]],
                                            vv[[met.idx]][[met.model.idx]])){
                        match.with.all.models = FALSE
                        break
                    }
                }
            }
            #debug
            #cat(paste(matched.models[matched.models>0], collapse='-'),
            #          met.model.idx, match.with.all.models, '\n')
            if (match.with.all.models){
                matched.models[[met.idx]] = met.model.idx
                model.scores[[met.idx]] =
                    get.model.score(vv[[met.idx]][[met.model.idx]])
                find.next.match(prim.idx, prim.model.idx, met.idx+1, 1,
                                matched.models, model.scores)
            }#else{
            if (length(vv[[met.idx]]) > met.model.idx){
                matched.models[[met.idx]] = 0
                model.scores[[met.idx]] = 0
                find.next.match(prim.idx, prim.model.idx,
                                met.idx, met.model.idx+1,
                                matched.models, model.scores)
            }
            #find.next.match(prim.idx)
            #}

        }
    }
    for (prim.model in 1:length(vv[[1]])){
        # cat('prim.model =', prim.model, '\n')
        matched.models = c(prim.model,rep(0,nSamples-1))
        model.scores = c(get.model.score(vv[[1]][[prim.model]]),
                         rep(0,nSamples-1))
        find.next.match(1, prim.model, 2, 1, matched.models, model.scores)
    }
    return(list(matched.models=matched, scores=scores))
}



#' Infer clonal structures and evolution models for multiple samples
#'
#' @description Infer clonal structures and evolution models for multi cancer
#' samples from a single patients (eg. primary tumors, metastatic tumors,
#' xenograft tumors, multi-region samples, etc.)
#'
#' @param c: clonality analysis data frame, consisting of N+1 columns. The first
#' column must be named 'cluster' and hold variant cluster number (ie. use number
#' to name cluster, starting from 1,2,3... 0 is reserved for normal cell clone).
#' The next N columns contain VAF estimated for the corresponding cluster
#' (values range from 0 to 0.5)
#' @param variants: data frame of the variants
#' @param cluster.col.name: column that holds the cluster identity
#' @param founding.cluster: the cluster of variants that are found in all
#' samples and is beleived the be the founding events. This is often
#' the cluster with highest VAF and most number of variants
#' @param ignore.clusters: clusters to ignore (not inluded in the models). Those
#' are the clusters that are thought of as outliers, artifacts, etc. resulted
#' from the error or bias of the sequencing and analysis. This is provided as
#' a debugging tool
#' @param model: cancer evolution model to use, c('monoclonal', 'polyclonal').
#' monoclonal model assumes the orginal tumor (eg. primary tumor) arises from
#' a single normal cell; polyclonal model assumes the original tumor can arise
#' from multiple cells (ie. multiple founding clones). In the polyclonal model,
#' the total VAF of the separate founding clones must not exceed 0.5
#' @param subclonal.test.method: 'bootstrap' = perform bootstrap subclonal test
#' 'none' = straight comparison of already estimated VAF for each cluster
#' provided in c
#' @param subclonal.test.model: What model to use when generating the bootstrap
#' Values are: c('non-parametric', 'normal', 'normal-truncated', 'beta',
#' 'beta-binomial')
#'
#'
infer.clonal.models <- function(c=NULL, variants=NULL,
                                cluster.col.name='cluster',
                                founding.cluster=NULL,
                                ignore.clusters=NULL,
                                vaf.col.names=NULL,
                                vaf.in.percent=TRUE,
                                sample.names=NULL,
                                model='monoclonal',
                                subclonal.test='none',
                                subclonal.test.model='non-parametric',
                                boot=NULL,
                                num.boots=1000,
                                p.value.cutoff=0.05,
                                alpha=0.05,
                                min.cluster.vaf=0,
                                verbose=TRUE){
    if (is.null(vaf.col.names)){
        # check format of input, find vaf column names
        if(!is.null(c)){
            cnames = colnames(c)
        }else if(!is.null(variants)){
            cnames = colnames(variants)
        }else{
            stop('ERROR: Need at least parameter c or variants\n')
        }
        if (!(cluster.col.name %in% cnames && length(cnames) >= 2)){
            stop('ERROR: No cluster column and/or no sample\n')
        }
        vaf.col.names = setdiff(cnames, cluster.col.name)
    }

    # convert cluster column to character
    if (!is.null(c)) {
        c[[cluster.col.name]] = as.character(c[[cluster.col.name]])
    }
    if (!is.null(variants)){
        variants[[cluster.col.name]] = as.character(variants[[cluster.col.name]])
    }


    if (is.null(sample.names)){
        sample.names = vaf.col.names
    }

    nSamples = length(sample.names)
    n.vaf.cols = length(vaf.col.names)

    if (nSamples != n.vaf.cols || nSamples == 0){
        stop('ERROR: sample.names and vaf.col.names have different length
         or both have zero length!\n')
    }

    if (nSamples >= 1 && verbose){
        for (i in 1:nSamples){
            cat('Sample ', i, ': ', sample.names[i], ' <-- ',
                vaf.col.names[i], '\n', sep='')
        }
    }

    # if polyclonal model, add normal as founding clone
    add.normal = NA
    if (model == 'monoclonal'){
        add.normal = FALSE
    }else if (model == 'polyclonal'){
        add.normal = TRUE
    }
    if (is.na(add.normal)){
        stop(paste0('ERROR: Model ', model, ' not supported!\n'))
    }
    if(verbose){cat('Using ', model, ' model\n', sep='')}

    # prepare cluster data and infer clonal models for individual samples
    if (is.null(c)){
        c = estimate.clone.vaf(variants, cluster.col.name,
                               vaf.col.names, vaf.in.percent=vaf.in.percent)
    }
    vv = list()
    for (i in 1:nSamples){
        s = vaf.col.names[i]
        v = make.clonal.data.frame(c[[s]], c[[cluster.col.name]], add.normal)
        if (subclonal.test == 'none'){
            #models = enumerate.clones.absolute(v)
            models = enumerate.clones(v, sample=s,
                                      founding.cluster=founding.cluster,
                                      ignore.clusters=ignore.clusters)
        }else if (subclonal.test == 'bootstrap'){
            if (is.null(boot)){
                #boot = generate.boot(variants, vaf.col.names=vaf.col.names,
                #                     vaf.in.percent=vaf.in.percent,
                #                     num.boots=num.boots)

                boot = generate.boot(variants, vaf.col.names=vaf.col.names,
                                     vaf.in.percent=vaf.in.percent,
                                     num.boots=num.boots,
                                     bootstrap.model=subclonal.test.model)
                bbb <<- boot
            }
            models = enumerate.clones(v, sample=s, variants, boot=boot,
                                      founding.cluster=founding.cluster,
                                      ignore.clusters=ignore.clusters,
                                      min.cluster.vaf=min.cluster.vaf,
                                      p.value.cutoff=p.value.cutoff,
                                      alpha=alpha)
        }

        if(verbose){cat(s, ':', length(models),
                        'clonal architecture model(s) found\n')}
        if (length(models) == 0){
            print(v)
            stop(paste('ERROR: No clonal models for sample:', s,
                       '\nCheck data or remove this sample, then re-run.
                       \nAlso check if founding.cluster was set correctly!'))
        }else{
            vv[[s]] = models
        }
    }

    # infer clonal evolution models,  an accepted model must satisfy
    # the clone-subclonal relationship
    matched = NULL
    scores = NULL
    if (nSamples == 1 && length(vv[[1]]) > 0){
        num.models = length(vv[[1]])
        matched = data.frame(x=1:num.models)
        colnames(matched) = c(sample.names[1])
        scores = data.frame(x=rep(0,num.models))
        colnames(scores) = c(sample.names[1])
        for (i in 1:num.models){
            scores[i,1] = get.model.score(vv[[1]][[i]])
        }
        scores$model.score = scores[, 1]
    }
    if (nSamples >= 2){
        z = find.matched.models(vv, sample.names)
        matched = z$matched.models
        scores = z$scores
        if (!is.null(matched)){
            rownames(matched) = seq(1,nrow(matched))
            colnames(matched) = sample.names
            matched = as.data.frame(matched)
            rownames(scores) = seq(1,nrow(matched))
            colnames(scores) = sample.names
            scores = as.data.frame(scores)
            scores$model.score = apply(scores, 1, prod)
        }
    }
    if (!is.null(matched)){
        # sort by score
        matched = matched[order(scores$model.score, decreasing=T), ,drop=F]
        scores = scores[order(scores$model.score, decreasing=T), , drop=F]
    }
    num.matched.models = ifelse(is.null(matched), 0, nrow(matched))
    if (verbose){ cat(paste0('Found ', num.matched.models,
                             ' compatible evolutional models\n'))}
    return (list(models=vv, matched=list(index=matched, scores=scores)))

}


#' Scale cellular fraction in a clonal architecture model
#'
#' @description Max VAF will be determined, and all vaf will be scaled such that
#' max VAF will be 0.5
#'
scale.cell.frac <- function(m, ignore.clusters=NULL){
    max.vaf = max(m$vaf[!(m$lab %in% as.character(ignore.clusters))])
    scale = 0.5/max.vaf
    m$vaf = m$vaf*scale
    m$free = m$free*scale
    m$free.mean = m$free.mean*scale
    m$free.lower = m$free.lower*scale
    m$free.upper = m$free.upper*scale
    m$occupied = m$occupied*scale
    return(m)
}


#' Plot evolution models (polygon plots and trees) for multi samples
#'
#' @description Plot evolution models inferred by infer.clonal.models function
#' Two types of plots are supported: polygon plot and tree plot
#'
#' @param models: list of model output from infer.clonal.models function
#' @param out.dir: output directory for the plots
#' @param matched: data frame of compatible models for multiple samples
#' @param scale.monoclonal.cell.frac: c(TRUE, FALSE); if TRUE, scale cellular
#' fraction in the plots (ie. cell fraction will be scaled by 1/purity =
#' 1/max(VAF*2))
#' @param width: width of the plots (in), if NULL, automatically choose width
#' @param height: height of the plots (in), if NULL, automatically choose height
#' @param out.format: format of the plot files ('png', 'pdf', 'pdf.multi.files')
#' @param resolution: resolution of the PNG plot file
#' @param overwrite.output: if TRUE, overwrite output directory, default=FALSE
#' @param max.num.models.to.plot: max number of models to plot; default = 10
#' @param ignore.clusters: cluster to ignore in VAF scaling, should be the ones
#' that are ignored in infer.clonal.models (TODO: automatically identify to what
#' cluster should the VAF be scaled from a model, and remove this param)
#' if NULL, unlimited
#'
plot.clonal.models <- function(models, out.dir, matched=NULL,
                               variants=NULL,
                               box.plot=FALSE,
                               box.plot.text.size=1.5,
                               cluster.col.name = 'cluster',
                               scale.monoclonal.cell.frac=TRUE,
                               ignore.clusters=NULL,
                               width=NULL, height=NULL, text.size=1,
                               tree.node.shape='circle',
                               tree.node.size = 50,
                               tree.node.text.size=1,
                               out.format='png', resolution=300,
                               overwrite.output=FALSE,
                               max.num.models.to.plot=10,
                               cell.frac.ci=FALSE,
                               show.score=T,
                               out.prefix='model')
{
    if (!file.exists(out.dir)){
        dir.create(out.dir)
    }else{
        if (!overwrite.output){
            stop(paste('ERROR: Output directory (', out.dir,
                       ') exists. Quit!\n'))
        }
    }
    nSamples = length(models)
    samples = names(models)
    w = ifelse(is.null(width), 7, width)
    h = ifelse(is.null(height), 3*nSamples, height)
    if(box.plot && is.null(variants)){
        box.plot = F
        message('box.plot = TRUE, but variants = NULL. No box plot!')
    }
    if (!is.null(matched$index)){
        scores = matched$scores
        matched = matched$index
        num.models = nrow(matched)
        if (num.models > max.num.models.to.plot &&
                !is.null(max.num.models.to.plot)){
            message(paste0(num.models,
               ' models requested to plot. Only plot the first ',
               max.num.models.to.plot,
               ' models. \nChange "max.num.models.to.plot" to plot more.\n'))
            matched = head(matched, n=max.num.models.to.plot)
        }
        if (out.format == 'pdf'){
            pdf(paste0(out.dir, '/', out.prefix, '.pdf'), width=w, height=h)
        }

        for (i in 1:nrow(matched)){
            combined.graph = NULL
            this.out.prefix = paste0(out.dir, '/', out.prefix, '-', i)
            if (out.format == 'png'){
                png(paste0(this.out.prefix, '.png'), width=w,
                    height=h, res=resolution, units='in')
            }else if (out.format == 'pdf.multi.files'){
                pdf(paste0(this.out.prefix, '.pdf'), width=w, height=h)
            }else if (out.format != 'pdf'){
                stop(paste0('ERROR: output format (', out.format,
                            ') not supported.\n'))
            }

            num.plot.cols = ifelse(box.plot, 3, 2)
            par(mfrow=c(nSamples,num.plot.cols), mar=c(0,0,0,0))
            mat = t(matrix(seq(1, nSamples*num.plot.cols), ncol=nSamples))
            #print(mat)
            ww = rep(1, num.plot.cols)
            ww[length(ww)] = 1
            if (box.plot){
                ww[1] = 1
            }
            hh = rep(1, nSamples)
            layout(mat, ww, hh)

            for (k in 1:length(samples)){
                s = samples[k]
                s.match.idx = matched[[s]][i]
                m = models[[s]][[matched[[s]][i]]]
                if (scale.monoclonal.cell.frac){
                    m = scale.cell.frac(m, ignore.clusters=ignore.clusters)
                }
                lab = s
                # turn this on to keep track of what model matched
                lab = paste0(s, ' (', s.match.idx, ')')
                if (show.score){
                    lab = paste0(s, '\n(prob=',
                                 sprintf('%0.3f', scores[[s]][i]), ')')
                }
                top.title = NULL
                if (k == 1 && show.score){
                    top.title = paste0('Model prob = ', scores$model.score[i])
                }
                if (box.plot){
                    current.mar = par()$mar
                    par(mar=c(3,5,3,3))
                    with(variants, boxplot(get(s) ~ get(cluster.col.name),
                                           cex.lab=box.plot.text.size,
                                           cex.axis=box.plot.text.size,
                                           cex.main=box.plot.text.size,
                                           cex.sub=box.plot.text.size,
                                           ylab=s))
                    par(mar=current.mar)
                }
                draw.sample.clones(m, x=2, y=0, wid=30, len=7,
                                   label=lab,
                                   text.size=text.size,
                                   cell.frac.ci=cell.frac.ci,
                                   top.title=top.title)

                gs = plot.tree(m, node.shape=tree.node.shape,
                               node.size=tree.node.size,
                               tree.node.text.size=tree.node.text.size,
                               cell.frac.ci=cell.frac.ci,
                               node.prefix.to.add=paste0(s,': '),
                               out.prefix=paste0(this.out.prefix, '__', s))
                if (is.null(combined.graph)){
                    combined.graph = gs
                }else{
                    combined.graph = graph.union(combined.graph, gs,
                                                 byname=TRUE)
                    # set color for all clones, if missing in 1st sample
                    # get color in other sample
                    V(combined.graph)$color <-
                        ifelse(is.na(V(combined.graph)$color_1),
                               V(combined.graph)$color_2,
                               V(combined.graph)$color_1)
                }
            }
            if (out.format == 'png' || out.format == 'pdf.multi.files'){
                dev.off()
            }
            write.graph(combined.graph,
                        file=paste0(this.out.prefix, '.graphml'),
                        format='graphml')
        }
        if (out.format == 'pdf'){
            dev.off()
        }

    }else{# plot all
        # TODO: plot all models for all samples separately.
        # This will serve as a debug tool for end-user when their models
        # from different samples do not match.
        message('No compatible multi-sample models provided.
                Individual sample models will be plotted!n')
        for (s in names(models)){
            draw.sample.clones.all(models[[s]],
                                   paste0(out.dir, '/', out.prefix, '-', s))
        }

    }
    cat(paste0('Output plots are in: ', out.dir, '\n'))
}


# Extra functionality:
# - estimate VAF from clusters of variants

#' Estimate VAFs of clones/clusters from clonality analysis result
#'
#' @description Estimate VAFs of clones/clusters by calculating the median of
#' the VAFs of all variants provided.
#'
#' @param v: data frame containing variants' VAF and clustering
#' @param cluster.col.name: name of the column that hold the cluster ID
#' @param vaf.col.names: names of the columns that hold VAFs
#' @param method: median or mean
#'
estimate.clone.vaf <- function(v, cluster.col.name='cluster',
                               vaf.col.names=NULL,
                               vaf.in.percent=TRUE,
                               method='median',
                               ref.count.col.names=NULL,
                               var.count.col.names=NULL,
                               depth.col.names=NULL){
    clusters = sort(unique(v[[cluster.col.name]]))
    clone.vafs = NULL

    if (is.null(vaf.col.names)){
        vaf.col.names = setdiff(colnames(v), cluster.col.name)
    }

    for (cl in clusters){
        #cat('cluster: ', cl, '\n')
        if (method == 'median'){
            median.vafs = apply(v[v[[cluster.col.name]]==cl,vaf.col.names],
                            2, median)
        }else if (method == 'mean'){
            median.vafs = apply(v[v[[cluster.col.name]]==cl,vaf.col.names],
                                2, mean)
        }
        median.vafs = as.data.frame(t(median.vafs))
        #median.vafs[[cluster.col.name]] = cl
        if (is.null(clone.vafs)){
            clone.vafs = median.vafs
        }else{
            clone.vafs = rbind(clone.vafs, median.vafs)
        }
    }
    clone.vafs = cbind(clusters, clone.vafs)
    colnames(clone.vafs)[1] = cluster.col.name
    clone.vafs = clone.vafs[order(clone.vafs[[cluster.col.name]]),]
    if (vaf.in.percent){
        clone.vafs[,vaf.col.names] = clone.vafs[,vaf.col.names]/100.00
    }
    return(clone.vafs)
}

#' Adjust clone VAF according to significant different test result
#' If two clones have close VAF, adjust the smaller VAF to the bigger
#' TODO: this test does not work yet, has to think more carefully about what
#' test to use, as well as test involving multiple samples
adjust.clone.vaf <- function(clone.vafs, var, cluster.col.name,
                             founding.cluster=1,
                             adjust.to.founding.cluster.only=TRUE,
                             p.value.cut=0.01){
    vaf.names = colnames(clone.vafs[2:length(colnames(clone.vafs))])
    founding.cluster.idx = which(clone.vafs$cluster == founding.cluster)
    base.clusters.idx = unique(c(founding.cluster.idx, 1:(nrow(clone.vafs)-1)))
    if (adjust.to.founding.cluster.only){
        base.clusters.idx = founding.cluster.idx
    }
    print(base.clusters.idx)
    for (vaf.name in vaf.names){
        #for (i in 1:(nrow(clone.vafs)-1)){
        for (i in base.clusters.idx){
            ci = clone.vafs$cluster[i]
            vaf.i = clone.vafs[clone.vafs[[cluster.col.name]]==ci, vaf.name]
            for (j in (i+1):nrow(clone.vafs)){
                cj = clone.vafs$cluster[j]
                vaf.j = clone.vafs[clone.vafs[[cluster.col.name]]==cj, vaf.name]
                if (!clone.vaf.diff(var[var[[cluster.col.name]]==ci,vaf.name],
                                    var[var[[cluster.col.name]]==cj,vaf.name])){

                    clone.vafs[clone.vafs[[cluster.col.name]]==cj, vaf.name] =
                        vaf.i
                }
            }
        }
    }
    return(clone.vafs)
}

#' Test if two clones have different VAFs
#' Deprecated!
#'
#' @description Test if the two clones/clusters have different VAFs, using a
#' Mann-Whitney U test
#'
#' @param clone1.vafs: VAFs of all variants in clone/cluster 1
#' @param clone2.vafs: VAFs of all variants in clone/cluster 2
#' @param p.value.cut: significant level (if the test produce p.value
#' less than or equal to p.val => significantly different)
#'
clone.vaf.diff <- function(clone1.vafs, clone2.vafs, p.value.cut=0.05){
    tst = wilcox.test(clone1.vafs, clone2.vafs)
    #print(tst)
    if (is.na(tst$p.value) || tst$p.value <= p.value.cut){
        return(TRUE)
    }else{
        return(FALSE)
    }
}


#a6cee3 light blue
#b2df8a light green
#cab2d6 light purple
#fdbf6f light orange
#fb9a99 pink/brown
#d9d9d9 light gray
#999999 gray
#33a02c green
#ff7f00 orange
#1f78b4 blue
#fca27e salmon/pink
#ffffb3 light yellow
#fccde5 light purple pink
#fb8072 light red
#b3de69 light green
#f0ecd7 light light brown/green
#e5f5f9 light light blue
#' Get the hex string of the preset colors optimized for plotting both
#' polygon plots and mutation scatter plots, etc.
get.clonevol.colors <- function(num.colors){
    colors = c('#a6cee3', '#b2df8a', '#cab2d6', '#fdbf6f', '#fb9a99',
               '#1f78b4','#999999', '#33a02c', '#ff7f00', '#bc80bd',
               '#fca27e', '#ffffb3', '#fccde5', '#fb8072', '#d9d9d9',
               '#f0ecd7', rep('#e5f5f9',100))
    if (num.colors > length(colors)){
        stop('ERROR: Not enough colors!\n')
    }else{
        return(colors[1:num.colors])
    }
}

plot.clonevol.colors <- function(num.colors=17){
    colors = get.clonevol.colors(num.colors)
    x = data.frame(hex=colors, val=1, stringsAsFactors=F)
    x$hex = paste0(sprintf('%02d', seq(1,num.colors)), '\n', x$hex)
    names(colors) = x$hex
    p = (ggplot(x, aes(x=hex, y=val, fill=hex))
         + geom_bar(stat='identity')
         + theme_bw()
         + scale_fill_manual(values=colors)
         + theme(legend.position='none')
         + ylab(NULL) + xlab(NULL)
         + theme(axis.text.y=element_blank())
         + theme(axis.ticks.y=element_blank())
         + ggtitle('Clonevol colors'))
    ggsave(p, file='clonevol.colors.pdf', width=num.colors*0.75, height=4)
}
