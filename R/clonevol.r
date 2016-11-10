# Clonevol: Inferring and visualizing clonal evolution in multi-sample cancer
# sequencing
# Created by: Ha Dang <hdangATgenomeDOTwustlDOTedu>
# Date: Dec. 25, 2014
# Last modified:
#   Dec. 27, 2014  -- more than 2 samples model inference and plotting
#   Feb. 02, 2015  -- polyclonal model supported, bugs fixed.
#   Mar. 07, 2015  -- subclonal bootstrap test implementation, bugs fixed.
#   Lots of other modifications (see github) --
#   Jun. 23, 2016  -- Clean up for initial release
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
                                    #normal.clone.color='#f0f0f0', # very light gray
                                    normal.clone.color='#e5f5f9', # very light blue
                                    founding.label=NULL, colors=NULL){
    v = data.frame(lab=as.character(labels), vaf=vafs, stringsAsFactors=F)
    if (is.null(colors)){
        #colors = c('#a6cee3', '#b2df8a', '#cab2d6', '#fdbf6f', '#fb9a99',
        #           '#d9d9d9','#999999', '#33a02c', '#ff7f00', '#1f78b4',
        #           '#fca27e', '#ffffb3', '#fccde5', '#fb8072', '#b3de69',
        #           'f0ecd7', rep('#e5f5f9',1))
        colors=get.clonevol.colors(nrow(v))
        # if normal clone added, set it color
        if (v$lab[1] == '0'){colors = c(normal.clone.color, colors)}
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


#' Calculate CI of CCF, pvals, etc.
#' @param vx: clonal evolution of a sample data frame
#' @param sample: name of vaf.col
#' @param i: row of vx where clone needs to be estimated
#' @param boot: pregenerated bootstrap
#' @param min.cluster.vaf: minimum cluster vaf to consider non-zero
#' @param alpha: alpha of CI
#' @param t: if subclonal.test was run already, provide the return object, if NULL
#' subclonal.test will be run
#' Return vx with additionally annotated columns related to CCF
estimate.ccf <- function(vx, sample, i, boot, min.cluster.vaf,
    alpha, t=NULL, sub.clusters=NULL){
    if (is.null(t)){
        t = subclonal.test(sample,
           as.character(vx[i,]$lab),
           sub.clusters=sub.clusters, boot=boot,
           cdf=vx,
           min.cluster.vaf=min.cluster.vaf,
           alpha=alpha)
    }
    vx$free.mean[i] = t$free.vaf.mean
    vx$free.lower[i] = t$free.vaf.lower
    vx$free.upper[i] = t$free.vaf.upper
    vx$p.value[i] = t$p.value
    vx$free.confident.level[i] =
        t$free.vaf.confident.level
    vx$free.confident.level.non.negative[i] =
        t$free.vaf.confident.level.non.negative
    return(vx)

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
# TODO: for sample with one clone, this returns NA as cell frac, fix it.
# TODO: this use a lazy recursive algorithm which is slow when clonal
# architecture is complex (eg. many subclones with low VAF). Improve.
enumerate.clones <- function(v, sample=NULL, variants=NULL,
                             founding.cluster = NULL,
                             ignore.clusters=NULL,
                             subclonal.test.method='bootstrap',
                             boot=NULL,
                             p.value.cutoff=0.05,
                             alpha=0.05,
                             min.cluster.vaf=0){
    cat(sample, ': Enumerating clonal architectures...\n')
    vv = list() # to hold list of output clonal models
    #cat('*********: p : ', p.value.cutoff, '\n')
    findParent <- function(v, i){
        #print(i)
        if (i > nrow(v)){
            #debug
            #print(v)
            
            v$is.zero = ifelse(v$free.lower >= 0, F, T)
            # determine subclone
            clone.stat = determine.subclone(v, v$lab[!is.na(v$parent)
                                            & v$parent == '-1'])

            v$is.subclone = clone.stat$is.sub[v$lab]
            v$is.founder = clone.stat$is.founder[v$lab]
            rownames(v) = v$lab
            vv <<- c(vv, list(v))
        }else{
            #print(head(v))
            vaf = v[i,]$vaf
            if (!is.na(v[i,]$parent) && v[i,]$parent == '-1'){# root
                vx = v
                # estimate CCF for root if it does not have subclones nested yet,
                # just in case there is no other clone to be nested
                if (vx$num.subclones[i] == 0){
                    vx = estimate.ccf(vx, sample, i, boot, min.cluster.vaf, alpha, t=NULL)
                }
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
                        # cat('Testing...', sample, '-', j, parent.cluster,
                          # 'sub clusters:', sub.clusters, '\n')
                        t = subclonal.test(sample, parent.cluster, sub.clusters,
                                           boot=boot,
                                           cdf=v,
                                           min.cluster.vaf=min.cluster.vaf,
                                           alpha=alpha)
                        # hdng: test direction changed to greater, so p = 1 - p, sign flipped to <
                        # if a clonal nesting do not violate sum rule, this is unresolvable
                        # so it will be recorded as a temporary solution, later, other samples
                        # come in, we may find one or a few models resolvable, then applying
                        # cross rule (matching between samples) will solve the model
                        if(t$p.value < 1 - p.value.cutoff){
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
                                #t = subclonal.test(sample,
                                #       as.character(vx[i,]$lab),
                                #       sub.clusters=NULL, boot=boot,
                                #       cdf=vx,
                                #       min.cluster.vaf=min.cluster.vaf,
                                #       alpha=alpha)
                                #vx$free.mean[i] = t$free.vaf.mean
                                #vx$free.lower[i] = t$free.vaf.lower
                                #vx$free.upper[i] = t$free.vaf.upper
                                #vx$p.value[i] = t$p.value
                                #vx$free.confident.level[i] =
                                #    t$free.vaf.confident.level
                                #vx$free.confident.level.non.negative[i] =
                                #    t$free.vaf.confident.level.non.negative
                                vx = estimate.ccf(vx, sample, i, boot,
                                            min.cluster.vaf, alpha, t=NULL)
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
    cat('Determining if cluster VAF is significantly positive...\n')
    if (is.null(min.cluster.vaf)){
        cat('No min.cluster.vaf provided. Using bootstrap test\n')
    }else{
        cat('Exluding clusters whose VAF < min.cluster.vaf=',
                min.cluster.vaf, '\n', sep='')
    }
    for (i in 1:nrow(v)){
        cl = as.character(v[i,]$lab)
        if (is.null(min.cluster.vaf)){
            # test if VAF of this cluster cl is > 0
            t = subclonal.test(sample, parent.cluster=cl, sub.clusters=NULL,
                           boot=boot, min.cluster.vaf=min.cluster.vaf,
                           alpha=alpha)
            # if not, exclude from analysis
            v[i,]$excluded = ifelse(t$p.value > p.value.cutoff, TRUE, FALSE)
        }else{
            # if the median/mean (estimated earlier) VAF < e, do
            # not consider this cluster in this sample
            v[i,]$excluded = ifelse(v[i,]$vaf < min.cluster.vaf, TRUE, FALSE)
        }

    }

    cat('Non-positive VAF clusters:',
        paste(v$lab[v$excluded], collapse=','), '\n')

    # also exlude clusters in the ignore.clusters list
    if (!is.null(ignore.clusters)){
        ignore.idx = v$lab %in% as.character(ignore.clusters)
        v$excluded[ignore.idx] = TRUE
        cat('User ignored clusters: ',
            paste(v$lab[ignore.idx], collapse=','), '\n')
    }
    #print(v)

    # if normal sample (0) is included, the normal sample
    # will be root (polyclonal model), otherwise find the
    # founding clone and place it first
    if (v[1,]$lab == 0 || v[1,]$lab == '0'){
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

#' Generate fill points for bell/polygon plots
generate.fill.points <- function(x, y, num.points=50){
    n = length(x)
    k = floor(n/2)
    #z = c(1,2,k,k+1,k+2,k+3,n)
    gen.points <- function(x1, x2, y11, y12, y21, y22, step=NULL){
        xmin = min(x1,x2)
        xmax = max(x1,x2)
        ymid = (y11 + y12)/2
        rx = c()
        ry = c()
        if (is.null(step)){step = (xmax-xmin)/10}
        xx = seq(xmin, xmax, step)
        xx = xx[-length(xx)]
        for (xi in xx){
            yi = y11 + (y21 - y11)/(xmax - xmin)*(xi-xmin)
            yr = runif(num.points, ymid-(yi-ymid), yi)
            xr = runif(num.points, xi, xi+step)
            rx = c(rx, xr)
            ry = c(ry, yr)
        }
        return(list(x=rx, y=ry))
    }

    t1 = gen.points(x[1], x[2], y[1], y[1], y[2], y[n])
    t2 = gen.points(x[2], x[k], y[2], y[k], y[k+2], y[k+3])
    return(list(x=c(t1$x, t2$x), y=c(t1$y, t2$y)))

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
#' @param shape: c("polygon", "triangle", "parabol")
#' @param wid: width of the polygon (representing cellular fraction)
#' @param len: length of the polygon
#' @param col: fill color of the polygon
#' @param label: name of the clone
#' @param cell.frac: cellular fraction of the clone
#' @param cell.frac.position: position for cell.frac =
#' c('top.left', 'top.right', 'top.mid',
#' 'right.mid', 'right.top', 'right.bottom',
#' 'side', 'top.out')
#' @param cell.frac.top.out.space: spacing between cell frac annotation when
#' annotating on top of the plot
#' @param cell.frac.side.arrow.width: width of the line and arrow pointing
#' to the top edge of the polygon from the cell frac annotation on top
#' @param variant.names: list of variants to highlight inside the polygon
#' @param border.color: color of the border
#' @param bell.curve.step: vertical distance between the end point of clone
#' bell curve, and its mid point (increase this will make the curve steeper,
#' set this equal to zero will give no curve; use case: sometimes bell of
#' subclone cannot fit parent clone bell, so decrease this will help fitting
draw.clone <- function(x, y, wid=1, len=1, col='gray',
                       clone.shape='bell',
                       bell.curve.step = 0.25,
                       label=NA, cell.frac=NA,
                       #cell.frac.position='top.out',
                       cell.frac.position='right.mid',
                       cell.frac.top.out.space = 0.75,
                       cell.frac.side.arrow.width=1.5,
                       cell.frac.angle=NULL,
                       cell.frac.side.arrow=TRUE,
                       cell.frac.side.arrow.col='black',
                       variant.names=NULL,
                       variant.color='blue',
                       variant.angle=NULL,
                       text.size=1,
                       border.color='black',
                       border.width=1
                       ){
    beta = min(wid/5, (wid+len)/20)
    gamma = wid/2

    if (clone.shape == 'polygon'){
        xx = c(x, x+beta, x+len, x+len, x+beta)
        yy = c(y, y+gamma, y+gamma, y-gamma, y-gamma)
        polygon(xx, yy, border=border.color, col=col, lwd=border.width)
    }else if(clone.shape == 'bell'){
        beta = min(wid/5, (wid+len)/10, len/3)
        xx0= c(x, x+beta, x+len, x+len, x+beta)
        yy0 = c(y, y+gamma, y+gamma, y-gamma, y-gamma)
        #polygon(xx, yy, border='black', col=col, lwd=0.2)

        gamma.shift = min(bell.curve.step, 0.5*gamma)
        
        # this is to prevent a coeff from being NaN when curve is generated below
        #if (beta <= 0.25){beta = 0.3}
        zeta = min(0.25, max(beta-0.1,0))

        x0=x+zeta; y0=0; x1=x+beta; y1 = gamma - gamma.shift
        n = 3; n = 1 + len/3
        a = ((y0^n-y1^n)/(x0-x1))^(1/n)
        b = y0^n/a^n - x0
        c = y
        #cat('a=', a, 'b=', b, 'c=', c, 'gamma=', gamma, 'len=', len, 'x0=',
        #    x0, 'x1=', x1, 'y0=', y0, 'y1=', y1,'\n')
        #curve(a*(x+b)^(1/n)+c, n=501, add=T, col=col, xlim=c(x0,x1))
        #curve(-a*(x+b)^(1/n)+c, n=501, add=T, col=col, xlim=c(x0,x1))

        beta0 = beta/5
        if (x0+beta0 > x1){beta0 = (x1-x0)/10}
        gamma0 = gamma/10
        
        xx = seq(x0+beta0,x1,(x1-x0)/100)
        yy = a*(xx+b)^(1/n)+c
        yy = c(y, yy, y+gamma, y-gamma, -a*(rev(xx)+b)^(1/n)+c)
        xx = c(x, xx, x+len, x+len, rev(xx))
        polygon(xx, yy, border=border.color, col=col, lwd=border.width)
        
        # generate some points to depict cells
        # buggy, does not work yet, and looks ugly
        if(F){
            cells = generate.fill.points(xx, yy)
            parNew = par('new')
            par(new=T)
            co = par('usr')
            plot(cells$x, cells$y, pch=20, axes=F, col='black', cex=0.1,
                xlim=co[1:2], ylim=co[3:4])
            par(new=parNew)
        }

        #xxx <<- xx; yyy <<- yy; ccc <<- cells
        #pdf('tmp.pdf');plot(xxx,yyy); co =par('usr'); par(new=T); plot(ccc$x, ccc$y, col='red', xlim=co[1:2], ylim=co[3:4], axes=F); dev.off(); dev.off(); dev.off()
        #stop()

    }else if (clone.shape == 'triangle'){
        #TODO: this does not work well yet. Implement!
        xx = c(x, x+len, x+len)
        yy = c(y/10, y+gamma, y-gamma)
        y = y/10
        polygon(xx, yy, border='black', col=col, lwd=0.2)
    }else if (clone.shape == 'parabol'){
        # TODO: Resovle overlapping (ie. subclone parabol expand outside of
        # parent clone parabol)
        x0=x; y0=0; x1=x+len; y1=gamma
        n = 3; n = 1 + len/3
        a = ((y0^n-y1^n)/(x0-x1))^(1/n)
        b = y0^n/a^n - x0
        c = y
        #cat('a=', a, 'b=', b, 'c=', c, 'gamma=', gamma, 'len=', len, 'x0=',
        #    x0, 'x1=', x1, 'y0=', y0, 'y1=', y1,'\n')
        curve(a*(x+b)^(1/n)+c, n=501, add=T, col=col, xlim=c(x0,x1))
        curve(-a*(x+b)^(1/n)+c, n=501, add=T, col=col, xlim=c(x0,x1))
        xx = seq(x0,x1,(x1-x0)/100)
        yy = a*(xx+b)^(1/n)+c
        yy = c(yy, -a*(rev(xx)+b)^(1/n)+c)
        xx = c(xx, rev(xx))
        #print(xx)
        #print(yy)
        polygon(xx, yy, col=col)
    }


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
            y.out <<- y.out + cell.frac.top.out.space
            x.out.shift <<- x.out.shift + 0.1
        }

        if (cell.frac.position == 'top.right' && clone.shape == 'bell'){
            angle = atan(gamma.shift/len*w2h.scale)*(180/pi)
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
            segments(x0, y0, x1, y1, col=cell.frac.side.arrow.col,
                     lwd=cell.frac.side.arrow.width)
            segments(x1, y1, x2, y2, col=cell.frac.side.arrow.col,
                     lwd=cell.frac.side.arrow.width)
            arrows(x0=x2, y0=y2, x1=x3, y1=y3,col=cell.frac.side.arrow.col,
                   length=0.025, lwd=cell.frac.side.arrow.width)
        }

        if (!is.null(variant.names)){
            if (is.null(variant.angle)){
                variant.angle = atan(gamma/beta)*(180/pi)# - 5
            }
            variant.x = x+beta/3+0.3*text.size
            variant.y = y+gamma/2-1.5*text.size
            variant.adj = c(0.5, 1)
            text(variant.x, variant.y, paste(variant.names, collapse='\n'),
                 cex=text.size*0.54, srt=variant.angle, adj=variant.adj,
                 col=variant.color)
        }
    }
}

#' Rescale VAF of subclones s.t. total VAF must not exceed parent clone VAF
#'
#' @description Rescale VAF of subclones s.t. total VAF must not exceed parent
#' clone VAF. When infered using bootstrap test, sometime the estimated
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
    #cat('Scaling called.\n')
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
    for (i in 1:nrow(v)){
        vi = v[i,]
        subs = v[!is.na(v$parent) & v$parent == vi$lab,]
        if (nrow(subs) == 0){next}
        vafs = subs$vaf
        margin = (vi$vaf - sum(vafs))/length(vafs)*scale
        sp = 0
        if (margin > 0){
            margin = margin*0.75
            sp = margin*0.25
        }
        spaces = rep(sp, length(vafs))
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

#' Determine which clones are subclone in a single sample, also determine what
#' clones are possible founder clones of the samples
#' @description: Determing which clones are subclone or founder clone or both or none
#' subclones are identified based on cellular
#' fractions of the ancestor clones. Eg. if a clone is a subclone, all of its
#' decendent clones are subclone. If a clone is not a subclone and has zero
#' cell frac, its only decendent clone is not a subclone
#' founder clones are identified as clones whose cell frac is non-zero and whose
#' parent has zero cell frac
#' @param v: data frame of subclonal structure as output of enumerate.clones
#' v must have row.names = v$lab
#' @param r: label of the clone where it and its decendent clones will be
#' evaluated. This is often the root of the tree
#'
# To do this, this function look at the root, and then flag all direct
# children of it as subclone/not subclone, then this repeats on all of
# its children
determine.subclone <- function(v, r){
    rownames(v) = v$lab
    next.clones = c(r)
    is.sub = rep(NA, nrow(v))
    names(is.sub) = v$lab
    is.founder = rep(NA, nrow(v))
    names(is.founder) = v$lab

    v$is.zero = ifelse(v$free.lower >= 0, F, T)

    # if no confidence interval estimated (no bootstrap model)
    if (all(is.na(v$free.lower))){
        v$is.zero = ifelse(v$free.mean > 0, T, F)
    }

    while (length(next.clones) > 0){
        cl = next.clones[1]
        children = v$lab[!is.na(v$parent) & v$parent == cl]
        next.clones = c(next.clones[-1], children)
        par = v[cl, 'parent'];
        if (!is.na(par) && (par == '-1' || par=='0')){
            # founding clone in monoclonal model, or clones
            # coming out of normal clone is not subclone
            is.sub[cl] = F
            is.founder[cl] = T
        }
        #if (v[cl, 'free.lower'] <= 0 && v[cl, 'num.subclones'] == 1){
        if (v[cl, 'is.zero'] && v[cl, 'num.subclones'] == 1){
            is.sub[children] = is.sub[cl]
            is.founder[children] = T
        }else{
            is.sub[children] = T
            if(v[cl, 'is.zero']){
                is.founder[children] = T
            }else{
                is.founder[children] = F
            }
        }
    }
    is.founder = is.founder & !v$is.zero
    return(list(is.sub=is.sub, is.founder=is.founder, is.zero=v$is.zero))
}

#' Get cellular fraction confidence interval
#'
#' @description Get cellular fraction confidence interval, and also determine
#' if it is greater than zero (eg. if it contains zero). This function return
#' a list with $cell.frac.ci = strings of cell.frac.ci, $is.zero.cell.frac =
#' tell if cell.frac.ci contains zero
#'
#' @param vi: clonal evolution tree data frame
#' @param include.p.value: include confidence level
#' @param sep: separator for the two confidence limits in output string
#'
get.cell.frac.ci <- function(vi, include.p.value=T, sep=' - '){
    cell.frac = NULL
    is.zero = NULL
    is.subclone = NULL 
    if('free.lower' %in% colnames(vi)){
        cell.frac.lower = ifelse(vi$free.lower == 0, '0',
                             gsub('\\.[0]+$|0+$', '',
                                  sprintf('%0.1f', 200*vi$free.lower)))
        cell.frac.upper = ifelse(vi$free.upper >= 0.5, '100%',
                             gsub('\\.[0]+$|0+$', '',
                                  sprintf('%0.1f%%', 200*vi$free.upper)))
        cell.frac = paste0(cell.frac.lower, sep , cell.frac.upper)
        if(include.p.value){
            cell.frac = paste0(cell.frac, '(',sprintf('%0.2f',
                            vi$free.confident.level),
                            #',p=', 1-vi$p.value,
                           ')')
            cell.frac = paste0(cell.frac, '/p=', sprintf('%0.3f', vi$p.value))
        }
        rownames(vi) = vi$lab
        
        # if only one clone as root, all cell.frac is NA
        # this is a dirty fix for output display
        # TODO: assign cell.frac for clone with zero subclone when
        # enumarating the models in enumerate.clones function.
        if(all(is.na(cell.frac.lower))){
            cell.frac[1] = '100-100%'
        }

    #if ('free.lower' %in% colnames(vi)){
        is.zero = ifelse(vi$free.lower >= 0, F, T)
        rownames(vi) = vi$lab
        names(is.zero) = vi$lab
        is.subclone = determine.subclone(vi,
            vi$lab[!is.na(vi$parent) & vi$parent == '-1'])$is.sub
    #}
    }

    #debug
    #print(vi$free.confident.level.non.negative)
    #if (vi$free.confident.level.non.negative == 0.687){
    #    print(vi$free.confident.level.non.negative)
    #    print(cell.frac)
    #    print(vi)
    #    vii <<- vi
    #}

    return(list(cell.frac.ci=cell.frac, is.zero.cell.frac=is.zero, is.subclone=is.subclone))
}

#' Draw clonal structures/evolution of a single sample
#'
#' @description Draw clonal structure for a sample with single or multiple
#' clones/subclones using polygon and tree plots
#'
#' @param v: clonal structure data frame (output of enumerate.clones)
#' @param clone.shape: c("bell", polygon"); shape of
#' the object used to present a clone in the clonal evolution plot.
#' @param adjust.clone.height: if TRUE, rescale the width of polygon such that
#' subclones should not have total vaf > that of parent clone when drawing
#' polygon plot
#' @param cell.frac.top.out.space: spacing between cell frac annotation when
#' annotating on top of the plot
#' @param cell.frac.side.arrow.width: width of the line and arrow pointing
#' to the top edge of the polygon from the cell frac annotation on top
#' @param color.border.by.sample.group: color border of bell plot based
#' on sample grouping
#'
#' @param variants.to.highlight: a data frame of 2 columns: cluster, variant.name
#' Variants in this data frame will be printed on the clone shape
#' @param bell.curve.step: see draw.clone function's bell.curve.step param
#' @param drop.zero.cell.frac.clone: c(T,F); if T, do not display zero cell
#' frac clones
#' @param clone.time.step.scale: scaling factor for distance between the tips
#' of the polygon/bell representing clone
#' @param zero.cell.frac.clone.color: color clone with zero cell fraction
#' in the sample with this color (default = NULL, color using matching color
#' auto-generated)
#' @param disable.cell.frac: disable cellular fraction display
draw.sample.clones <- function(v, x=2, y=0, wid=30, len=8,
                               clone.shape='bell',
                               bell.curve.step=0.25,
                               bell.border.width=1,
                               clone.time.step.scale=1,
                               label=NULL, text.size=1,
                               cell.frac.ci=F,
                               disable.cell.frac=F,
                               zero.cell.frac.clone.color=NULL,
                               zero.cell.frac.clone.border.color=NULL,
                               top.title=NULL,
                               adjust.clone.height=TRUE,
                               cell.frac.top.out.space=0.75,
                               cell.frac.side.arrow.width=1.5,
                               variants.to.highlight=NULL,
                               variant.color='blue',
                               variant.angle=NULL,
                               show.time.axis=TRUE,
                               color.node.by.sample.group=FALSE,
                               color.border.by.sample.group=TRUE){
    v = v[!v$excluded,]
    if (adjust.clone.height){
        #cat('Will call rescale.vaf on', label, '\n')
        #print(v)
        v = rescale.vaf(v)

    }
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
                # for bell curve, needs to shift x further to make sure
                # bell of subclone falls completely in its parent bell
                x.shift = 1 * ifelse(clone.shape=='bell', 1.2, 1)

                if (vi$y.shift + vi$vaf >= high.vaf && vi$vaf < low.vaf){
                    x.shift = 2*x.shift
                }
                if (clone.shape=='triangle'){
                    x.shift = x.shift + 1
                }
                par = v[v$lab == vi$parent,]

                if (vi$vaf < 0.05 && par$num.subclones > 1){x.shift = x.shift*2}
                x.shift = x.shift*clone.time.step.scale
                xi = par$x + x.shift
                
                yi = par$y - wid*par$vaf/2 + wid*vi$vaf/2 + vi$y.shift*wid
                leni = par$len - x.shift
            }
            #cell.frac.position = ifelse(vi$free.lower < 0.05 & vi$vaf > 0.25, 'side', 'top.right')
            #cell.frac.position = ifelse(vi$free.lower < 0.05, 'top.out', 'top.right')
            cell.frac.position = ifelse(vi$free < 0.05, 'top.out', 'top.right')
            #cell.frac.position = ifelse(vi$free < 0.05, 'top.out', 'right.mid')
            #cell.frac.position = ifelse(vi$free < 0.05, 'top.out', 'top.out')
            #cell.frac.position = ifelse(vi$num.subclones > 0 , 'right.top', 'right.mid')
            #cell.frac.position = 'top.mid'
            cell.frac = paste0(gsub('\\.[0]+$|0+$', '',
                                    sprintf('%0.2f', vi$free.mean*2*100)), '%')
            if(cell.frac.ci && !disable.cell.frac){
                cell.frac = get.cell.frac.ci(vi, include.p.value=T)$cell.frac.ci
            }else if (disable.cell.frac){
                cell.frac = NA
            }
            variant.names = variants.to.highlight$variant.name[
                variants.to.highlight$cluster == vi$lab]
            if (length(variant.names) == 0) {
                variant.names = NULL
            }
            clone.color = vi$color
            border.color='black'
            if (color.border.by.sample.group){
                border.color = vi$sample.group.color
            }else if (color.node.by.sample.group){
                clone.color = vi$sample.group.color
            }
            if (!is.null(zero.cell.frac.clone.color) & vi$is.zero){
                clone.color = zero.cell.frac.clone.color
            }
            if (!is.null(zero.cell.frac.clone.border.color) & vi$is.zero){
                border.color = zero.cell.frac.clone.border.color
                if (border.color == 'fill'){border.color = clone.color}
            }
            draw.clone(xi, yi, wid=wid*vi$vaf, len=leni, col=clone.color,
                       clone.shape=clone.shape,
                       bell.curve.step=bell.curve.step,
                       border.width=bell.border.width,
                       label=vi$lab,
                       cell.frac=cell.frac,
                       cell.frac.position=cell.frac.position,
                       cell.frac.side.arrow.col=clone.color,
                       text.size=text.size,
                       cell.frac.top.out.space=cell.frac.top.out.space,
                       cell.frac.side.arrow.width=cell.frac.side.arrow.width,
                       variant.names=variant.names,
                       variant.color=variant.color,
                       variant.angle=variant.angle,
                       border.color=border.color)
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
        # draw time axis
        if (show.time.axis && i==1){
            axis.y = -9
            arrows(x0=x,y0=axis.y,x1=10,y1=axis.y, length=0.05, lwd=0.5)
            text(x=10, y=axis.y-0.75, label='time', cex=1, adj=1)
            segments(x0=x,y0=axis.y-0.2,x1=x, y1=axis.y+0.2)
            text(x=x,y=axis.y-0.75,label='Cancer initiated', cex=1, adj=0)
            segments(x0=x+len,y0=axis.y-0.2,x1=x+len, y1=axis.y+0.2)
            text(x=x+len, y=axis.y-0.75, label='Sample taken', cex=1, adj=1)
        }
    }
    plot(c(0, 10),c(-10,10), type = "n", xlab='', ylab='', xaxt='n',
         yaxt='n', axes=F)
    if (!is.null(label)){
        text(x-1, y, label=label, srt=90, cex=text.size, adj=c(0.5,1))
    }
    if (!is.null(top.title)){
        text(x, y+10, label=top.title, cex=(text.size), adj=c(0,0.5))
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
make.graph <- function(v, cell.frac.ci=TRUE, node.annotation='clone', node.colors=NULL){
    library(igraph)
    #v = v[!is.na(v$parent),]
    #v = v[!is.na(v$parent) | v$vaf != 0,]
    v = v[!is.na(v$parent),]
    #rownames(v) = seq(1,nrow(v))
    rownames(v) = v$lab
    g = matrix(0, nrow=nrow(v), ncol=nrow(v))
    rownames(g) = rownames(v)
    colnames(g) = rownames(v)
    if (nrow(v) == 0){# single sample, no pruned tree
        return(NULL)
    }
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
        cell.frac = get.cell.frac.ci(v, include.p.value=F, sep=' -\n')$cell.frac.ci
    }
    labels = v$lab
    colors = v$color
    if (!is.null(node.colors)){
        colors = node.colors[labels]
    }
    if (!is.null(cell.frac) && !all(is.na(cell.frac))){
        labels = paste0(labels,'\n', cell.frac)
    }
    
    # add sample name
    # trick to strip off clone having zero cell.frac and not a founding clone of a sample
    # those samples prefixed by 'o*'
    remove.founding.zero.cell.frac = F
    if (node.annotation == 'sample.with.cell.frac.ci.founding.and.subclone'){
        node.annotation = 'sample.with.cell.frac.ci'
        remove.founding.zero.cell.frac = T
    }
    if (node.annotation != 'clone' && node.annotation %in% colnames(v)){
        # this code is to add sample to its terminal clones only, obsolete
        #leaves = !grepl(',', v$sample)
        #leaves = v$is.term
        #if (any(leaves)){
        #    #labels[leaves] = paste0('\n', labels[leaves], '\n', v$leaf.of.sample[leaves])
        #}
        has.sample = !is.na(v[[node.annotation]])
        samples.annot = v[[node.annotation]][has.sample]
        if (!cell.frac.ci){
            # strip off cell.frac.ci
            #tmp = unlist(strsplit(',', samples.annot))
            #tmp = gsub('\\s*:\\s*[^:].+', ',', tmp)
            samples.annot = gsub('\\s*:\\s*[^:]+(,|$)', ',', v[[node.annotation]][has.sample])
            samples.annot = gsub(',$', '', samples.annot)
        }
        if (remove.founding.zero.cell.frac){
            samples.annot = gsub('o\\*[^,]+(,|$)', '', samples.annot)
        }
        labels[has.sample] = paste0(labels[has.sample], '\n', samples.annot)
    }
    V(g)$name = labels
    V(g)$color = colors
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
#' @param show.sample: show sample names in node
#' @param node.annotation: c('clone', 'sample.with.cell.frac.ci',
#' 'sample.with.nonzero.cell.frac.ci', 'sample.with.cell.frac.ci',
#' 'sample.with.cell.frac.ci.founding.and.subclone')
#' Labeling mode for tree node. 'sample.with.cell.frac.ci' = all samples where clone's signature
#' mutations detected togeter with cell.frac.ci if cell.frac.ci=T; 'sample.nonzero.cell.frac'
#' = samples where clones detected with nonzero cell fraction determined by subclonal.test
#' 'sample.with.cell.frac.ci.founding.and.subclone': show all execpt samples where cell frac is
#' zero and not subclone
#' if cell.frac is diff from zero; default = 'clone' = only clone label is shown
#' 'sample.term.clone' = samples where clones are terminal clones (ie. clones that do not
#' have subclones in that sample)
#' @param node.label.split.character: replace this character by "\n" to allow multi
#' line labeling for node; esp. helpful with multi samples sharing a
#' node being plotted.
#' @param node.colors: named vector of colors to plot for nodes, names = clone/cluster
#' labels; if NULL (default), then use v$color
#' @param color.node.by.sample.group: if TRUE, and if column sample.group and
#' sample.group.color exists, then color node by these; also add legend
#' @param color.border.by.sample.group: if TRUE, color border of clones in tree
#' or bell plot based on sample grouping
#' @param show.legend: show sample group legends, etc.
#'
plot.tree <- function(v, node.shape='circle', display='tree',
                      node.size=50,
                      node.colors=NULL,
                      color.node.by.sample.group=FALSE,
                      color.border.by.sample.group=TRUE,
                      show.legend=T,
                      tree.node.text.size=1,
                      cell.frac.ci=T,
                      node.prefix.to.add=NULL,
                      title='',
                      #show.sample=FALSE,
                      node.annotation='clone',
                      node.label.split.character=NULL,
                      node.num.samples.per.line=NULL,
                      out.prefix=NULL,
                      graphml.out=FALSE,
                      out.format='graphml'){
    library(igraph)
    grps = NULL
    grp.colors = 'black'
    if (color.border.by.sample.group){
        color.node.by.sample.group = F #disable coloring node by group if blanket is used
        #grps = list()
        #for (i in 1:nrow(v)){
        #    grps = c(grps, list(i))
        #}
        grp.colors = v$sample.group.color
        # get stronger color for borders
        #uniq.colors = unique(grp.colors)
        #border.colors = get.clonevol.colors(length(uniq.colors), T)
        #names(border.colors) = uniq.colors
        #grp.colors = border.colors[grp.colors]
        #v$sample.group.border.color = grp.colors
    }else if (color.node.by.sample.group){
        node.colors = v$sample.group.color
        names(node.colors) = v$lab
    }
    #x = make.graph(v, cell.frac.ci=cell.frac.ci, include.sample.in.label=show.sample, node.colors)
    x = make.graph(v, cell.frac.ci=cell.frac.ci, node.annotation=node.annotation, node.colors)
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

    vertex.labels = V(g)$name
    #vlabs <<- vertex.labels
    if (!is.null(node.label.split.character)){
        num.splits = sapply(vertex.labels, function(l)
            nchar(gsub(paste0('[^', node.label.split.character, ']'), '', l)))

        # only keep the node.label.split.char in interval of node.num.samples.per.line
        # such that a block of node.num.samples.per.line samples will be grouped and
        # kept in one line
        if (!is.null(node.num.samples.per.line)){
            for (i in 1:length(vertex.labels)){
                vl = unlist(strsplit(vertex.labels[i], node.label.split.character))
                sel = seq(min(node.num.samples.per.line, length(vl)),
                    length(vl),node.num.samples.per.line)
                vl[sel] = paste0(vl[sel], node.label.split.character)
                vl[-sel] = paste0(vl[-sel], ';')
                vertex.labels[i] = paste(vl, collapse='')
            }
            num.splits = length(sel) + 1
        }
        extra.lf = sapply(num.splits, function(n) paste(rep('\n', n), collapse=''))
        vertex.labels = paste0(extra.lf, gsub(node.label.split.character, '\n',
            vertex.labels))
    }

    plot(g, edge.color='black', layout=layout, main=title,
         edge.arrow.size=0.75, edge.arrow.width=0.75,
         vertex.shape=node.shape, vertex.size=node.size,
         vertex.label.cex=tree.node.text.size,
         #vertex.label.color=sample(c('black', 'blue', 'darkred'), length(vertex.labels), replace=T),
         vertex.label=vertex.labels,
         #mark.groups = grps,
         #mark.col = 'white',
         #mark.border = grp.colors,
         vertex.frame.color=grp.colors)
         #, vertex.color=v$color, #vertex.label=labels)
    if ((color.node.by.sample.group || color.border.by.sample.group) & show.legend &
            'sample.group' %in% colnames(v)){
        vi = unique(v[!v$excluded & !is.na(v$parent),
            c('sample.group', 'sample.group.color')])
        vi = vi[order(vi$sample.group),]
        if (color.border.by.sample.group){
            legend('topright', legend=vi$sample.group, pt.cex=3, cex=1.5,
                 pch=1, col=vi$sample.group.color)
        }else{
            legend('topright', legend=vi$sample.group, pt.cex=3, cex=1.5,
                 pch=16, col=vi$sample.group.color)
        }
        legend('bottomleft', legend=c('*  sample founding clone',
                                    '°  zero cellular fraction',
                                   '°* ancestor of sample founding clone'
                                  ),
                                  pch=c('', '', ''))
        # events on each clone legend
        if ('events' %in% colnames(v)){
            ve = v[v$events != '',]
            ve = ve[order(as.integer(ve$lab)),]
            # only print 5 events per-line
            ve$events = insert.lf(ve$events, 5, ',')
            legend('topleft', legend=paste0(sprintf('%2s', ve$lab), ': ', ve$events),
                    pt.cex=2, cex=1, pch=19, col=ve$color)
        }
    }

    # remove newline char because Cytoscape does not support multi-line label
    V(g)$name = gsub('\n', ' ', V(g)$name, fixed=T)
    if (!is.null(node.prefix.to.add)){
        V(g)$name = paste0(node.prefix.to.add, V(g)$name)
    }

    if (!is.null(out.prefix)){
        out.file = paste0(out.prefix, '.', out.format)
        #cat('Writing tree to ', out.file, '\n')
        if (graphml.out){
            write.graph(g, file=out.file, format=out.format)
        }
    }

    return(g)

}

#' Write clonal evolution tree to file
#' TODO: do!
write.tree <- function(v, out.file, out.format='tabular'){
    v = v[, c('parent', 'lab', '')]
}

get.model.score <- function(v){
    #return(prod(v$p.value[!is.na(v$p.value)]))
    return(max(v$p.value[!is.na(v$p.value)]))
}

#' Get all set of subclones for all clone across samples
#' together with p value
get.subclones.across.samples <- function(x, matched.model.index){
    samples = names(x$models)
    tree = x$matched$merged.trees[[matched.model.index]]
    labs = tree$lab[!tree$excluded]
    # look for subclones in each sample for each clone
    subs = NULL
    for (s in samples){
        m = x$models[[s]][[x$matched$index[matched.model.index, s]]]
        for (cl in labs){
            sc = m$lab[!is.na(m$parent) & m$parent == cl & !m$excluded]
            p = m$p.value[m$lab == cl]
            if (length(sc) > 0){
                sc = paste(sort(sc), collapse=',')
                r = data.frame(lab=cl, sample=s, subclones=sc, p=p,
                    stringsAsFactors=F)
                if(is.null(subs)){subs = r}else{subs = rbind(subs,r)}
            }
        }
    }

    subs = subs[order(subs$lab),]
    return(subs)
}

#' Apply cross rule to all clones in all matched models using p-value combination
#' @description: For each model, each clone will receive a score that
#' is equal to the  combined p-value of the test that the CCF of the
#' clone is >= 0
#' @param x: output of infer.clonal.models
#' @param meta.p.method: method for combining p-values across samples
#' values = c('fisher', 'z'), default = 'fisher'
#' @param exhaustive.mode: placeholder for exhaustive.mode, not implemented yet.
#
cross.rule.score <- function(x, meta.p.method='fisher', exhaustive.mode=F, boot=NULL){
    if (!is.null(x$matched) && x$num.matched.models > 0 && ncol(x$matched$index) > 1){
        samples = names(x$models)
        num.models = nrow(x$matched$index)
        x$matched$scores$max.clone.ccf.combined.p = NA
        x$matched$clone.ccf.pvalues = list()
        # foreach matched model, recalc score by combining p values across
        # samples for each clone
        for (i in 1:num.models){
            trees = NULL
            t = x$match$merged.trees[[i]]
            p = NULL
            for (s in samples){
                mi = x$models[[s]][[x$match$index[i,s]]]
                mi = mi[!mi$excluded & !is.na(mi$parent), c('lab', 'p.value')]
                colnames(mi) = c('lab', s)
                if (is.null(p)){
                    p = mi
                }else{
                    p = merge(p, mi, all=T)
                }
            }
            ppp <<- p
            if (ncol(p) == 2){#single sample
                p$cmb.p = apply(p[,c(2,2)], 1, combine.p, method=meta.p.method)
            }else{
                p$cmb.p = apply(p[,-1], 1, combine.p, method=meta.p.method)
            }
            # model score = max (combined p of each clone)
            x$matched$scores$max.clone.ccf.combined.p[i] = max(p$cmb.p)
            x$matched$merged.trees[[i]]$clone.ccf.combined.p = p$cmb.p
            # save the whole pvalue matrix
            x$matched$clone.ccf.pvalues[[i]] = p
        }
        # order matched models by new score
        idx = order(x$matched$scores$max.clone.ccf.combined.p)
        x$matched$index = x$matched$index[idx,, drop=F]
        x$matched$scores = x$matched$scores[idx,, drop=F]
        # order merged trees
        tmp = list()
        for (i in idx){
            tmp = c(tmp, list(x$matched$merged.trees[[i]]))
        }
        x$matched$merged.trees = tmp
        # order merged traces
        tmp = list()
        for (i in idx){
            tmp = c(tmp, list(x$matched$merged.traces[[i]]))
        }
        x$matched$merged.traces = tmp
        # order pvalues
        tmp = list()
        for (i in idx){
            tmp = c(tmp, list(x$matched$clone.ccf.pvalues[[i]]))
        }
        x$matched$clone.ccf.pvalues = tmp

        # remove previous model scores (which was very small probability)
        x$matched$scores$model.score = NULL
    }
    return(x)
}

#' Merge clonnal evolution trees from multiple samples into a single tree
#' 
#' @description Merge a list of clonal evolution trees (given by the clonal
#' evolution tree data frames) in multiple samples into a single clonal
#' evolution tree, and label the leaf nodes (identified in individual tree)
#' with the corresponding samples in the merged tree.
#' 
#' @param trees: a list of clonal evolution trees' data frames
#' @param samples: name of samples that will be used in node labels
#' @param sample.groups: named vector of sample grouping
#' @param merge.similar.samples: drop a sample if there is already
#' another sample with the same tree
#' 
merge.clone.trees <- function(trees, samples=NULL, sample.groups=NULL, merge.similar.samples=F){
    # to keep track of what samples merged with what samples
    merged.trace = NULL

    if (merge.similar.samples){
        # remove similar trees
        z = trim.clone.trees(trees, samples, remove.sample.specific.clones=F)
        merged.trace = z$merged.trace
        trees = z$unique.trees
        samples = names(trees)
        sample.groups = sample.groups[samples]
    }
    n = length(trees)
    merged = NULL
    if (is.null(samples)){samples = seq(1,n)}
    #leaves = c()
    lf = NULL
    ccf.ci = NULL
    ccf.ci.nonzero = NULL
    subclones = NULL #subclonal status
    cgrp = NULL #grouping clones based on sample groups
    #let's group all samples in one group if sample groups not provided
    if (is.null(sample.groups)){
        sample.groups = rep('group1', length(samples))
        names(sample.groups) = samples
    }

    #TODO: there is a bug in infer.clonal.models that did not give
    # consistent ancestors value across samples, let's discard this
    # column now for merging, but later need to fix this.
    #v = trees[[i]][, c('lab', 'color', 'parent', 'ancestors', 'excluded')]
    key.cols = c('lab', 'color', 'parent', 'excluded')
 
    for (i in 1:n){
        v = trees[[i]] 
        s = samples[i]

        # get cell.frac
        v = v[!v$excluded & !is.na(v$parent),]
        if (nrow(v) == 0){stop('ERROR: Something wrong. No clones left after filter. They might have been excluded.\n')}
        # TODO: scale.cell.frac here works independent of plot.clonal.models
        # which has a param to ask for scaling too. Make them work together nicely.
        # also, clonal tree data frame v is now have two more column indicating
        # if a clone is.subclone or is.zero cell frac, so getting these info via
        # get.cell.frac.ci is redundant and potentially create inconsistency if code changes
        # TODO: utilize is.subclone and is.zero columns in v
        cia = get.cell.frac.ci(scale.cell.frac(v), sep='-')
        #ci = data.frame(lab=v$lab, sample.with.cell.frac.ci=paste0(ifelse(cia$is.subclone,
        #    '', '*'), s, ' : ', cia$cell.frac.ci), stringsAsFactors=F)
        ci = data.frame(lab=v$lab, sample.with.cell.frac.ci=paste0(ifelse(v$is.founder,
            '*', ''), s, ' : ', cia$cell.frac.ci), stringsAsFactors=F)
        #ci.nonzero = ci[!is.na(cia$is.zero.cell.frac) & !cia$is.zero.cell.frac,]
        ci.nonzero = ci[!cia$is.zero.cell.frac,]
        ci$sample.with.cell.frac.ci[cia$is.zero.cell.frac] = paste0('°',
            ci$sample.with.cell.frac.ci[cia$is.zero.cell.frac])
        if (is.null(ccf.ci)){ccf.ci = ci}else{ccf.ci = rbind(ccf.ci, ci)}
        if (is.null(ccf.ci.nonzero)){ccf.ci.nonzero = ci.nonzero}else{
            ccf.ci.nonzero = rbind(ccf.ci.nonzero, ci.nonzero)}
        v$sample = s
        #v$sample[!cia$is.subclone] = paste0('*', v$sample[!cia$is.subclone])
        v$sample[v$is.founder] = paste0('*', v$sample[v$is.founder])
        v$sample[cia$is.zero.cell.frac] = paste0('°', v$sample[cia$is.zero.cell.frac])
        this.leaves = v$lab[!is.na(v$parent) & !(v$lab %in% v$parent)]
        this.lf = data.frame(lab=this.leaves, leaf.of.sample=s, stringsAsFactors=F)
        if (is.null(lf)){lf = this.lf}else{lf = rbind(lf, this.lf)}
        #leaves = c(leaves, this.leaves)
 
        #clone grouping only non.zero cell.frac clones
        vz = v[!cia$is.zero.cell.frac,]
        if (!is.null(sample.groups)){#this is uneccesary if given default grouping above
            cg = data.frame(lab=vz$lab, sample.group=sample.groups[s],
                stringsAsFactors=F, row.names=NULL)
            if (is.null(cgrp)){cgrp = cg}else{cgrp = rbind(cgrp, cg)}
        }
       
        # keep only key.cols and sample cols for merging
        v = v[, c(key.cols, 'sample')]
        if (is.null(merged)){merged = v}else{merged = rbind(merged, v)}
    }
    merged = merged[!is.na(merged$parent),]

    #merged = unique(merged)
    merged = aggregate(sample ~ ., merged, paste, collapse=',')
    #leaves = unique(leaves)
    lf = aggregate(leaf.of.sample ~ ., lf, paste, collapse=',')
    lf$is.term = T
    merged = merge(merged, lf, all.x=T)

    ccf.ci = aggregate(sample.with.cell.frac.ci ~ ., ccf.ci, paste, collapse=',')
    merged = merge(merged, ccf.ci, all.x=T)

    ccf.ci.nonzero = aggregate(sample.with.cell.frac.ci ~ ., ccf.ci.nonzero,
                                paste, collapse=',')
    colnames(ccf.ci.nonzero) = c('lab', 'sample.with.nonzero.cell.frac.ci')
    merged = merge(merged, ccf.ci.nonzero, all.x=T)


    if (!is.null(cgrp)){
        #print(cgrp)
        cgrp = unique(cgrp)
        cgrp = cgrp[order(cgrp$sample.group),]
        cgrp = aggregate(sample.group ~ ., cgrp, paste, collapse=',')
        #print(cgrp)
        sample.grps = unique(cgrp$sample.group)
        sample.group.colors = get.clonevol.colors(length(sample.grps), strong.color=T)
        names(sample.group.colors) = sample.grps
        cgrp$sample.group.color = sample.group.colors[cgrp$sample.group]
        merged = merge(merged, cgrp, all.x=T)
    }

    
    #leaves = unique(lf$lab)
    #merged$is.term = F
    #merged$is.term[merged$lab %in% leaves] = T
    merged$is.term[is.na(merged$is.term)] = F
    merged$num.samples = sapply(merged$sample, function (l)
        length(unlist(strsplit(l, ','))))
    merged$leaf.of.sample.count = sapply(merged$leaf.of.sample, function (l)
        length(unlist(strsplit(l, ','))))
    merged$num.samples[is.na(merged$num.samples)] = 0
    merged$leaf.of.sample.count[is.na(merged$leaf.of.sample)] = 0
    rownames(merged) = merged$lab
    return (list(merged.tree=merged, merged.trace=merged.trace))
}


#' Compare two clonal evolution trees
#' 
#' @description Compare two clonal evolution trees to see if they differ
#' assuming excluded nodes are already excluded, labels are ordered
#' 
#' Return F if they are different, T if they are the same
#' 
#'
#' @params v1: clonal evolution tree data frame 1
#' @params v2: clonal evolution tree data frame 2
#'
#'

compare.clone.trees <- function(v1, v2){
    res = F
    if (nrow(v1) == nrow(v2)){
        if (all(v1$parent == v2$parent) & all(v1$lab == v2$lab)){
            res = T
        }
    }
    return(res)
}

#' Reduce clone trees to core trees after excluding clones
#' that are present in only a single sample
#'
#' @description Reduce clone trees produced by infer.clonal.models
#' to the core models, ie. ones that do not affect how clonal evolve
#' and branch across samples. This function will strip off the clones
#' that involve only one sample, any excluded nodes, and compare the trees
#' and keep only ones that are different. Return a list of merged.trees
#' 
#' @param merged.trees: output of infer.clonal.models()$matched$merged.trees
#' or any list of trees produced by infer.clonal.models()$models
#' @param remove.sample.specific.clones: remove clones that are found in
#' only one sample before comparing, default=T (used when reducing
#' merged trees. Set this to F to compare and reduce multiple samples
#' with the same tree to one sample
# TODO: strip off info in cell.frac (currently keeping it for convenient
# plotting
trim.clone.trees <- function(merged.trees, remove.sample.specific.clones=T, samples=NULL){

    n = length(merged.trees)
    if (n == 0){
        return(list(unique.trees=NULL, merged.trace=NULL))
    }
    # trim off sample specific clones (if requested) and excluded nodes, sort by label
    for (i in 1:n){
        v = merged.trees[[i]]
        v = v[!v$excluded,]
        if (remove.sample.specific.clones){
            v = v[v$num.samples > 1,]
        }
        v = v[order(v$lab),]
        merged.trees[[i]] = v
    }

    #ttt2 <<- merged.trees

    # compare and reduce
    idx = seq(1,n)
    i = 1;
    merged.trace = c(1,1)
    while (i < n){
        j = i + 1
        if (i > 1){merged.trace = rbind(merged.trace, c(idx[i],idx[i]))}
        while (j <= n){
            #cat(samples[idx[i]], samples[idx[j]], '\t')
            if(compare.clone.trees(merged.trees[[i]], merged.trees[[j]])){
                #cat('Drop tree', j, '\n')
                merged.trees[[j]] = NULL
                n = n - 1
                merged.trace = rbind(merged.trace, c(idx[i], idx[j]))
                idx = idx[-j]
                #cat('equal\n')
            }else{
                j = j + 1
                #cat('diff\n')
            }
        }
        i = i + 1
    }
    #print(idx)
    #print(samples)
    
    # today debug
    # mtr <<- merged.trace
    if (is.null(dim(merged.trace))){# one row, 1 tree to merge,
        #need to make matrix for followed code to work
        merged.trace = matrix(merged.trace, nrow=1)
    }

    colnames(merged.trace) = c('sample', 'similar.sample')
    # last sample that were not merged with any other
    merged.trace = as.data.frame.matrix(merged.trace)
    rownames(merged.trace) = NULL
    if (!(idx[n] %in% c(merged.trace$sample, merged.trace$similar.sample))){
        merged.trace = rbind(merged.trace, c(idx[n], idx[n]))
        #cat(idx[n], '\n')
        #stop('HMMM')
    }

    if (!is.null(samples)){
        names(merged.trees) = samples[idx]
        merged.trace$sample = samples[merged.trace$sample]
        merged.trace$similar.sample = samples[merged.trace$similar.sample]
    }


    return(list(unique.trees=merged.trees, merged.trace=merged.trace))
}

#' Compare two merged clonal evolution trees
#' 
#' @description Compare two clonal evolution trees to see if they differ
#' and if they also differ if all termial node (leaves) are removed.
#' This is useful to determine the number of unique models ignoring
#' sample specific clones which often give too many models when they
#' are present and low frequency. This function return 0 if tree match
#' at leaves, 1 if not match at leave but match when leaves are removed
#' 2 if not matching at internal node levels
#' 
#'
#' @params v1: clonal evolution tree data frame 1
#' @params v2: clonal evolution tree data frame 1
#'
#'

compare.clone.trees.removing.leaves <- function(v1, v2){
    res = 2
    v1 = v1[!v1$excluded,]
    v2 = v2[!v1$excluded,]
    v1 = v1[order(v1$lab),]
    v2 = v2[order(v2$lab),]

    if (nrow(v1) == nrow(v2)){
        if (all(v1$parent == v2$parent)){
            res = 0
        }
    }
    if (res !=0){
        # remove leaves
        v1 = v1[!is.na(v1$parent) & (v1$lab %in% v1$parent),]
        v2 = v2[!is.na(v2$parent) & (v2$lab %in% v2$parent),]
        if (all(v1$parent == v2$parent)){
            res = 1
        }
    }

    
    return(res)
}



#' Find matched models between samples
#' infer clonal evolution models, given all evolve from the 1st sample
#' @description Find clonal evolution models across samples that are
#' compatible, given the models for individual samples
#' @param merge.similar.samples: see merge.clone.trees
# TODO: recursive algorithm is slow, improve.
find.matched.models <- function(vv, samples, sample.groups=NULL, merge.similar.samples=F){
    cat('Finding matched clonal architecture models across samples...\n')
    nSamples = length(samples)
    matched = NULL
    scores = NULL
    # for historical reason, variables were named prim, met, etc., but
    # it does not mean samples are prim, met.
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
    num.models.found = ifelse(is.null(matched), 0, nrow(matched))
    cat('Found ', num.models.found, 'compatible model(s)\n')

	# merge clonal trees
    merged.trees = list()
    merged.traces = list()
    if (num.models.found > 0){
        cat('Merging clonal evolution trees across samples...\n')
        for (i in 1:num.models.found){
            m = list()
            for (j in 1:nSamples){
                m = c(m, list(vv[[j]][[matched[i, j]]]))
            }
            
            zz = merge.clone.trees(m, samples=samples, sample.groups, merge.similar.samples=merge.similar.samples)
            mt = zz$merged.tree
            trace = zz$merged.trace
            # after merged, assign sample.group and color to individual tree
            #print(mt)
            for (j in 1:nSamples){
                vv[[j]][[matched[i, j]]] = merge(vv[[j]][[matched[i, j]]],
                    mt[, c('lab', 'sample.group', 'sample.group.color')], all.x=T)
            }
 
            merged.trees = c(merged.trees, list(mt))
            merged.traces = c(merged.traces, list(trace))
        }
    }
	
    return(list(models=vv, matched.models=matched, merged.trees=merged.trees,
        merged.traces=merged.traces, scores=scores))
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
#' @param cluster.col.name: column that holds the cluster identity, overwriting
#' the default 'cluster' column name
#' @param founding.cluster: the cluster of variants that are found in all
#' samples and is beleived the be the founding events. This is often
#' the cluster with highest VAF and most number of variants
#' @param ignore.clusters: clusters to ignore (not inluded in the models). Those
#' are the clusters that are thought of as outliers, artifacts, etc. resulted
#' from the error or bias of the sequencing and analysis. This is provided as
#' a debugging tool
#' @param sample.groups: named vector of sample groups, later clone will be
#' colored based on the grouping of shared samples, eg. clone specific to
#' primaries, metastasis, or shared between them. Default = NULL
#' @param model: cancer evolution model to use, c('monoclonal', 'polyclonal').
#' monoclonal model assumes the orginal tumor (eg. primary tumor) arises from
#' a single normal cell; polyclonal model assumes the original tumor can arise
#' from multiple cells (ie. multiple founding clones). In the polyclonal model,
#' the total VAF of the separate founding clones must not exceed 0.5
#' @param subclonal.test: 'bootstrap' = perform bootstrap subclonal test
#' 'none' = straight comparison of already estimated VAF for each cluster
#' provided in c
#' @param subclonal.test.model: What model to use when generating the bootstrap
#' Values are: c('non-parametric', 'normal', 'normal-truncated', 'beta',
#' 'beta-binomial')
#' @param min.cluster.vaf: the minimum cluster VAF to be considered positive
#' (detectable cluster); default=NULL, this will be used in two places:
#' (i): detection of positive VAF cluster, if mean/median cluster VAF is
#' greater; if not provided (NULL), bootstrap test will be used instead
#' to determine if cluster VAF is significantly greater/smaller than zero (using
#' the sum.p.cutoff param)
#' (ii): when no bootstrap model used, any cluster VAF falling below this
#' is considered non-existed/non-detectable cluster
#' @param cluster.center: median or mean
#' @param random.seed: a random seed to bootstrap generation.
#' @param merge.similar.samples: if a latter sample has the same tree
#' compared with a sample processed earlier (appear first in vaf.col.names)
#' then that sample will be removed from the tree when merging clonal
#' evolution trees across samples. An output file *.sample-reduction.tsv
#' will be created when plot.clonal.models is called later.
#' @param clone.colors: vector of colors that will be used for the clone
#' drawing in the results
infer.clonal.models <- function(c=NULL, variants=NULL,
                                cluster.col.name='cluster',
                                founding.cluster=NULL,
                                ignore.clusters=NULL,
                                vaf.col.names=NULL,
                                vaf.in.percent=TRUE,
                                sample.names=NULL,
                                sample.groups=NULL,
                                model='monoclonal',
                                subclonal.test='none',
                                cluster.center='median',
                                subclonal.test.model='non-parametric',
                                merge.similar.samples=F,
                                clone.colors=NULL,
                                random.seed=NULL,
                                boot=NULL,
                                num.boots=1000,
                                p.value.cutoff=NULL,
                                sum.p.cutoff=0.01,
                                cross.p.cutoff=NULL,
                                alpha=NULL,
                                min.cluster.vaf=NULL,
                                verbose=TRUE){
    # backward compatible with old p.value.cutoff
    if (!is.null(p.value.cutoff)){sum.p.cutoff = p.value.cutoff}
    if (is.null(alpha)){alpha = sum.p.cutoff}
    if (is.null(cross.p.cutoff)){cross.p.cutoff = sum.p.cutoff}

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
        founding.cluster = '0'
        # add a faked normal clone with VAF = norm(mean=50, std=10)
        # TODO: follow the distribution chosen by user
        if (add.normal){
            tmp = variants[rep(1,100),]
            tmp[[cluster.col.name]] = founding.cluster
            vaf50 = matrix(rnorm(100, 50, 10), ncol=1)[,rep(1, length(vaf.col.names))]
            tmp[, vaf.col.names] = vaf50
            variants = rbind(tmp, variants)
        }

    }
    if (is.na(add.normal)){
        stop(paste0('ERROR: Model ', model, ' not supported!\n'))
    }
    if(verbose){cat('Using ', model, ' model\n', sep='')}

    # prepare cluster data and infer clonal models for individual samples
    if (is.null(c)){
        c = estimate.clone.vaf(variants, cluster.col.name,
                               vaf.col.names, vaf.in.percent=vaf.in.percent,
                               method=cluster.center)
        #print(c)
    }
    vv = list()
    for (i in 1:nSamples){
        s = vaf.col.names[i]
        sample.name = sample.names[i]
        v = make.clonal.data.frame(c[[s]], c[[cluster.col.name]],
            colors=clone.colors)
        if (subclonal.test == 'none'){
            #models = enumerate.clones.absolute(v)
            models = enumerate.clones(v, sample=s,
                                      founding.cluster=founding.cluster,
                                      min.cluster.vaf=min.cluster.vaf,
                                      ignore.clusters=ignore.clusters)
        }else if (subclonal.test == 'bootstrap'){
            if (is.null(boot)){
                #boot = generate.boot(variants, vaf.col.names=vaf.col.names,
                #                     vaf.in.percent=vaf.in.percent,
                #                     num.boots=num.boots)

                boot = generate.boot(variants, vaf.col.names=vaf.col.names,
                                     vaf.in.percent=vaf.in.percent,
                                     num.boots=num.boots,
                                     bootstrap.model=subclonal.test.model,
                                     cluster.center.method=cluster.center,
                                     random.seed=random.seed)
                #bbb <<- boot
            }

            models = enumerate.clones(v, sample=s, variants, boot=boot,
                                      founding.cluster=founding.cluster,
                                      ignore.clusters=ignore.clusters,
                                      min.cluster.vaf=min.cluster.vaf,
                                      p.value.cutoff=sum.p.cutoff,
                                      alpha=alpha)
        }

        if(verbose){cat(s, ':', length(models),
                        'clonal architecture model(s) found\n\n')}
        if (length(models) == 0){
            print(v)
            message(paste('ERROR: No clonal models for sample:', s,
                       '\nCheck data or remove this sample, then re-run.
                       \nAlso check if founding.cluster was set correctly!'))
            return(NULL)
        }else{
            vv[[sample.name]] = models
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
        merged.trees = list()
        merged.traces = NULL
        for (i in 1:num.models){
            scores[i,1] = get.model.score(vv[[1]][[i]])
            merged.trees = c(merged.trees, list(vv[[1]][[i]]))
        }
        scores$model.score = scores[, 1]
    }
    if (nSamples >= 2){
        z = find.matched.models(vv, sample.names, sample.groups,
            merge.similar.samples=merge.similar.samples)
        matched = z$matched.models
        scores = z$scores
        merged.trees = z$merged.trees
        merged.traces = z$merged.traces
        vv = z$models
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
        # sort models by score
        idx = order(scores$model.score, decreasing=T)
        matched = matched[idx, ,drop=F]
        scores = scores[idx, , drop=F]
        merged.trees = merged.trees[idx]
        merged.traces = merged.traces[idx]
    }
    num.matched.models = ifelse(is.null(matched), 0, nrow(matched))
    if (verbose){ cat(paste0('Found ', num.matched.models,
                             ' compatible evolution models\n'))}
    # trim and remove redundant merged.trees
    cat('Pruning merged clonal evolution trees....\n')
    trimmed.merged.trees = trim.clone.trees(merged.trees)$unique.trees
    cat('Number of unique pruned trees:', length(trimmed.merged.trees), '\n')
    results = list(models=vv, matched=list(index=matched,
        merged.trees=merged.trees, merged.traces=merged.traces,
        scores=scores, trimmed.merged.trees=trimmed.merged.trees),
        num.matched.models=num.matched.models)
    cat('Scoring models...\n')
    results = cross.rule.score(results)
    num.sig.models = sum(results$matched$scores$max.clone.ccf.combined.p <= cross.p.cutoff)
    cat(num.sig.models, 'model(s) with p-value <=', cross.p.cutoff, '\n')
    if(num.sig.models == 0){
        message('\n***WARN: Inra-tumor heterogeneity could result in a clone (eg. founding)
         that is not present (no cells) in any samples, although  detectable via
         clonal marker variants due to that its subclones are distinct across
         samples. Therefore, a model with a higher p-value for the CCF of such
         a clone can still be biologically consistent, interpretable, and
         interesting! Manual investigation of those higher p-value models
         is recommended\n\n')
    }

    # record data and params used
    results$variants = variants
    results$params = list(cluster.col.name=cluster.col.name,
                          sum.p.cutoff=sum.p.cutoff,
                          cross.p.cutoff=cross.p.cutoff,
                          alpha=alpha,
                          min.cluster.vaf=min.cluster.vaf,
                          vaf.col.names=vaf.col.names,
                          vaf.in.percent=vaf.in.percent,
                          sample.groups=sample.groups,
                          num.boots=num.boots,
                          bootstrap.model=subclonal.test.model,
                          cluster.center.method=cluster.center,
                          merge.similar.samples=merge.similar.samples,
                          random.seed=random.seed
                          )
    results$boot=boot

    return(results)
}


#' Scale cellular fraction in a clonal architecture model
#'
#' @description Root clone VAF will be determined, and all vaf will be scaled such that
#' roo clone VAF will be 0.5
#'
scale.cell.frac <- function(m, ignore.clusters=NULL){
    #max.vaf = max(m$vaf[!m$excluded & !(m$lab %in% as.character(ignore.clusters))])
    max.vaf = m$vaf[!is.na(m$parent) & m$parent=='-1']
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
#' @param models.to.plot: row numbers of the models to plot in matched$index
#' @param scale.monoclonal.cell.frac: c(TRUE, FALSE); if TRUE, scale cellular
#' fraction in the plots (ie. cell fraction will be scaled by 1/purity =
#' 1/max(VAF*2))
#' @param width: width of the plots (in), if NULL, automatically choose width
#' @param height: height of the plots (in), if NULL, automatically choose height
#' @param out.format: format of the plot files ('png', 'pdf', 'pdf.multi.files')
#' @param resolution: resolution of the PNG plot file
#' @param overwrite.output: if TRUE, overwrite output directory, default=FALSE
#' @param max.num.models.to.plot: max number of models to plot; default = 10
#' @param individual.sample.tree.plot: c(TRUE, FALSE); plot individual sample trees
#' if TRUE, combined graph that preserved sample labels will be produced in graphml
#' output
#' @param merged.tree.plot: Also plot the merged clonal evolution tree across
#' samples
#' @param merged.tree.node.annotation: see plot.tree's node.annotation param;
#' default = 'sample.with.nonzero.cell.frac.ci'
#' @param merged.tree.cell.frac.ci: Show cell fraction CI for samples in merged tree
#' @param tree.node.label.split.character: sometimes the labels of samples are long,
#' so to display nicely many samples annotated at leaf nodes, this parameter
#' specify the character that splits sample names in merged clonal evolution
#' tree, so it will be replaced by line feed to display each sample in a line, 
#' @param tree.node.num.samples.per.line: Number of samples displayed on each line
#' in the tree node; default NULL (do not attempt to distribute samples on lines)
#' @param trimmed.tree.plot: Also plot the trimmed clonal evolution trees across
#' samples in a separate PDF file
#' @param color.node.by.sample.group: color clones by grouping found in sample.group.
#' based on the grouping, clone will be stratified into different groups according
#' to what sample group has the clone signature variants detected. This is useful
#' when analyzing primary, metastasis, etc. samples and we want to color the clones
#' based on if it is primary only, metastasis only, or shared, etc. etc.
#' @param merged.tree.node.size.scale: scale merged tree node size by this value,
#' compared to individual tree node size
#' @param merged.tree.node.text.size.scale: scale merged tree node text size
#' by this value, compared to  individual tree text size
#' @param clone.time.step.scale: scaling factor for distance between the tips
#' of the polygon/bell representing clone
#' @param zero.cell.frac.clone.color: color clone with zero cell fraction
#' in the sample with this color (default = NULL, color using matching color
#' @param zero.cell.frac.clone.border.color: border color of the bell for clones that
#' are not found in sample; if equal "fill", the fill color of bell is used
#' auto-generated)
#' @param cell.frac.ci: Display cell frac CI
#' @param disable.cell.frac: Completely remove cell frac CI info. in plots
#' @param samples: samples to plot (ordered), equal vaf.col.names used in
#' infer.clonal.models
#' @param bell.border.width: border with of bell curve
plot.clonal.models <- function(models, out.dir,
                               matched=NULL,
                               samples=NULL,
                               models.to.plot=NULL,
                               variants=NULL,
                               clone.shape='bell',
                               bell.curve.step=0.25,
                               bell.border.width=1,
                               clone.time.step.scale=1,
                               zero.cell.frac.clone.color=NULL,
                               zero.cell.frac.clone.border.color=NULL,
                               box.plot=FALSE,
                               fancy.boxplot=FALSE,
                               box.plot.text.size=1.5,
                               cluster.col.name = 'cluster',
                               ignore.clusters=NULL,# this param is now deprecated
                               scale.monoclonal.cell.frac=TRUE,
                               adjust.clone.height=TRUE,
                               individual.sample.tree.plot=FALSE,
                               merged.tree.plot=TRUE,
                               merged.tree.node.annotation='sample.with.nonzero.cell.frac.ci',
                               merged.tree.cell.frac.ci=FALSE,
                               trimmed.merged.tree.plot=TRUE,
                               tree.node.label.split.character=',',
                               tree.node.num.samples.per.line=NULL,
                               color.node.by.sample.group=FALSE,
                               color.border.by.sample.group=TRUE,
                               variants.to.highlight=NULL,
                               variant.color='blue',
                               variant.angle=NULL,
                               width=NULL, height=NULL, text.size=1,
                               panel.widths=NULL,
                               panel.heights=NULL,
                               tree.node.shape='circle',
                               tree.node.size = 50,
                               tree.node.text.size=1,
                               merged.tree.node.size.scale=0.5,
                               merged.tree.node.text.size.scale=1,
                               out.format='png', resolution=300,
                               overwrite.output=FALSE,
                               max.num.models.to.plot=10,
                               cell.frac.ci=TRUE,
                               cell.frac.top.out.space=0.75,
                               cell.frac.side.arrow.width=1.5,
                               disable.cell.frac=FALSE,
                               show.score=TRUE,
                               show.matched.index=FALSE,
                               show.time.axis=T,
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
    if (is.null(samples)){
        samples = names(models)
    }
    nSamples = length(samples)
    w = ifelse(is.null(width), 7, width)
    h = ifelse(is.null(height), 3*nSamples, height)
    w2h.scale <<- h/w/nSamples*ifelse(box.plot, 2, 1.5)
    if(box.plot && is.null(variants)){
        box.plot = F
        message('box.plot = TRUE, but variants = NULL. No box plot!')
    }

    if (!is.null(matched$index)){
        scores = matched$scores
        merged.trees = matched$merged.trees
        merged.traces = matched$merged.traces
        trimmed.trees = matched$trimmed.merged.trees
        # for historical reason, use 'matched' variable to indicate index of matches here
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
            pdf(paste0(out.dir, '/', out.prefix, '.pdf'), width=w, height=h,
                useDingbat=F, title='')
        }

        for (i in 1:nrow(matched)){
            if (!is.null(models.to.plot)){
                if (!(i %in% models.to.plot)){next}
            }
            combined.graph = NULL
            this.out.prefix = paste0(out.dir, '/', out.prefix, '-', i)
            if (out.format == 'png'){
                png(paste0(this.out.prefix, '.png'), width=w,
                    height=h, res=resolution, units='in')
            }else if (out.format == 'pdf.multi.files'){
                pdf(paste0(this.out.prefix, '.pdf'), width=w, height=h,
                    useDingbat=F, title='')
            }else if (out.format != 'pdf'){
                stop(paste0('ERROR: output format (', out.format,
                            ') not supported.\n'))
            }

            #num.plot.cols = ifelse(box.plot, 3, 2)
            #num.plot.cols = 2 + box.plot + merged.tree.plot
            num.plot.cols = 1 + box.plot + individual.sample.tree.plot
            par(mfrow=c(nSamples,num.plot.cols), mar=c(0,0,0,0))
            mat = t(matrix(seq(1, nSamples*num.plot.cols), ncol=nSamples))
            if (merged.tree.plot){mat = cbind(mat, rep(nSamples*num.plot.cols+1,nrow(mat)))}
            #print(mat)
            if (is.null(panel.widths)){
                ww = rep(1, num.plot.cols)
                if (merged.tree.plot){ww = c(ww , 1.5)}
                #ww[length(ww)] = 1
                if (box.plot){
                    ww[1] = 1
                }
            }else{
                if (length(panel.widths) != num.plot.cols){
                    stop(paste0('ERROR: panel.widths does not have ',
                                num.plot.cols, ' elements\n'))
                }else{
                    ww = panel.widths
                }
            }

            hh = rep(1, nSamples)
            layout(mat, ww, hh)

            # TODO: Make this ggplot work together with R base plots of polygons
            # and igraph trees
            #var.box.plots = variant.box.plot(var, vaf.col.names = vaf.col.names,
            #                 variant.class.col.name=NULL,
            #                 highlight='is.cancer.gene',
            #                 highlight.note.col.name='gene_name',
            #                 violin=F,
            #                 box=F,
            #                 jitter=T, jitter.shape=1, jitter.color='#80b1d3',
            #                 jitter.size=3,
            #                 jitter.alpha=1,
            #                 jitter.center.method='median',
            #                 jitter.center.size=1.5,
            #                 jitter.center.color='#fb8072',
            #                 display.plot=F)


            for (k in 1:length(samples)){
                s = samples[k]
                s.match.idx = matched[[s]][i]
                m = models[[s]][[matched[[s]][i]]]
                merged.tree = merged.trees[[i]]
                if (scale.monoclonal.cell.frac){
                    #TODO: auto identify ignore.clusters
                    m = scale.cell.frac(m, ignore.clusters=NULL)
                }
                lab = s
                # turn this on to keep track of what model matched
                if (show.matched.index){
                    lab = paste0(s, ' (', s.match.idx, ')')
                }
                if (show.score){
                    lab = paste0(s, '\n(max.p=',
                                 sprintf('%0.3f', scores[[s]][i]), ')')
                }
                top.title = NULL
                if (k == 1 && show.score){
                    top.title = paste0('Max (clone cross-sample p) = ', scores$max.clone.ccf.combined.p[i])
                }
                if (box.plot){
                    current.mar = par()$mar
                    par(mar=c(3,5,3,3))
                    if (fancy.boxplot){
                        # TODO: Make this ggplot work together with R base
                        # plots of polygons
                        # and igraph trees
                        #print(var.box.plots[[i]])
                        stop('ERROR: fancy.plot is not yet available. You
                             can use variant.box.plot function to plot
                             separately\n')
                    }else{
                        with(variants, boxplot(get(s) ~ get(cluster.col.name),
                                           cex.lab=box.plot.text.size,
                                           cex.axis=box.plot.text.size,
                                           cex.main=box.plot.text.size,
                                           cex.sub=box.plot.text.size,
                                           ylab=s))
                    }

                    par(mar=current.mar)
                }
                draw.sample.clones(m, x=2, y=0, wid=30, len=7,
                                   clone.shape=clone.shape,
                                   bell.curve.step=bell.curve.step,
                                   clone.time.step.scale=clone.time.step.scale,
                                   bell.border.width=bell.border.width,
                                   zero.cell.frac.clone.color=zero.cell.frac.clone.color,
                                   zero.cell.frac.clone.border.color=zero.cell.frac.clone.border.color,
                                   label=lab,
                                   text.size=text.size,
                                   cell.frac.ci=cell.frac.ci,
                                   disable.cell.frac=disable.cell.frac,
                                   top.title=top.title,
                                   adjust.clone.height=adjust.clone.height,
                                   cell.frac.top.out.space=cell.frac.top.out.space,
                                   cell.frac.side.arrow.width=cell.frac.side.arrow.width,
                                   variants.to.highlight=variants.to.highlight,
                                   variant.color=variant.color,
                                   variant.angle=variant.angle,
                                   show.time.axis=show.time.axis,
                                   color.node.by.sample.group=color.node.by.sample.group,
                                   color.border.by.sample.group=color.border.by.sample.group)

                if (individual.sample.tree.plot){
                    gs = plot.tree(m, node.shape=tree.node.shape,
                               node.size=tree.node.size,
                               tree.node.text.size=tree.node.text.size,
                               cell.frac.ci=cell.frac.ci,
                               color.node.by.sample.group=color.node.by.sample.group,
                               color.border.by.sample.group=color.border.by.sample.group,
                               show.legend=F,
                               node.prefix.to.add=paste0(s,': '),
                               out.prefix=paste0(this.out.prefix, '__', s))
                }
                
                # plot merged tree
                if (merged.tree.plot && k == nSamples){
                    current.mar = par()$mar
                    par(mar=c(3,3,3,3))

                    # determine colors based on sample grouping
                    node.colors = NULL
                    if ('sample.group' %in% colnames(merged.tree)){
                        node.colors = merged.tree$sample.group.color
                        names(node.colors) = merged.tree$lab
                    }

                    gs2 = plot.tree(merged.tree,
                               node.shape=tree.node.shape,
                               node.size=tree.node.size*merged.tree.node.size.scale,
                               tree.node.text.size=tree.node.text.size*merged.tree.node.text.size.scale,
                               node.annotation=merged.tree.node.annotation,
                               node.label.split.character=tree.node.label.split.character,
                               node.num.samples.per.line=tree.node.num.samples.per.line,
                               cell.frac.ci=merged.tree.cell.frac.ci,
                               #title='\n\n\n\n\n\nmerged\nclonal evolution\ntree\n|\n|\nv',
                               node.prefix.to.add=paste0(s,': '),
                               #node.colors=node.colors,
                               color.node.by.sample.group=color.node.by.sample.group,
                               color.border.by.sample.group=color.border.by.sample.group,
                               out.prefix=paste0(this.out.prefix, '__merged.tree__', s))
                    par(mar=current.mar)
                }

                if (individual.sample.tree.plot){
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
            }
            if (out.format == 'png' || out.format == 'pdf.multi.files'){
                dev.off()
            }
            if (individual.sample.tree.plot){
                write.graph(combined.graph,
                        file=paste0(this.out.prefix, '.graphml'),
                        format='graphml')
            }
        }
        if (out.format == 'pdf'){
            #plot(combined.graph)
            dev.off()
        }
        
        # plot trimmed trees
        if (nrow(trimmed.trees[[1]]) == 0){#single sample, no trimmed merged tree
            trimmed.merged.tree.plot = F
        }
        if (trimmed.merged.tree.plot){
            cat('Plotting trimmed merged trees...\n')
            pdf(paste0(out.dir, '/', out.prefix, '.trimmed-trees.pdf'),
                width=w/num.plot.cols*1.5, height=h/nSamples*7, useDingbat=F, title='')
            for (i in 1:length(trimmed.trees)){
                gs3 = plot.tree(trimmed.trees[[i]],
                           node.shape=tree.node.shape,
                           node.size=tree.node.size*0.5,
                           tree.node.text.size=tree.node.text.size,
                           node.annotation=merged.tree.node.annotation,
                           node.label.split.character=tree.node.label.split.character,
                           node.num.samples.per.line=tree.node.num.samples.per.line,
                           color.border.by.sample.group=color.border.by.sample.group,
                           #cell.frac.ci=cell.frac.ci,
                           cell.frac.ci=F, 
                           node.prefix.to.add=paste0(s,': '),
                           out.prefix=paste0(this.out.prefix, '__trimmed.merged.tree__', s))

            }
            dev.off()

        }

        # write trace file, telling what samples are eliminated from tree
        # due to its similar clonal architecture to another sample
        if (!is.null(merged.traces[[1]])){
            all.traces = NULL
            for (i in 1:nrow(matched)){
                tmp = merged.traces[[i]]
                tmp = cbind(model=i, tmp)
                if (is.null(all.traces)){all.traces = tmp}
                else{all.traces = rbind(all.traces, tmp)}
            }
            all.traces$sample = factor(all.traces$sample, levels=unique(all.traces$sample))
            all.traces = aggregate(similar.sample ~ model+sample, all.traces,
                paste, collapse=',')
            all.traces = all.traces[order(all.traces$model),]
            write.table(all.traces, paste0(out.dir, '/', out.prefix,
                '.sample-reduction.tsv'), sep='\t', row.names=F, quote=F)
        }

    }else{# of !is.null(matched$index); plot all
        # TODO: plot all models for all samples separately.
        # This will serve as a debug tool for end-user when their models
        # from different samples do not match.
        message('No compatible multi-sample models provided.
                Individual sample models will be plotted!\n')
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
        #print(str(v))
        #print(str(clusters))
        #print(str(cl))
        #print(length(vaf.col.names))
        is.one.sample = length(vaf.col.names) == 1
        if (method == 'median'){
            if (is.one.sample){
                median.vafs = median(v[v[[cluster.col.name]]==cl,vaf.col.names])
                names(median.vafs) = vaf.col.names
            }else{
                median.vafs = apply(v[v[[cluster.col.name]]==cl,vaf.col.names],
                                    2, median)
            }
        }else if (method == 'mean'){
            if (is.one.sample){
                median.vafs = mean(v[v[[cluster.col.name]]==cl,vaf.col.names])
                names(median.vafs) = vaf.col.names
            }else{
                median.vafs = apply(v[v[[cluster.col.name]]==cl,vaf.col.names],
                                    2, mean)
            }
        }
        #print(str(median.vafs))
        median.vafs = as.data.frame(t(median.vafs))
        #print(str(median.vafs))
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
    #debug
    #print(base.clusters.idx)
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


# identify if a model is polyclonal (ie. when more than 1 clone
# arose from the normal clone)
is.poly <- function(v){
    res = F
    if (nrow(v[!is.na(v$parent) & v$parent == '0' & !v$excluded,]) > 1){
        res = T
    }
    return(res)
}

# summarize polyclonal models
# return a data frame of polyclonal model status (TRUE if poly)
# the row index of the data frame is the same as row index of
# x$matched$index matrix
sum.polyclonal <- function(x){
    samples = colnames(x$matched$index)
    poly = matrix(rep(NA, length(samples)), nrow=1)
    colnames(poly) = samples
    for (i in 1:nrow(x$matched$index)){
        p = c()
        for (s in samples){
           p = c(p, is.poly(x$models[[s]][[x$matched$index[i, s]]]))
        }
        poly = rbind(poly, p)
    }
    poly = poly[!is.na(poly[,1]),]
    rownames(poly) = seq(1, nrow(poly))
    poly = as.data.frame.matrix(poly)
    return(poly)
}

#' Merge multi region samples into a meta sample
# v1 = x$models$C[[1]]; v2=x$models$M1197[[1]]; mt = x$matched$merged.trees[[1]]
# z = merge.samples(x, 1, c('C', 'M1197'), 'new', 'P', c('C_ref', 'M1197_ref'), c('C_var', 'M1197_var'))
merge.samples <- function(x, samples, new.sample, new.sample.group, ref.cols=NULL, var.cols=NULL){
    if (!all(samples %in% names(x$models))){
        stop('ERROR: Sample not found when merging!')
    }
    
    # merge trees from samples (we need to do this to make sure only clones
    # from the samples' trees are included
    x$models[[new.sample]] = list()
    for (mid in 1:x$num.matched.models){
        #print(mid)
        v = NULL
        for (i in 1:length(samples)){
            s = samples[i]
            vi = x$models[[s]][[x$matched$index[mid,s]]]
            rownames(vi) = vi$lab
            #print(vi[, c( 'lab', 'excluded', 'parent', 'free.mean')])
            if (is.null(v)){
                v = vi
            }else{
                if (nrow(v) != nrow(vi)){
                    stop('ERROR: Number of clusters/clones not equal while merging\n')
                }
                vi = vi[as.character(v$lab),]
                v$vaf = v$vaf + vi$vaf
                vi.only = !vi$excluded & v$excluded
                v$parent[vi.only] = vi$parent[vi.only]
                v$excluded[vi.only] = F
                v$ancestors[vi.only] = vi$ancestors[vi.only]
            }
        }
        v$vaf = v$vaf/length(samples)
        rownames(v) = v$lab
        #print(v[, c( 'lab', 'excluded', 'parent', 'free.mean')])
        x$models[[new.sample]][[mid]] = v
    }

    # update model index in matched
    x$matched$index[[new.sample]] = seq(1,x$num.matched.models)

    # recalculate VAF using read counts combined from all samples
    # if ref and var counts available
    if (!is.null(ref.cols) && !is.null(var.cols)){
        new.ref.col = paste0(new.sample, '.ref')
        new.var.col = paste0(new.sample, '.var')
        new.depth.col = paste0(new.sample, '.depth')
        
        x$variants[[new.ref.col]] = rowSums(x$variants[, ref.cols], na.rm=T)
        x$variants[[new.var.col]] = rowSums(x$variants[, var.cols], na.rm=T)
        x$variants[[new.depth.col]] = x$variants[[new.ref.col]] +
                                        x$variants[[new.var.col]]
        
        # recalc vaf of variants
        x$variants[[new.sample]] = x$variants[[new.var.col]]/x$variants[[new.depth.col]]
        if (x$params$vaf.in.percent){
            x$variants[[new.sample]] = 100*x$variants[[new.sample]]
        }
    }
    
    # recalc center vaf of clusters
    c = estimate.clone.vaf(x$variants, x$params$cluster.col.name,
                            vaf.col.names=new.sample,
                            vaf.in.percent=x$params$vaf.in.percent,
                            method=x$params$cluster.center.method)
    tmp = make.clonal.data.frame(c[[new.sample]], c[[x$params$cluster.col.name]],
        colors=unique(x$models[[1]][[1]]$color))
    rownames(tmp) = tmp$lab
    for (mid in 1:x$num.matched.models){
        tmp = tmp[as.character(x$models[[new.sample]][[mid]]$lab),]
        x$models[[new.sample]][[mid]]$vaf = tmp$vaf
        x$models[[new.sample]][[mid]]$vaf = tmp$free
        x$models[[new.sample]][[mid]]$free.mean = 0
    }

    # update bootstrap samples for the merged sample
    boot = generate.boot(x$variants, vaf.col.names=new.sample,
                        vaf.in.percent=x$params$vaf.in.percent,
                        num.boots=x$params$num.boots,
                        bootstrap.model=x$params$bootstrap.model,
                        cluster.center.method=x$params$cluster.center.method,
                        random.seed=x$params$random.seed)
    # add new sample' bootstraps, and push zero.means down to bottom in boot list
    x$boot[[new.sample]] = boot[[new.sample]]
    tmp = x$boot[['zero.means']]
    x$boot[['zero.means']] = NULL
    x$boot[['zero.means']] = tmp
    tmp = NULL

    # reestimate CCF
    cat('Estimating CCF of clones for merged sample...\n')
    for (mid in 1:x$num.matched.models){
        v =  x$models[[new.sample]][[mid]]
        rownames(v) = v$lab
        for (cl in v$lab[!v$excluded]){
            sub.clusters = v$lab[v$parent == cl & !is.na(v$parent)]
            if(length(sub.clusters) == 0){sub.clusters = NULL}
            v = estimate.ccf(v, new.sample, which(v$lab == cl), x$boot,
                x$params$min.cluster.vaf, x$params$alpha,
                t=NULL, sub.clusters=sub.clusters)
            cat('ccf: ', v[cl,'lab'], paste(sub.clusters, collapse=','), '\n')
        }
        clone.stat = determine.subclone(v, v$lab[!is.na(v$parent)
                                         & v$parent == '-1'])

        v$is.subclone = clone.stat$is.sub[v$lab]
        v$is.founder = clone.stat$is.founder[v$lab]
        v$is.zero = clone.stat$is.zero
        x$models[[new.sample]][[mid]] = v
    }

    # add sample group
    x$params$sample.groups[new.sample] = new.sample.group

    # cleanup: remove models of samples merged to ensure data structure
    # consistency
    x$params$vaf.col.names = c(setdiff(x$params$vaf.col.names, samples), new.sample)
    for (s in samples){
        x$params$sample.groups = x$params$sample.groups[names(x$params$sample.groups) != s]
        x$models[[s]] = NULL
        x$matched$index[[s]] = NULL
    }

    x = merge.all.matched.clone.trees(x)

    return(x)
}

#' Recreate merged trees for matched models, given output of infer.clonal.models
#' 
merge.all.matched.clone.trees <- function(x){
    if (x$num.matched.models == 0 || is.null(x$matched)){
        message('WARN: No matched model to merge.\n')
        return(x)
    }
    samples = names(x$params$sample.groups)
    merged.trees = list()
    merged.traces = list()
    cat('Merging clonal evolution trees across samples...\n')
    x$params$sample.groups = 
    for (i in 1:x$num.matched.models){
        m = list()
        for (j in 1:length(samples)){
            # the order of samples should be the same
            m = c(m, list(x$models[[j]][[x$matched$index[i, j]]]))
        }

        zz = merge.clone.trees(m, samples=samples, x$params$sample.groups,
               merge.similar.samples=x$params$merge.similar.samples)
        mt = zz$merged.tree
        trace = zz$merged.trace
        # after merged, assign sample.group and color to individual tree
        #print(mt)
        for (j in 1:length(samples)){
            x$models[[j]][[x$matched$index[i, j]]] = merge(x$models[[j]][[x$matched$index[i, j]]],
                mt[, c('lab', 'sample.group', 'sample.group.color')], all.x=T)
        }

        # if events already mapped to old merged trees, take over from one of them
        if ('events' %in% colnames(x$matched$merged.trees[[1]])){
            mt = merge(mt, x$matched$merged.trees[[1]][, c('lab', 'events')])
        }

        merged.trees = c(merged.trees, list(mt))
        merged.traces = c(merged.traces, list(trace))
    }
    x$matched$merged.trees = merged.trees
    x$matched$merged.traces = merged.traces

    return(x)   
}

#' Assign events to clones based on presence/absence of them across samples
#' @description: Many events can be clustered correctly due to VAF can not
#' be estimated (eg. indels, copy number, copy altered SNVs). This function
#' try a best guess of what clone should the events belong to. It will look
#' at the "sample" column of the merged tree data frame and match with the
#' samples that we see the event. The clone that have the max number of
#' samples matching will be chosen.
#' @param tree: a merged tree
#' @param events: a data frame containing events and presence/absence
#' status in each sample (required colums: "event", followed by sample
#' columns representing VAF, each in a column
#' @param samples: samples (equivalent to vaf.col.names when calling
#' infer.clonal.models
#' @param cutoff: VAF cutoff to determine presence/absence of an event
#' in a sample
#' 
assign.events.to.clones.of.a.tree <- function(tree, events, samples, cutoff=0){
    rownames(tree) = tree$lab
    if(nrow(tree) == 0 || nrow(events) == 0){return(NULL)}

    # strip off sample note (eg. zero cell frac)
    tree$samples = gsub('°|\\*', '', tree$sample)

    # make binary based on vaf cutoff
    events[, samples][events[, samples] < cutoff] = 0
    events[, samples][events[, samples] >= cutoff] = 1
    events = events[apply(events[, samples] > 0, 1, sum, na.rm=T) > 0,]
    rownames(events) = NULL
    tree$events = ''

    # map events to clones
    # for each event, find the 1st clone that have the max ratio of
    # samples carrying the event
    for (i in 1:nrow(events)){# each event
        e = events[i,samples] > 0
        event.samples = colnames(e)[e[1,]]
        # find clone that shared the most number of samples
        # with the event
        best.match.idx = NULL
        best.match.clone = NULL
        max.match.rate = 0
        max.match.num.samples = 0
        for (j in 1:nrow(tree)){ # each clone
            clone = tree$lab[j]
            clone.samples = unlist(strsplit(tree$samples[j], ','))
            num.match.samples = sum(clone.samples %in% event.samples)
            match.rate = num.match.samples^2/length(clone.samples)
            #match.rate = num.match.samples*length(clone.samples)
            if (match.rate > max.match.rate ){
                max.match.rate = match.rate
                best.match.clone = clone
                best.match.idx = j
            }

        }
        tree$events[best.match.idx] = paste0(tree$events[best.match.idx],
                                                events$event[i], ',')
    }
    tree$events = gsub(',$', '', tree$events)
    tree$samples = NULL
    return(tree)

}

#' Wrapper to assign events to clones in all merged trees
#' This function calls assign.events.to.clones.of.a.tree
#' on each merged tree in list x$matched$merged.trees
#' @param x: output of infer.clonal.models
#' @param events: see assign.events.to.clones.of.a.tree
#' @param samples: see assign.events.to.clones.of.a.tree
#' @param cutoff: see assign.events.to.clones.of.a.tree
assign.events.to.clones <- function(x, events, samples, cutoff=0){
    if (x$num.matched.models > 0){
        for (i in 1:x$num.matched.models){
            x$matched$merged.trees[[i]] = assign.events.to.clones.of.a.tree(
                x$matched$merged.trees[[i]], events, samples, cutoff)
        }
    }
    return(x)
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
get.clonevol.colors <- function(num.colors, strong.color=F){
    colors = c('#a6cee3', '#b2df8a', '#cab2d6', '#fdbf6f', '#fb9a99',
               '#1f78b4','#999999', '#33a02c', '#ff7f00', '#bc80bd',
               '#fca27e', '#ffffb3', '#fccde5', '#fb8072', '#d9d9d9',
               '#f0ecd7', rep('#e5f5f9',10000))
    if(strong.color){
        colors = c('#e41a1c', '#377eb8', '#4daf4a', '#984ea3',
            '#ff7f00', 'black', 'darkgray', rep('lightgray',10000))
        colors[1:3] = c('red', 'blue', 'green')
    }
    if (num.colors > length(colors)){
        stop('ERROR: Not enough colors!\n')
    }else{
        return(colors[1:num.colors])
    }
}

plot.clonevol.colors <- function(num.colors=17){
    library(ggplot2)
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



#'
insert.lf <- function(ss, n, split.char=','){
        if (!is.null(split.char)){
            num.splits = sapply(ss, function(l)
                nchar(gsub(paste0('[^', split.char, ']'), '', l)))

            # only keep the node.label.split.char in interval of node.num.samples.per.line
            # such that a block of node.num.samples.per.line samples will be grouped and
            # kept in one line
            if (!is.null(n)){
                for (i in 1:length(ss)){
                    vl = unlist(strsplit(ss[i], split.char))
                    sel = seq(min(n, length(vl)),
                        length(vl),n)
                    vl[sel] = paste0(vl[sel], split.char)
                    vl[-sel] = paste0(vl[-sel], ';')
                    ss[i] = paste(vl, collapse='')
                }
                num.splits = length(sel) + 1
            }
            extra.lf = sapply(num.splits, function(n) paste(rep('\n', n), collapse=''))
            extra.lf = ''
            ss = paste0(extra.lf, gsub(split.char, '\n',ss))
         }
         return(ss)
}

