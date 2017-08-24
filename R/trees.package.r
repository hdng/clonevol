# This file contains the parital code adapted from the trees package:
# https://github.com/johnbaums/trees
# Various tweaks were made by Ha Dang
# to allow for clonal evolution tree visualization in ClonEvol

#############################################################################
##################### copied from germinate.R file  #########################
#############################################################################
#' Germinate a seed.
#'
#' Observe the miracle of life as your seed germinates and spews forth a glorious
#' tree.
#'
#' @param x Either a \code{seed} object returned by \code{\link{seed}}, or a
#'   named \code{list} containing: \describe{\item{\code{trunk.height}}{Exactly
#'   how glorious will this tree be?} \item{\code{branches}}{A binary coded
#'   vector of addresses of branches to be included. Branches can branch to the
#'   left or the right from the end of the trunk, or from the end of another
#'   branch included in \code{branches}. Elements of \code{branches} may only
#'   contain the characters given in args \code{left} and \code{right}, and all
#'   parent branches of each element must also be provided. E.g. if \code{left}
#'   and \code{right} are \code{'0'} and \code{'1'}, respectively, then
#'   \code{'0'} is the first branch to the left stemming from the top of the
#'   trunk, while \code{'1'} is the first branch to the right, stemming from the
#'   top of the trunk; \code{'01'} would be a branch forking to the right from
#'   the end of the first branch that forked left off the trunk.}
#'   \item{\code{lengths}}{A vector of branch lengths corresponding to the
#'   elements of \code{branches}. Should be the same length as \code{branches}.}}
#' @param angle The angle of branches relative to their parent branch (or
#'   relative to the trunk). This angle is negated for left-pointing branches.
#' @param trunk.width The line width of the trunk. Widths are then scaled down
#'   for successive child branches, to a minimum of 1.
#' @param left The character used to represent left-turning branches in the
#'   \code{branches} vector (default is \code{'0'}). Must not be \code{'Y'}.
#' @param right The character used to represent right-turning branches in the
#'   \code{branches} vector (default is \code{'1'}). Must not be \code{'Y'}.
#' @param plot Should the tree be plotted? (logical).
#' @param ... Further arguments to \code{plot.plant}.
#' @return a \code{plant} object, which is a \code{data.frame} comprising branch
#'   addresses, depths, lengths, angles, base coordinates, and tip coordinates.
#' @seealso \code{\link{seed}} \code{\link{foliate}} \code{\link{prune}}
# @examples
#
# # Motivating example from http://stackoverflow.com/q/28163979/489704.
# # Pass a named list (describing the seed) to germinate.
# germinate(list(trunk.height=32,
#                branches=c('0', '1', '00', '01', '010', '011'),
#                lengths=c(21, 19, 5, 12, 6, 2)),
#            angle=40)
#
# plot(g, trunk.width=15, col='peachpuff4')
#
# Dec.7,2016: code altered by hdng to allow for germinating tree for clonal evolution
#
germinate <- function(x, angle=15, trunk.width=20, middle='0', left='1', right='2',
                      left2 = '3', right2 = '4', left3 = '5', right3 = '6',
                      left4 = '7', right4 = '8',
                      plot=FALSE, ...) {

    #if(is(x, 'seed')) {
    trunk.color = x$branch.colors[1]
    trunk.node.color = x$node.colors[1]
    trunk.node.border.color = x$node.border.colors[1]
    trunk.node.border.width = x$node.border.widths[1]
    trunk.border.color = x$branch.border.colors[1]
    trunk.border.linetype = x$branch.border.linetypes[1]
    trunk.border.width = x$branch.border.widths[1]
    trunk.node.label = x$node.labels[1]
    trunk.node.text = x$node.texts[1]
    trunk.text = x$branch.texts[1]
    x <- list(trunk.height=x$length[1],
              branches=x$branches[-1],
              lengths=x$lengths[-1],
              branch.colors=x$branch.colors[-1],
              branch.border.colors=x$branch.border.colors[-1],
              branch.border.linetypes=x$branch.border.linetypes[-1],
              branch.border.widths=x$branch.border.widths[-1],
              node.colors=x$node.colors[-1],
              node.border.colors=x$node.border.colors[-1],
              node.border.widths=x$node.border.widths[-1],
              node.labels=x$node.labels[-1],
              node.texts=x$node.texts[-1],
              branch.texts=x$branch.texts[-1])
    #}
    if ('Y' %in% c(left, right))
        stop('"Y" is reserved for the trunk.')
    if (any(nchar(c(left, right))) != 1 | left==right)
        stop('left and right must be single, distinct alphanumeric characters.')
    ord = order(x$branches)
    x$lengths <- x$lengths[ord]
    x$branch.colors <- x$branch.colors[ord]
    x$branch.border.colors <- x$branch.border.colors[ord]
    x$branch.border.linetypes <- x$branch.border.linetypes[ord]
    x$branch.border.widths <- x$branch.border.widths[ord]
    x$node.colors <- x$node.colors[ord]
    x$node.border.colors <- x$node.border.colors[ord]
    x$node.border.widths <- x$node.border.widths[ord]
    x$node.labels <- x$node.labels[ord]
    x$node.texts <- x$node.texts[ord]
    x$branch.texts <- x$branch.texts[ord]
    x$branches <- sort(x$branches)
    x$angles <- sapply(sapply(as.character(x$branches), strsplit, ''), function(x) {
        tab <- table(x)
        #sum(c(tab[left]*-angle, tab[right]*angle), na.rm=TRUE)
        sum(c(tab[left]*-angle, tab[right]*angle,
              tab[left2]*-0.5*angle, tab[right2]*angle*0.5,
              tab[left3]*-1.5*angle, tab[right3]*angle*1.5,
              tab[left4]*-2*angle, tab[right4]*angle*2
        ), na.rm=TRUE)
    }, USE.NAMES=FALSE)
    y1 <- x1 <- y0 <- x0 <- rep(NA, length(x$branches))
    for (i in seq_len(length(x$branches))) {
        if(x$branches[i] %in% c(middle, left, right, left2, right2, left3, right3, left4, right4)) {
            x0[i] <- 0
            y0[i] <- x$trunk.height
        } else {
            parent <- substr(x$branches[i], 1, nchar(x$branches[i])-1)
            if (is.null(parent) || length(parent) == 0){
                message('\nERROR: Cannot generate tree, probably due to too many branches (max degree of node = 8)\n')
                return(NULL)
            }
            #cat('dbg:', parent, '---', x$branches[i], '\n')
            x0[i] <- x1[which(x$branches==parent)]
            y0[i] <- y1[which(x$branches==parent)]
        }
        tip <- get.xy(x$angles[i], x$lengths[i], x0[i], y0[i])
        x1[i] <- tip[, 1]
        y1[i] <- tip[, 2]
    }
    d <- data.frame(branches=x$branches,
                    depth=ifelse(x$branches=='Y', 0, nchar(x$branches)),
                    length=x$lengths,
                    angles=x$angles, x0, y0, x1, y1,
                    branch.colors=x$branch.colors,
                    branch.border.colors=x$branch.border.colors,
                    branch.border.linetypes=x$branch.border.linetypes,
                    branch.border.widths=x$branch.border.widths,
                    node.colors=x$node.colors,
                    node.border.colors=x$node.border.colors,
                    node.border.widths=x$node.border.widths,
                    node.labels=x$node.labels,
                    node.texts=x$node.texts,
                    branch.texts=x$branch.texts,
                    stringsAsFactors=FALSE)
    d <- rbind(setNames(
        data.frame('Y', 0, x$trunk.height, 0, 0, 0, 0, x$trunk.height,
                   trunk.color, trunk.border.color, trunk.border.linetype, trunk.border.width,
                   trunk.node.color, trunk.node.border.color, trunk.node.border.width,
                   trunk.node.label, trunk.node.text, trunk.text,
                   stringsAsFactors=FALSE),
        names(d)), d)
    class(d) <- c('plant', 'data.frame')
    if(isTRUE(plot)) plot(d, trunk.width=trunk.width, ...)
    return(d)
}



#############################################################################
##################### copied from plot.plant.R file  #########################
#############################################################################
#' Plot method for plant objects.
#'
#' @param x The \code{plant} object to be plotted.
#' @param trunk.width The plotting line width for the trunk. Successive child
#'   branches are plotted with increasingly finer \code{lwd}, to a minimum of
#'   \code{1}.
#' @param add If \code{TRUE}, the plant will be added to the current plot.
#' @param ... Additional arguments passed to \code{plot} and \code{segments}.
#' @return \code{NULL}
#' @seealso \code{\link{germinate}}
#' @export
#' @param tree.rotation: tree rotation
#' values = c(0, 90, 180) ~ (bottom-up, left-right, top-down)
#' @param text.angle: text angle, if NULL then auto determine
plot.plant <- function(x, trunk.width=20, add=FALSE,
                       tree.rotation=180, text.angle=NULL,
                       branch.width=1, branch.text.size=0.5,
                       node.size=2, node.label.size=0.75,
                       node.text.size=0.5, tree.label=NULL, ...) {

    if(is.null(x)){message('ERROR: No tree to plot\n');return()}

    # distance from node to its label
    l = max(abs(c(x$y0, x$y1)))/30
    r = 1
    #print(tree.rotation)
    if (tree.rotation == 90){
        if (is.null(text.angle)){text.angle=90}
        tmp = x$x0; x$x0=x$y0; x$y0=-tmp
        tmp = x$x1; x$x1=x$y1; x$y1=-tmp
    }else if (tree.rotation == 180){
        r = -1
        x$y0 = r*x$y0; x$y1 = r*x$y1
    }
    if (is.null(text.angle)){text.angle = 0}
    #cat('text.angle = ', text.angle, '\n')

    if(isTRUE(add)) {
        with(x, segments(x0, y0, x1, y1, col=colors, lwd=pmax(trunk.width/nchar(x$branches), 1), ...))
    } else {
        if (tree.rotation == 90){
            plot(c(x$x0-2*r, x$x1+1*r), c(x$y0, x$y1), type='n', asp=1, axes=FALSE, xlab='', ylab='', ...)
        }else{
            plot(c(x$x0, x$x1), c(x$y0, x$y1+2*r), type='n', asp=1, axes=FALSE, xlab='', ylab='', ...)
        }
        #with(x, segments(x0, y0, x1, y1, col=colors, lwd=pmax(trunk.width/nchar(x$branches), 1), ...))
        #with(x, segments(x0, y0, x1, y1, col=branch.colors, lwd=branch.width, ...))
        with(x, draw.branch(x0, y0, x1, y1, fill.color=branch.colors,
                            border.color=branch.border.colors, border.linetype=branch.border.linetypes,
                            border.width=branch.border.widths,
                            w=branch.width, ...))
        with(x, points(x1, y1, pch=21, col=node.border.colors,
                       lwd=node.border.widths,
                       bg=node.colors, cex=node.size, ...))
        with(x, text(x1, y1, labels=node.labels, col='black', cex=node.label.size, srt=text.angle, ...))
        if (tree.rotation == 90){
            with(x, text(x1+l*r, y1, labels=node.texts, col='black', cex=node.text.size,
                         srt=text.angle, ...))
            if(!is.null(tree.label)){
                text(x$x0[1]-2*l, x$y0[1], label=tree.label, cex=node.label.size, srt=text.angle)
            }
        }else{
            with(x, text(x1, y1+l*r, labels=node.texts, col='black', cex=node.text.size,
                         srt=text.angle,...))
        }
        with(x, text((x0+x1)/2, (y0+y1)/2, labels=branch.texts, col='black',
                     cex=branch.text.size, srt=text.angle, ...))
    }
}

#' Draw tree branch using polygon that allows for
#' choosing border style and fill
#' @examples
#' plot(x=c(0,30), y=c(0, 30));
#' \dontrun{
#' draw.branch(c(5, 20, 12, 5, 5, 20),
#'  c(5, 20, 10, 20, 25, 25),
#'  c(10, 12, 20, 12, 10, 20),
#'  c(12, 15, 5, 12, 25, 29),
#'  border.linetype=c(1, 2, 3, 4,1,2),
#'  border.width=c(1,2,1,2,2,2))
#' }
draw.branch <- function(x0, y0, x1, y1, w=1, border.color='black', border.linetype='solid', border.width=0.5, fill.color=NULL, ...){
    X0 = x0; Y0=y0; X1=x1; Y1=y1
    #print(w)
    w = w/2
    #cat('w===', w, '\n')
    #arrows(x0, y0, x1, y1, col='red')
    n = length(x0)
    if (length(border.color) == 1){border.color = rep(border.color, n)}
    if (length(fill.color) == 1){fill.color = rep(fill.color, n)}
    if (length(border.linetype) == 1){border.linetype = rep(border.linetype, n)}
    if (length(border.width) == 1){border.width = rep(border.width, n)}
    for (i in 1:length(X0)){
        x0 = X0[i]; y0 = Y0[i]; x1 = X1[i]; y1 = Y1[i]
        # don't need to catch for atan(Inf), still works
        #if(x0 == x1){
        #    a = pi/2
        #}else{
        a = atan((y1-y0)/(x1-x0))
        #}
        x01 = x0 - sin(a)*w
        x02 = x0 + sin(a)*w
        y01 = y0 + cos(a)*w
        y02 = y0 - cos(a)*w
        x11 = x1 + sin(a)*w
        x12 = x1 - sin(a)*w
        y11 = y1 - cos(a)*w
        y12 = y1 + cos(a)*w
        xx = c(x01, x02, x11, x12, x01)
        yy =  c(y01, y02, y11, y12, y01)
        #print(xx)
        #print(yy)
        polygon(xx, yy, col=fill.color[i], border=border.color[i], lty=border.linetype[i], lwd=border.width[i])

        if (x0 > x1){
            #rect

            # nice
            #polygon(c(x0, x0-b, x1, x1, x1+b, x0, x0), c(y0, y0, y1-b, y1, y1, y0+b, y0), ...)
        }else{
            # nice
            #polygon(c(x0, x0+b, x1, x1, x1-b, x0, x0), c(y0, y0, y1-b, y1, y1, y0+b, y0), ...)
        }
    }

}



#############################################################################
##################### copied from utils.R file  #########################
#############################################################################
# These functions written by Brodie Gaslam (with minor modification by John
# Baumgartner). <http://stackoverflow.com/a/30781090/489704>.

fertilise <- function(size, n) {
    size <- as.integer(size)
    n <- as.integer(n)
    if(size > 25L || size < 3L) stop("Size out of valid range. ",
                                     "'max.depth' must be in the range [3, 25].")

    # Generate integer pool and weights

    size0 <- size - 1L
    pool.raw <- seq.int(2L ^ size) - 1L
    pool.raw.len <- valid.unique <- length(pool.raw)

    # weights are a function of how many trailing zeroes each number has, for
    # example `1000` has three trailing zeroes and represnts `1000`, `100`,
    # `10`, and `1`, so it should be weighed 4x

    weights <- rep(1L, pool.raw.len)
    for(i in seq.int(size0))
        weights[seq.int(from=1L, to=pool.raw.len, by=2 ^ i)] <- i + 1L

    # Create indices to map from the "weighted" vectors to the original
    # vectors

    pool.vals <- rep(pool.raw, weights)
    pool.len <- length(pool.vals)

    # For each repeated value, what count of trailing zeros does it correspond
    # to (equivalent to: `unlist(lapply(weights, seq.int))`, but faster)

    z <- integer(pool.len)
    z[c(1L, cumsum(head(weights, -1L)) + 1L)] <- 1L
    w <- cumsum(!z)
    t <- cummax(z * w)
    zeros.imp <- w - t + 1L
    pad.imp <- weights[pool.vals + 1L] - zeros.imp

    # Generate our encoded vectors by right padding with enough zeros and then
    # adding as a value the number of zeros to the padded area

    zero.pad <- as.integer(2L ^ ceiling(log2(size)))
    vals.enc.init <- pool.vals * zero.pad + zeros.imp - 1L

    # Results tracking

    res <- matrix(0L, nrow=n, ncol=2L)
    res[, 2L] <- size      # leads to "" if not changed
    max.allowed <- size0   # how padded a number to pick can be
    free <- free.init <- rep(TRUE, pool.len)

    # Pre compute frequently used sequences and number patterns

    zero.mx <- as.integer(2 ^ (size - seq(size))) *
        !lower.tri(matrix(ncol=size, nrow=size))
    seqs <- lapply(1L:size, seq.int)
    seqs0 <- lapply(seqs, `-`, 1L)
    seq.rev <- rev(seq.int(size))
    seq.rev0 <- seq.rev - 1L
    ones <- rep(1L, size)
    weights.cs <- cumsum(weights)
    pool.lu <- c(1L, head(weights.cs, -1L) + 1L)

    # Setup our pool to draw blindly (i.e. without checking disqualification);

    blind.i <- 1L
    blind.len <- min(n * 3L, pool.len)        # should most likely get enough values in this draw
    blind.pool <- sample(vals.enc.init, blind.len)

    # Loop through the `n` requested samples

    for(i in seq.int(n)) {
        # Check for completeness, and remove values that would lead to incomplete
        # pools.  We only remove padded values so `valid.unique` is unchanged

        if(max.allowed) {
            rem.pow <- which(n - i >= valid.unique - 2L ^ seqs[[max.allowed]])
            for(j in rev.default(rem.pow)) {
                to.rem <- which(pad.imp == max.allowed)
                free[to.rem] <- FALSE
                max.allowed <- max.allowed - 1L
            }
            if(!max.allowed && n - i >= valid.unique)
                stop(
                    "Logic Error: pool is not large enough to support complete samples"
                ) }
        # Pick from our shuffled pool; after we pick, if turns out value is
        # disqualified, pick again until we find a non-picked value.  Infer from
        # encoding where that value is in our original sorted list

        repeat {
            if(blind.i > blind.len) {  # Initial sample set not enough, add more
                which.free <- which(free)
                if(!length(which.free)) stop("Error: ran out of pool to sample")
                blind.len <- min(length(which.free), blind.len)
                blind.pool <- sample(vals.enc.init[which.free], blind.len)
                blind.i <- 1L
            }
            val.enc <- blind.pool[[blind.i]]
            val <- val.enc %/% zero.pad
            enc <- val.enc %% zero.pad
            blind.idx <- pool.lu[[val + 1L]] + enc
            if(free[[blind.idx]]) break
            blind.i <- blind.i + 1L
        }
        # Figure out how many trailing zeros our number has (recall, this is
        # encoded in the least significant bits of our number); note, zeros is a bit
        # misleading, it means: "how many digits after initial digit are explicilty
        # specied".  The name `zeros` comes from numbers like `1` that would need to
        # add zeros to be specified (e.g. `1000`, which has three zeros)

        weight <- weights[[val + 1L]]
        zeros <- size - weight + enc
        pad <- size0 - zeros
        res[i, ] <- c(val, pad)

        # Based on number of zeros, we can figure out up to what value we need
        # to disqualify (NOTE: different than withbin, here we get the next value
        # greater than our range because `free` is always same size)

        disq.hi.enc <- as.integer((val + 2L ^ pad)) * zero.pad

        # Incremental disqualification of smaller patterns by computing the
        # decimal value from a sequantially truncated bit matrix

        disq.loc.extra <- if(zeros) {
            seq.z <- seqs[[zeros]]
            disqual.more.tmp <- as.integer(
                ones[seq.z] %*% (
                    as.integer(intToBits(val)[seq.rev])[seq.z] *
                        zero.mx[seq.z, seq.z, drop=F]
                ) )
            ws <- weights[disqual.more.tmp + 1L]
            offset <- seqs0[[zeros]] + ws - size
            disq.loc <- pool.lu[disqual.more.tmp + 1L] + offset
            disqualifiable <- which(disqual.more.tmp < val)
            valid.unique <- valid.unique - sum(!(ws - offset - 1L)[disqualifiable])
            unique.default(disq.loc[disqualifiable])
        } else integer()

        # Find values to remove, first with the range of values disqualified by our
        # pick

        lo <- val.enc %/% zero.pad + 1L
        hi <- disq.hi.enc %/% zero.pad + 1L

        free[
            seq.int(
                from=pool.lu[[lo]],
                to=max(pool.lu[[lo]], pool.lu[[hi - 1L]] + weights[[hi - 1L]] - 1L)
            ) ] <- FALSE

        # Now remove any parent values

        free[disq.loc.extra] <- FALSE
        valid.unique <- valid.unique - 2 ^ pad
    }
    # Now convert to binary representation; note we assume ints are 32 bits

    res.raw <- matrix(as.integer(intToBits(res[, 1L])), nrow=32L)[seq.rev, ]
    substr(do.call(paste0, split(res.raw, row(res.raw))), 0L, size - res[, 2L])
}

# Calculate x and y coordinates in a cartesian plane, given distance and angle
# from a given origin.
get.xy <- function(a, d, x0, y0) {
    a <- ifelse(a <= 90, 90 - a, 450 - a)
    data.frame(x = x0 + d * cos(a / 180 * pi),
               y = y0 + d * sin(a / 180 * pi))
}
