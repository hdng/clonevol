# get ggplot2 colors
get.ggplot2.colors <- function(num.colors) {
    hues = seq(15, 375, length=num.colors+1)
    hcl(h=hues, l=65, c=100)[1:num.colors]
}


get.n <- function(x){
    return(c(y = mean(x), label = length(x)))
}

get.median <- function(x){
    return(c(y = median(x), label = as.numeric(sprintf('%0.3f',
        median(x, na.rm=TRUE)))))
}

get.mean <- function(x){
    return(c(y = mean(x), label = as.numeric(sprintf('%0.3f',
        mean(x, na.rm=TRUE)))))
}



# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot
# objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
# e=0.15, # extra height needed for last plot (vertical layout),
# or extra width for first plot (horizontal layout)
multiplot <- function(..., plotlist=NULL, file, cols=1,
                      layout=NULL, horizontal=FALSE, e=0.15) {
    require(grid)

    # Make a list from the ... arguments and plotlist
    plots = c(list(...), plotlist)

    numPlots = length(plots)
    #message(paste0('>>>>>>>INFO: num plots 2 = ', numPlots), '\n')

    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout = matrix(seq(1, cols * ceiling(numPlots/cols)),
                        ncol = cols, nrow = ceiling(numPlots/cols))
    }

    if (numPlots==1) {
        print(plots[[1]])

    } else {

        ## set up heights/widths of plots

        # extra height needed for last plot (vertical layout),
        # or extra width for first plot (horizontal layout)
        hei = rep(1, numPlots)
        # bottom plot is taller
        hei[numPlots] = hei[numPlots]*(1+e)
        wid = rep(1, numPlots)
        # first left plot is wider
        wid[1] = wid[1]*(1+e)
        # Set up the page
        grid.newpage()
        if(horizontal){
            pushViewport(viewport(layout = grid.layout(nrow(layout),
                                                    ncol(layout), widths=wid)))
        }else{
            pushViewport(viewport(layout = grid.layout(nrow(layout),
                                                    ncol(layout), heights=hei)))

        }

        # Make each plot, in the correct location
        for (i in 1:numPlots) {
            # Get i,j matrix positions of the regions containing this subplot
            matchidx = as.data.frame(which(layout == i, arr.ind = TRUE))
            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
        }
    }
}



#' Replace NA with some value
replaceNA <- function(x, by=0){
    idx = is.na(x)
    if (any(idx)){
        cat('WARN: Some NA replaced!\n')
        x[idx] = by
    }
    return(x)
}

#' Replace values bigger than a cutoff by a value
cutBigValue <- function(x, maxValue){
    idx = x > maxValue
    if (any(idx)){
        x[idx] = maxValue
    }
    return(x)
}


# Create opposite hjust values for variants that have VAF
# close (ie. next) to each other
randomizeHjust <- function(df.hi, cluster.col.name='cluster',
                           vaf.name, hjust=0.5){
    df.hi$vaf = df.hi[[vaf.name]]
    df.hi$cluster = df.hi[[cluster.col.name]]
    df.hi = df.hi[with(df.hi, order(cluster, vaf)),]
    df.hi$hjust = 0
    df.hi$newX = df.hi[[cluster.col.name]]
    for (c in unique(df.hi$cluster)){
        x = df.hi[df.hi$cluster == c,]
        for (i in 1:nrow(x)){
            x[i,]$hjust = 1 - (hjust*2*(i%%2-0.5) + 0.5)
            x[i,]$newX = x[i,]$newX - hjust*(i%%2-0.5)/6
        }
        df.hi[df.hi$cluster == c,] = x
    }
    #print(df.hi$hjust)
    return(df.hi)
}

#backward compatible for variant.box.plot
#' @export variant.box.plot
variant.box.plot <- function(...){
    return(plot.variant.clusters(...))
}

#' Plot variant clustering using combination of box, violin, and jitter plots
#' @description Plot variant clustering using combination of box, violin,
#' and jitter plots. Value in columns vaf.col.names variants data.frame, grouped
#' group by cluster.col.name column
#' @usage variant.box.plot(df, cluster.col.name, vaf.col.names, ...)
#' @param df Variants data frame, containing at least 'cluster' and VAF or CCF
#' columns
#' @param cluster.col.name Name of cluster column in df (default='cluster')
#' @param vaf.col.names Names of VAF columns in df
#' @param horizontal if TRUE, samples are laid out horizontally,
#' otherwise vertically
#' @param show.cluster.size show cluster size in the box, jitter, violin plots
#' @param jitter.center.display.value: display cluster center value in the box
#'  when show.cluster.size == FALSE
#' @param jitter.center.display.value.text.size: text size of the value annotated
#' for the center of the box, jitter, violin plots
#' @param width,height,w1,h1 (w/h = with/height of whole plot,
#' w1/h1 = width/height of component plot)
#' w1/h1 will orverwite w/h if they are non-zero. If all w/h/w1/h1 are
#' zero (default), auto scale the plots
#' @param vaf.suffix Suffix to add to vaf.col.names and display in plot axis
#' @param show.cluster.label Show cluster label in axis (set to FALSE to get
#' more space when too many samples are aligned)
#' @param panel.border.colors Colors of panel.border of sample plot
#' @param panel.border.sizes Line sizes of panel.border of sample plot
#' @param panel.border.linetypes Line types of panel.border of sample plot
#' @param vaf.limits Vector (or scalar) of max VAF in
#' @param variant.class.col.name Name of the column containing classification of
#' variants
#' @param show.cluster.size Show cluster size (default=FALSE)
#' @param cluster.size.text.color Cluster size text color (default='blue')
#' @param cluster.axis.name Name to place before cluster IDs in the cluster axis
#' (default = 'cluster:')
#' @param show.cluster.axis.label Show label of the cluster axis (deafault=TRUE)
#' @param sample.title.size Text size of the label of the sample axis (default
#' defined relative to param base_size)
#' @param cluster.title.size Text size of the label of the cluster axis (default
#' defined relative to param base_size)
#' @param base_size base_size Paramter to passed to in ggplot2 theme_bw(base_size)
#' @param axis.ticks.length Length of the axes' ticks (default=1)
#' @param axis.text.angle Angle of the axis text label (default=0)
#' @param plot.margin Margin of the individual sample plot (default=0.1)
#' @param box Plot the box plot (default = TRUE)
#' @param box.line.type Linetype of the boxplot (eg. 'solid', 'dotted',...)
#' @param box.line.size Size of the border/line of the boxplot (default=0.5)
#' @param box.outlier.shape Shape of outliers of boxplot (default=1)
#' @param box.alpha Alpha (transparency) level of the box (default=0.5)
#' @param violin Plot violin plot (default=TRUE)
#' @param violin.line.type Line type for violin plot (default='dotted')
#' @param violin.line.size Linesize for violin plot (default=0.5)
#' @param violin.fill.color Fill color for violin plot (default='grey80')
#' @param violin.alpha Alpha (transparency) level of the violin (default=0.5)
#' @param jitter Plot jitter plot (default=FALSE)
#' @param jitter.width Width of the jitter plot (default=0.5)
#' @param jitter.color Jitter color (default='lightblue')
#' @param jitter.alpha Alpha level of jitter points (default=0.5)
#' @param jitter.size Size of jitter points (default=1)
#' @param jitter.shape Shape of jitter points (default=3)
#' @param jitter.center.method Method to calculate the center of jitter plot
#' (default='median', can also be 'median')
#' @param jitter.center.color Color of the line placed an the center of the jitter
#' plot (default='black')
#' @param jitter.center.size Size of the jitter center line (default=1)
#' @param jitter.center.linetype Line type of the jitter center line
#' (default='solid')
#' @param jitter.center.display.value String indicating the value to show at the
#' center of the jitter points (default='none', can also be 'mean', 'median')
#' @param jitter.center.display.value.text.size Size of text showing jitter
#' center value (default=5)
#' @param highlight A TRUE/FALSE index vector of the variants in df to higlight
#' and plot using a different color (default=NULL)
#' @param highlight.color Color of the higlighted variants in jitter
#' (default='darkgray')
#' @param highlight.fill.color Fill color of the jitter points (default='red'),
#' only highlight.shape a is fillable R shape
#' @param highlight.shape Shape of the highlighted variant jitter points
#' (default=21)
#' @param highlight.size Size of the highligted variant jitter points (default=1)
#' @param highlight.color.col.name Name of a column of colors used to highlight
#' the border of the variant points, to be used when we want to highlight different
#' driver variants using different colors, eg.
#' APC by red, TP53 variants by orange, etc. (default=NULL)
#' @param highlight.fill.col.names Names of columns that store the values to be filled
#' in variant jitter points (when its R shape is fillable). An example case is to
#' fill points with copy number values stored in highlight.fill.col.names for
#' in individual samples. The length(highlight.fill.col.names) must be
#' equal length(vaf.col.names). (default=NULL).
#' @param highlight.fill.min Lower bound of fill value highlight.fill.col.names
#' for gradient fill. Any value lower than this will be fill using the same
#' color (default=1)
#' @param highlight.fill.mid Center value in highlight.fill.col.names columns
#' for generating fill gradient (default=2)
#' @param highlight.fill.max Higher bound of fill value highlight.fill.col.names
#' for gradient fill. Any value higher than this will be fill using the same
#' color (default=3)
#' @param highlight.fill.min.color Low value fill color (default='green'), used
#' in colorpanel(low, mid, hi)
#' @param highlight.fill.mid.color Median value fill color (default='black'), used
#' in colorpanel(low, mid, hi)
#' @param highlight.fill.max.color High value fill color (default='red'), used
#' in colorpanel(low, mid, hi)
#' @param highlight.size.names Column names of the value to be used to scale the
#' size of the variant jitters (default=NULL). Example use case is to draw variants
#' with higher depth bigger.
#' @param max.highlight.size.value Max value to provide to size scaling, higher
#' than this does not make the size bigger (default=500)
#' @param size.breaks Values used to break sizes to show in legends (default=NULL)
#' @param highlight.size.legend.title Title of the size legend (default='depth')
#' @param highlight.note.col.name Column name of the label of the highligthed
#' variants (default=NULL), can use gene name, variant details, etc.
#' @param highlight.note.color Color used for highligted note (default='blue')
#' @param highlight.note.size Size of highlighted note text (default= 3)
#' @param ordered.x The order that the clusters are plot (default=NULL)
#' @param order.by.total.vaf Order clusters by their total VAF across samples
#' (default=TRUE)
#' @param display.plot Show plot (default=TRUE)
#' @param ccf Plot CCF (calculated as 2xVAF) instead of VAF (default=FALSE)
#' @param founding.cluster Founding cluster (default=NULL)
#' @param show.cluster.label Show cluster label in axis (default=TRUE)
#' @param highlight.note.angle Angle of the text annotation for higlighted
#' variants (default = NULL => auto)
#' @import ggplot2
#' @export plot.variant.clusters
#' @examples
#' data(aml1)
#' pp <- plot.variant.clusters(aml1$variants,
#'          cluster.col.name = 'cluster',
#'          show.cluster.size = FALSE,
#'          cluster.size.text.color = 'blue',
#'          vaf.col.names = aml1$params$vaf.col.names,
#'          sample.title.size = 20,
#'          violin = FALSE,
#'          box = FALSE,
#'          jitter = TRUE,
#'          jitter.color=c('#999793', '#8d4891', '#f8e356',
#'                      '#fe9536', '#d7352e'),
#'          display.plot=FALSE)
#'
plot.variant.clusters <- function(df,
                             cluster.col.name='cluster',
                             vaf.col.names=NULL,
                             panel.border.colors='black',
                             panel.border.sizes=1,
                             panel.border.linetypes='solid',
                             vaf.suffix='',
                             vaf.limits=100,
                             variant.class.col.name=NULL,
                             show.cluster.size2=FALSE,
                             show.cluster.size=FALSE,
                             cluster.size.text.color='blue',
                             cluster.axis.name='cluster:',
                             show.cluster.axis.label=TRUE,

                             sample.title.size=NULL,
                             cluster.title.size=NULL,
                             #panel.border.linetype='solid',
                             #panel.border.linesize=1,
                             base_size=18, width=0, height=0,
                             width1=0, height1=0, hscale=1, vscale=1,
                             axis.ticks.length=1,
                             axis.text.angle=0,
                             plot.margin=0.1,
                             horizontal=FALSE,

                             box=TRUE,
                             box.line.type = 'solid',
                             box.line.size=0.5,
                             box.outlier.shape=1,
                             box.alpha=0.5,
                             violin=TRUE,
                             violin.line.type = 'dotted',
                             violin.line.size=0.5,
                             violin.fill.color='grey80',
                             violin.alpha=0.5,
                             jitter=FALSE,
                             jitter.width=0.5,
                             jitter.color='lightblue',
                             jitter.alpha=0.5,
                             jitter.size=1,
                             jitter.shape=1,
                             jitter.center.method='median',
                             jitter.center.color='black',
                             jitter.center.size=1,
                             jitter.center.linetype='solid',
                             jitter.center.display.value='none', #'mean', 'median'
                             jitter.center.display.value.text.size=5,

                             highlight=NULL,
                             highlight.color='blue',
                             highlight.fill.color='red',
                             highlight.shape=21,
                             highlight.size=1,
                             highlight.color.col.name=NULL,
                             highlight.fill.col.names=NULL,
                             highlight.fill.min=1,
                             highlight.fill.mid=2,
                             highlight.fill.max=3,
                             highlight.fill.min.color='green',
                             highlight.fill.mid.color='black',
                             highlight.fill.max.color='red',
                             highlight.size.names=NULL,
                             max.highlight.size.value=500,
                             size.breaks=NULL,
                             highlight.size.legend.title='depth',
                             highlight.note.col.name = NULL,
                             highlight.note.color = 'blue',
                             highlight.note.size = 3,
                             highlight.note.angle = NULL,

                             ordered.x = NULL,
                             order.by.total.vaf=TRUE,

                             display.plot=TRUE,
                             ccf=FALSE, # cancer cell fraction plot
                             founding.cluster=NULL,
                             show.cluster.label=TRUE

){
    library(ggplot2)
    library(grid)
    library(gridExtra)

    # make sure factor is converted to string first to avoid factor cluster
    # being treated as number
    df[[cluster.col.name]] = as.character(df[[cluster.col.name]])
    founding.cluster = as.character(founding.cluster)

    # prepare panel.border params
    n.samples = length(vaf.col.names)
    if (length(panel.border.colors) != n.samples){
        panel.border.colors = rep(panel.border.colors[1], n.samples)
    }

    if (length(panel.border.sizes) != n.samples){
        panel.border.sizes = rep(panel.border.sizes[1], n.samples)
    }
    if (length(panel.border.linetypes) != n.samples){
        panel.border.linetypes = rep(panel.border.linetypes[1], n.samples)
    }

    # scale ccf
    if (ccf){
        founding.vaf = colMeans(df[df[[cluster.col.name]]==founding.cluster, vaf.col.names])
        ccf.scale = t(as.matrix(100/founding.vaf))
        ccf.scale = ccf.scale[rep(1,nrow(df)),]
        df[, vaf.col.names] = df[, vaf.col.names]*ccf.scale
        vaf.limits = 2*vaf.limits
    }

    # order variants by decreasing total vafs
    if (!is.null(ordered.x)){
        cluster.orders = seq(1,length(ordered.x))
        names(cluster.orders) = ordered.x
        df = df[order(cluster.orders[df$cluster]),]
    }else if (order.by.total.vaf){
        #TODO: there is a potential bug here, cluster.col.name should be used!!!
        # fixed, but haven't tested
        df$total.vaf = apply(df[,vaf.col.names], 1, sum)
        df1 = df; df1$cluster = df1[[cluster.col.name]]
        mean.total.vafs = aggregate(total.vaf ~ cluster, df1, mean)
        mean.total.vafs = mean.total.vafs[order(mean.total.vafs$total.vaf,
                                                decreasing=TRUE),]
        rownames(mean.total.vafs) = mean.total.vafs$cluster
        mean.total.vafs.names = mean.total.vafs$cluster
        mean.total.vafs = mean.total.vafs$total.vaf
        names(mean.total.vafs) = mean.total.vafs.names
        df = df[order(mean.total.vafs[df[[cluster.col.name]]], decreasing=TRUE),]
    }


    # change cluster id to continous values to enable adjustment of postions
    # of highlighed genes
    cluster.levels = unique(df[[cluster.col.name]])
    df$cluster.label = as.character(df[[cluster.col.name]])
    df[[cluster.col.name]] = as.integer(factor(df[[cluster.col.name]],
                                               levels=cluster.levels))

    clusters = unique(df[[cluster.col.name]])
    x.axis.breaks = clusters
    if (cluster.axis.name != ''){
        x.axis.breaks = c(0, clusters)
    }
    cluster.labels = unique(df$cluster.label)


    nPlots = length(vaf.col.names)
    plots = list()
    plotCnt = 0
    clusterSizes = table(df[[cluster.col.name]])

    sumCnts = NULL
    if (!is.null(variant.class.col.name)){
        sumCnts = as.data.frame.matrix(table(df[[cluster.col.name]],
                                             df[[variant.class.col.name]]))
    }
    #sumCnts$total = apply(sumCnts, 1, sum)
    nClusters = length(clusterSizes)
    #cat('Number of clusters:', nClusters, '\n', sep='')
    boxColor = 'black'
    if (length(vaf.limits) == 1){
        vaf.limits = rep(vaf.limits, length(vaf.col.names))
    }
    for (ii in 1:length(vaf.col.names)){
        yName = vaf.col.names[ii]
        size.col.name = NULL
        fill.col.name = NULL
        # get color and size columns for highlight variants
        if (!is.null(highlight.size.names)){
            size.col.name = highlight.size.names[ii]
        }
        highlight.fill.col.name=NULL
        if (!is.null(highlight.fill.col.names) && length(highlight.fill.col.names) > 0){
            highlight.fill.col.name = highlight.fill.col.names[ii]
            df[[highlight.fill.col.name]] = cutBigValue(df[[highlight.fill.col.name]],
                highlight.fill.max)
        }


        # since violin plot will throw an error if there is zero variance
        # we'll add a very small number to first value of each cluster
        for (cl in clusters){
            df[df[[cluster.col.name]] == cl,][[yName]][1] =
                df[df[[cluster.col.name]] == cl,][[yName]][1] + 0.0001
        }


        plotCnt = plotCnt + 1
        p = ggplot(data=df, aes_string(x = cluster.col.name, y = yName,
                                       group=cluster.col.name))


        if (jitter){
            if (is.null(jitter.color)){jitter.color='lightblue'}
            if (length(jitter.color) > 1){
                names(jitter.color) = cluster.labels
                #print(jitter.color)
                p = (p + geom_jitter(aes(color=cluster.label),
                                     height = 0,
                                     width=jitter.width,
                                     size=jitter.size,
                                     alpha=jitter.alpha, shape=jitter.shape,
                                     stroke=jitter.size,
                                     )
                     + scale_color_manual(values=jitter.color, guide=FALSE)
                )
            }else{
                p = p + geom_jitter(height = 0, color=jitter.color,
                                    size=jitter.size,
                                    alpha=jitter.alpha, shape=jitter.shape,
                                    stroke=jitter.size,
                                    width=jitter.width)
            }

            # mean or median
            if (jitter.center.method %in% c('median', 'mean')){
#                p = p + geom_errorbar(#stat = "hline",
#                                  yintercept = jitter.center.method,
#                                  width=0.8, size=jitter.center.size,
#                                  linetype=jitter.center.linetype,
#                                  color=jitter.center.color,
#                                  aes(ymax=..y..,ymin=..y..))
                p = p + stat_summary(fun.y=jitter.center.method,
                    aes(ymin=..y.., ymax=..y..), geom='errorbar',
                    width=0.5, size=jitter.center.size,
                    linetype=jitter.center.linetype,
                    color=jitter.center.color)

            }

        }

        if (box && violin){
            p = (p + geom_violin(scale='width', color=boxColor,
                                 fill = violin.fill.color,
                                 alpha=violin.alpha,
                                 linetype=violin.line.type,
                                 size=violin.line.size)
                 + geom_boxplot(color=boxColor, width=0.25,
                                linetype=box.line.type, size=box.line.size,
                                outlier.shape=box.outlier.shape,
                                alpha=box.alpha)
            )

        }else if(box){
            p = p + geom_boxplot(color=boxColor,
                                outlier.shape=box.outlier.shape,
                                alpha=box.alpha)
        }else if(violin){
            p = p + geom_violin(scale='width', color=boxColor,
                                fill = violin.fill.color,
                                alpha=violin.alpha,
                                linetype=violin.line.type,
                                size=violin.line.size)
        }else if (!jitter){
            stop('Must specify at least boxplot, violin, or jitter plot\n')
        }

        if (!is.null(highlight)){

            df.hi = df[df[[highlight]],]

            if (nrow(df.hi) > 0){
                df.hi = randomizeHjust(df.hi, cluster.col.name=cluster.col.name,
                                       vaf.name=yName, hjust=0.7)
                if (!is.null(size.col.name)){
                    df.hi[[size.col.name]] = cutBigValue(df.hi[[size.col.name]],
                                                    max.highlight.size.value)
                }
                #df.hi$note = paste0(df.hi$gene_name, '\n(',
                #    df.hi$amino_acid_change, ')')
                if (!is.null(highlight.fill.col.name)){
                    p = p + geom_point(data=df.hi,
                                       aes_string(x = 'newX', y=yName,
                                                  #size=size.col.name,
                                                  fill=highlight.fill.col.name),
                                                  size=highlight.size,
                                                  shape=highlight.shape,
                                                  color=highlight.color)
                    p = p + scale_fill_gradient2(low=highlight.fill.min.color,
                        mid=highlight.fill.mid.color, high='red', midpoint=2,
                        limits=c(highlight.fill.min, highlight.fill.max),
                        guide=FALSE)
                }else if(!is.null(size.col.name)){
                    p = p + geom_point(data=df.hi,
                                       aes_string(x = 'newX', y=yName,
                                                  size=size.col.name),
                                       color=highlight.color,
                                       shape=highlight.shape,
                                       fill=highlight.fill.color,
                                       show_guide=TRUE)

                }else{
                    p = p + geom_point(data=df.hi,
                                       aes_string(x = 'newX', y=yName),
                                       size=highlight.size,
                                       color=highlight.color,
                                       shape=highlight.shape,
                                       fill=highlight.fill.color,
                                       stroke=highlight.size/5,
                                       #stroke=0.01,
                                       show_guide=TRUE)
                }
                if (!is.null(size.col.name)){
                    #size.breaks = c(0, 50, 100, 200, 300, 500)
                    if (is.null(size.breaks)){
                        size.breaks = seq(0, max.highlight.size.value, max.highlight.size.value/4)
                    }
                    #size.breaks = c(0,1,2,3,4)
                    size.breaks = size.breaks[size.breaks <=
                                                  max.highlight.size.value]
                    size.labels = size.breaks
                    size.labels[length(size.labels)] =
                        paste0('>', size.labels[length(size.labels)])
                    p = (p + scale_radius(name=highlight.size.legend.title,
                                               limits=c(0,max.highlight.size.value),
                                               breaks=size.breaks,
                                               labels=size.labels,
                                               #max_size=1
                                               )
                         + theme(legend.position=c(0.7,0.9))
                         #+ theme(legend.position=c(0.5,0.5))
                    )

                }
                if (!is.null(highlight.note.col.name)){
                    if (is.null(highlight.note.angle)){
                        highlight.note.angle = 0
                        if (horizontal){highlight.note.angle = 45}
                    }
                    p = p + geom_text(data=df.hi,
                                      aes_string(x=cluster.col.name, y=yName,
                                                 label=highlight.note.col.name,
                                                 angle=highlight.note.angle,
                                                 hjust='hjust'),
                                      size=highlight.note.size,
                                      color=highlight.note.color)
                }

            }

        }

        p = (
            p + scale_y_continuous(limits = c(0,vaf.limits[plotCnt]))
            + theme_bw(base_size=base_size)
            + theme(panel.border=element_rect(linetype=panel.border.linetypes[ii],
                                              size=panel.border.sizes[ii],
                                              color=panel.border.colors[ii]))
            + theme(plot.margin = unit(x = c(plot.margin, plot.margin,
                                             plot.margin, plot.margin),
                                       units = "mm"))
            + theme(axis.ticks.length = unit(axis.ticks.length, units = "mm"))
            #+ theme(axis.text.x = element_text(angle=axis.text.angle, hjust=1))
            #+ theme(axis.text.y = element_text(angle=axis.text.angle, hjust=1))
        )
        #show.cluster.label = FALSE
        if (!show.cluster.label){
            p = p + theme(axis.text.x=element_blank())
        }
        if (show.cluster.size){
            p = p + stat_summary(fun.data = get.n, geom = "text",
                                 position = position_dodge(#height = 0,
                                                           width = 0.75),
                                 size = 5, color=cluster.size.text.color)
        } else {
            if (jitter.center.display.value == 'mean'){
                 p = p + stat_summary(fun.data = get.mean, geom = "text",
                                     position = position_dodge(#height = 0,
                                                               width = 0.75),
                                     size = jitter.center.display.value.text.size,
                                     color=cluster.size.text.color)

            } else if (jitter.center.display.value == 'median'){
                 p = p + stat_summary(fun.data = get.median, geom = "text",
                                     position = position_dodge(#height = 0,
                                                               width = 0.75),
                                     size = jitter.center.display.value.text.size,
                                     color=cluster.size.text.color)

            }
        }

        if (vaf.suffix != ''){p = p + ylab(paste0(yName, vaf.suffix))}

        if (horizontal){
            if (plotCnt > 1){
                #p = p + theme(axis.title.y = element_blank())
                p = p + xlab(NULL)
                p = p + scale_x_continuous(breaks = seq(1,nClusters),
                                           #labels=clusterSizes,
                                           limits=c(0, nClusters+1))

            }else{
                p = p + scale_x_continuous(breaks = seq(1,nClusters),
                                           #labels=clusterSizes,
                                           limits=c(0, nClusters+1))
                if (!show.cluster.axis.label){ p = p + xlab(NULL)}
            }
            p = p + coord_flip()
        }else{
            if (plotCnt < nPlots){
                labs = cluster.labels
                if (cluster.axis.name != ''){
                    labs = c(cluster.axis.name, cluster.labels)
                }
                p = (p + theme(axis.title.x = element_blank())
                     + scale_x_continuous(breaks = x.axis.breaks,
                                          labels=labs,
                                          limits=c(0.4, nClusters+0.75))
                )
            }else{
                x.title = cluster.col.name
                if (!is.null(sumCnts)){
                    if (cluster.axis.name == ''){
                        message('Empty cluster.axis.name parameter was reset!')
                        cluster.axis.name = 'cluster:'
                    }
                    z = sumCnts
                    zNames = colnames(z)
                    z$summary = apply(z, 1, paste, collapse="\n")
                    strSummary = paste(apply(sumCnts, 1, sum), '--',
                                       z$summary, sep='\n')
                    labs = paste(rownames(z), '\n--\n', strSummary , sep='')
                    labs = paste(cluster.labels, '\n--\n', strSummary,
                                 sep='')
                    zName = paste(zNames, collapse=":\n", sep='')
                    labs = c(paste(paste0(cluster.axis.name,
                                        '\n--\ntotal:\n--\n'), zName, ':',
                                   sep=''), labs)
                    x.title = paste(cluster.col.name,'(w/ sum of ',
                                    variant.class.col.name, ')', sep='')
                }else{
                    labs = cluster.labels
                    if (cluster.axis.name != ''){
                        labs = c(cluster.axis.name, cluster.labels)
                    }
                }

                #if (!show.cluster.axis.label){
                #    p = p + xlab(NULL)
                #    p = p + scale_x_continuous(breaks = x.axis.breaks,
                #                                labels=labs,
                #                                limits=c(0.4, nClusters+0.75))

                #}else{

                    p = p + scale_x_continuous(x.title,
                                               breaks = x.axis.breaks,
                                               labels=labs,
                                               limits=c(0.4, nClusters+0.75))
                #}

            }

        }
        if (!is.null(sample.title.size)){
            if(horizontal){
                p = p + theme(axis.title.x = element_text(size=sample.title.size))
            }else{
                p = p + theme(axis.title.y = element_text(size=sample.title.size))
            }
        }
        if (!is.null(cluster.title.size)){
            if(horizontal){
                p = p + theme(axis.title.y = element_text(size=cluster.title.size))
            }else{
                p = p + theme(axis.title.x = element_text(size=cluster.title.size))
            }
        }
        if (!show.cluster.axis.label) {p = p + theme(axis.title.x=element_blank())}



        plots = c(plots, list(p))
    }

    # adjust width and height of plot if not given
    w = width
    h = height
    w1 = width1
    h1 = height1
    if (w1 > 0 & horizontal){ w = w1*length(vaf.col.names)}
    if (h1 > 0 & !horizontal){ h = h1*length(vaf.col.names)}
    e = ifelse(is.null(sumCnts), 0.1, 0.125*(ncol(sumCnts) + 4))
    #print(e)
    if ((w == 0 | h == 0) & (w1 == 0 | h1 == 0))
    {
        w = 1+0.5*nClusters
        h = 2.5*length(vaf.col.names) + 2*e
        if (horizontal){
            w = 15
            h = 5
        }
    }

    if (violin){
        #w = 2*w
    }

    w = w*hscale
    h = h*vscale
    #message('\nINFO:num. plots = ', nPlots, '\n')

    if (horizontal){
        if (display.plot){
            multiplot(plotlist=plots, cols=nPlots, horizontal=TRUE)
        }
    }else{
        if (display.plot){
            multiplot(plotlist=plots, cols=1, horizontal=FALSE, e=e)
        }
    }
    return(plots)
}

#' Plot the mean/median of the clusters of variants across samples
#' @description Plot the mean or median of the clusters VAF across
#' samples in a single plot
#' @usage plot.cluster.flow(variants, cluster.col.name, vaf.col.names, ...)
#' @param variants Variant data frame
#' @param cluster.col.name Name of cluster column (default='cluster')
#' @param sample.names Names of samples, corresponding to vaf.col.names
#' (default=NULL)
#' @param vaf.in.percent VAF is in percent (default=TRUE)
#' @param center.measure Method used to determine the center of VAFs of variants
#' within a cluster (default='median', can also be 'mean')
#' @param x.title Title of x axis (default="Variant Allele Frequency (\%)")
#' @param y.title Title of y axis
#' @param line.size Size of lines (default=1)
#' @param colors Colors of the clusters' variant data points
#' @param shapes Shapes of the clusters' variant data points
#' @param width Width of the output file
#' @param height Height of the output file
#' @param out.file Output file (can be pdf, png, etc.) (default=NULL). If equal
#' NULL, this function return the plot that can be print
#' @import ggplot2
#' @export plot.cluster.flow
#' @examples
#' data(aml1)
#' plot.cluster.flow(aml1$variants,
#'                   vaf.col.names = aml1$params$vaf.col.names,
#'                   sample.names = c('Primary', 'Relapse'),
#'                   out.file = paste0(tempfile(), '.pdf'),
#'                   colors = c('#999793', '#8d4891', '#f8e356',
#'                      '#fe9536', '#d7352e'))
#'
plot.cluster.flow <- function(variants,
                              cluster.col.name='cluster',
                              ignore.clusters=NULL,
                              vaf.col.names=NULL,
                              sample.names=NULL,
                              vaf.in.percent=TRUE,
                              center.measure='median',
                              low.vaf.no.line=FALSE,
                              min.cluster.vaf=0,
                              line.size=1,
                              shape.size=5,
                              colors=NULL,
                              shapes=NULL,
                              x.title=NULL,
                              y.title='Variant Allele Frequency (%)',
                              out.file=NULL,
                              width=7,
                              height=5){
    library(grid)
    library(gridExtra)
    library(ggplot2)
    library(reshape2)

    var = variants
    #var[[cluster.col.name]] = as.character(var[[cluster.col.name]])
    cluster.names = unique(var[[cluster.col.name]])
    sorted.cluster.names = cluster.names[order(as.integer(cluster.names))]
    num.clusters = length(sorted.cluster.names)
    if (is.null(colors)){
        #colors = get.ggplot2.colors(num.clusters)
        colors = get.clonevol.colors(num.clusters)
    }
    names(colors) = sorted.cluster.names
    if (is.null(shapes)){
        shapes = seq(0, num.clusters)
    }
    names(shapes) = sorted.cluster.names
    if (!is.null(ignore.clusters)){
        ignore.clusters = as.character(ignore.clusters)
        var = var[!(var[[cluster.col.name]] %in% ignore.clusters),]
        colors = colors[!(names(colors) %in% ignore.clusters)]
        shapes = shapes[!(names(shapes) %in% ignore.clusters)]
    }
    clone.vafs = estimate.clone.vaf(var, cluster.col.name=cluster.col.name,
                       vaf.col.names=vaf.col.names,
                       vaf.in.percent=vaf.in.percent,
                       method=center.measure)
    #print(clone.vafs)
    if (vaf.in.percent){
        clone.vafs[,vaf.col.names] = clone.vafs[,vaf.col.names] * 100
    }
    if(!is.null(sample.names)){
        colnames(clone.vafs) = c(cluster.col.name, sample.names)
    }
    x = melt(clone.vafs, id.var=cluster.col.name)
    colnames(x) = c(cluster.col.name, 'sample', 'VAF')
    if (low.vaf.no.line){
        x = x[x$VAF >= min.cluster.vaf,]
    }
    x[[cluster.col.name]] = factor(x[[cluster.col.name]], levels=sorted.cluster.names)
    p = (ggplot(x, aes_string(x='sample', y='VAF'))
         + geom_point(aes_string(shape=cluster.col.name,
                                 color=cluster.col.name),
                      size=shape.size)
         + geom_line(aes_string(group=cluster.col.name,
                                color=cluster.col.name,
                                linetype=cluster.col.name),
                     size=line.size)
         + theme_bw(base_size=18)
         + scale_color_manual(values=colors)
         + scale_shape_manual(values=shapes)
         + scale_y_continuous(limits=c(0,max(x$VAF)*1.2))
         + theme(legend.key.width=unit(15,'mm'))
         + theme(panel.border=element_rect(linetype='solid',
                                           color='black',
                                           size=1))
         + theme(panel.grid.major=element_line(linetype='dotted',
                                               color='darkgray',
                                               size=0.5))
         + xlab(x.title)
         + ylab(y.title)
    )
    if (!is.null(out.file)){
        ggsave(p, file=out.file, width=width, height=height, useDingbats=FALSE)
    }else{
        return(p)
    }
}
#example
# plot.cluster.flow(var, vaf.col.names=vaf.col.names, out.file='tmp.pdf')

#boxplot.example()

#stop()



#' Plot values of columns pairwise
#' @description Plot value (eg. VAF) between samples pairwise-ly, annotated
#' and grouped by (eg.) clusters
#' @param data Variant data frames
#' @param col.names The columns to plot all pairwise, eg. VAFs of the samples
#' @param group.col.name The column used  to determine the category
#' when plotting to give different shapes,colors. eg. the cluster identity
#' of the variants
#' @param suffix String used to add to col.names as suffix in pairwise plot's
#' axis labels.
#' @param onePage Plot all pairwise plots in a page (default=TRUE)
#' @param multiPages Plot each pairwise plot in a separate page (default=FALSE)
#' of a multi-page PDF output file
#' @param xMin,xMax Min and max values for the x axis in multi-page plots (default=0-100)
#' @param yMin,yMax Min and max values for the y axis in multi-page plots (default=0-100)
#' @param xMinSmall,xMaxSmall Min and max values for the x axis in single-page plot (default=0,70)
#' @param yMinSmall,yMaxSmall Min and max values for the y axis in single-page plot (default=0,70)
#' @param out.prefix Output files' prefix
#' @import ggplot2
#' @export plot.pairwise
#' @examples
#' data(aml1)
#' plot.pairwise(aml1$variants, col.names = aml1$params$vaf.col.names,
#'               out.prefix = 'variants.pairwise.plot')
#'
plot.pairwise <- function(data,
                         col.names=c(),
                         suffix='',
                         group.col.name='cluster',
                         group.col.is.integer=TRUE,
                         colors=NULL,
                         shapes=NULL,
                         show.legend.title=FALSE,
                         sharedCategoryColName='',
                         onePage=TRUE, multiPages=FALSE,
                         xMin=0, xMax=100,
                         yMin=0, yMax=100,
                         xMinSmall=0, xMaxSmall=70,
                         yMinSmall=0, yMaxSmall=70,
                         show.none.zero.count=FALSE,
                         out.prefix=''){
    library(ggplot2)
    library(grid)
    library(gridExtra)
    n = length(col.names)
    nPlots = as.integer(n*(n-1)/2)
    smallPlots = list()
    if (group.col.is.integer){
        data[[group.col.name]] = as.integer(data[[group.col.name]])
        data[[group.col.name]] = factor(data[[group.col.name]],
            levels=unique(sort(data[[group.col.name]])))
    }else{
        data[[group.col.name]] = as.character(data[[group.col.name]])
    }
    num.groups = length(unique(data[[group.col.name]]))
    if (is.null(colors)){
        colors = get.ggplot2.colors(num.groups)
    }
    if(is.null(shapes)){
        shapes = seq(0, num.groups)
    }
    for (i in 1:(n-1)){
        for (j in (i+1):n){
            x = col.names[i]
            y = col.names[j]

            z = data[[x]] + data[[y]]
            nMutations = length(z[z > 0])

            # CI column names
            #xmax = paste(x, '_CI_hi', sep='')
            #xmin = paste(x, '_CI_lo', sep='')
            #ymax = paste(y, '_CI_hi', sep='')
            #ymin = paste(y, '_CI_lo', sep='')

            # big scatter plot
            # small scatter plot
            pSmall = (
                ggplot(data=data, aes_string(x=x, y=y))
                + geom_point(aes_string(shape=group.col.name,
                                        color=group.col.name))

                + theme_bw()
                + scale_shape_manual(values=shapes)
                + scale_color_manual(values=colors)
                + theme(panel.border=element_rect(linetype='solid',
                                                  size=1, color='black'))
                + scale_x_continuous(limits=c(xMinSmall,xMaxSmall))
                + scale_y_continuous(limits=c(yMinSmall,yMaxSmall))
                + theme(legend.position=c(0.925, 0.5))
                + theme(legend.key.height=unit(5,'mm'),
                        legend.key.width=unit(5,'mm'))
                + theme(legend.background = element_blank())
                + theme(panel.grid.major = element_line(linetype='dotted',
                                                        color='darkgray',
                                                        size=0.25))
                + theme(plot.margin=unit(c(1,1,1,1),"mm"))
                + xlab(paste0(x, suffix))
                + ylab(paste0(y, suffix))
            )
            if (show.none.zero.count){
                pSmall = pSmall + annotate("text", x=xMaxSmall*0.8,
                                           y=yMaxSmall*0.9,
                                           label=paste0('N=', nMutations))

            }
            if(!show.legend.title){
                pSmall = pSmall + theme(legend.title=element_blank())
            }
            #print(pSmall)
            smallPlots = c(smallPlots, list(pSmall))
        }
    }

    # Plot all scatter plots in one page
    nCols = ceiling(sqrt(nPlots))
    nRows = ceiling(nPlots/nCols)
    pdfOutFile = paste(out.prefix, '.scatter.1-page.pdf', sep='')
    pdf(file=pdfOutFile, width=3.5*nCols, height=3*nRows, useDingbats=FALSE)
    multiplot(plotlist=smallPlots, cols=nCols, horizontal=TRUE, e=0)
    dev.off()
    #system(paste('convert -density 200', pdfOutFile,
    #             paste(out.prefix,'.scatter.1-page.png', sep='')))
}


#' Merge variants and mapped events onto the same data frame that
#' is more convenient for plotting fancy boxplot of variants with
#' driver events highlighted
#'
#' @description Merge variants and mapped.events data frame to a single
#' data frame for fancy boxplot. Generally the columns of the two data frames
#' should be named the same, including (cluster.col.name, vaf.col.names,
#' cn.col.names, loh.col.names, column to highlight (driver), cluster, event).
#' event column should match between the two and events in mapped.events will
#' be excluded from variants data frame before merging them to prevent duplicate
#' @param variants: variant data frame
#' @param mapped.events: data frame of events mapped onto cluster/clone
#' that is output assign.events.to.clones()$events
#'
merge.variants.and.events <- function(variants, mapped.events=NULL,
                                       cluster.col.name='cluster',
                                       event.col.name='event',
                                       vaf.col.names,
                                       cn.col.names=c(),
                                       loh.col.names=c(),
                                       other.col.names=c()){

    cols = c(cluster.col.name, event.col.name, vaf.col.names,
        cn.col.names, loh.col.names, other.col.names)
    if (!all(cols %in% colnames(variants))){
        stop(paste0('ERROR: one of the following columns: ',
                    paste(cols, collapse=','),
                    ' not found in variants data frame\n'))

    }
    va = variants[, cols]
    va[[cluster.col.name]] = as.character(va[[cluster.col.name]])

    # add mapped.events to data.frame
    e = mapped.events
    if (!is.null(e)){
        if (!all(cols %in% colnames(variants))){
            stop(paste0('ERROR: one of the following columns: ',
                    paste(cols, collapse=','),
                    ' not found in mapped.events data frame\n'))
        }
        e = e[, cols]
        e[[cluster.col.name]] = as.character(e[[cluster.col.name]])
        va = va[!(va$event %in% e$event),] # remove duplicated event
        va = merge(e, va, all=TRUE)
    }
    # order by cluster # to make sure 1 is plotted before 2
    va = va[order(as.integer(as.character(va[[cluster.col.name]]))),]
    return(va)
}

# scale vaf to ccf
# method: mean or median
vaf2ccf <- function(df, founding.cluster, cluster.col.name='cluster',
                    vaf.col.names, method='mean'){
    if (method == 'mean'){
        founding.vaf = colMeans(df[df[[cluster.col.name]]==founding.cluster,
                                   vaf.col.names])
    }else if (method == 'median'){
        founding.vaf = apply(df[df[[cluster.col.name]]==founding.cluster,
                                vaf.col.names], 2, median)
    }
    ccf.scale = t(as.matrix(100/founding.vaf))
    ccf.scale = ccf.scale[rep(1,nrow(df)),]
    df[, vaf.col.names] = df[, vaf.col.names]*ccf.scale
    return(df)
}


