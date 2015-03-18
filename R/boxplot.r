
get.n <- function(x){
    return(c(y = mean(x), label = length(x)))
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
                      layout=NULL, horizontal=F, e=0.15) {
    require(grid)

    # Make a list from the ... arguments and plotlist
    plots = c(list(...), plotlist)

    numPlots = length(plots)

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
            x[i,]$hjust = hjust*2*(i%%2-0.5) + 0.5
            x[i,]$newX = x[i,]$newX - hjust*(i%%2-0.5)/5
        }
        df.hi[df.hi$cluster == c,] = x
    }
    #print(df.hi$hjust)
    return(df.hi)
}

# boxplot selected columns (names given by vaf.col.names) in data frame df
# group by cluster.col.name column
# eg. usage: boxPlot(t, 'cluster', vafColNames, 5, F, 'ppp.pdf')
# horizontal = T ==> all samples are lay out horizontally
# showClusterSize ==> show cluster size in the box
# Output: both pdf and png files
# width=0, height=0, w1=0, h1=0 (w/h = with/height of whole plot,
# w1/h1 = width/height of component plot)
# w1/h1 will orverwite w/h if they are non-zero. If all w/h/w1/h1 are
# zero, auto scale
# hightlight = index vector to select/subset df to hightlight using geom_point()
# eg highlight = df$tier == 'tier1' ==> hightlight tier1
# sizeName included to plot depth, etc.
# variant.class.col.name='tier' => summary based on this column
# if outPlotPrefix='', do not print output plot, return plot list
# hscale=1, vscale=1, ==> scale up width, height of the plot
variant.box.plot <- function(df,
                             cluster.col.name='cluster',
                             vaf.col.names=NULL,
                             vaf.limits=70,
                             variant.class.col.name='tier',

                             showClusterSize=F,
                             cluster.axis.name='cluster:',

                             sample.title.size=NULL,
                             panel.border.linetype='solid',
                             panel.border.linesize=1,
                             base_size=18, width=0, height=0,
                             width1=0, height1=0, hscale=1, vscale=1,
                             horizontal=F,

                             box=T,
                             box.line.type = 'solid',
                             box.line.size=0.5,
                             box.outlier.shape=1,
                             box.alpha=0.5,
                             violin=T,
                             violin.line.type = 'dotted',
                             violin.line.size=0.5,
                             jitter=F,
                             jitter.color='lightblue',
                             jitter.alpha=0.5,
                             jitter.size=1,
                             jitter.shape=3,


                             highlight=NULL,
                             highlight.color='red',
                             highlight.color.col.name=NULL,
                             highlight.size.names=NULL,
                             max.highlight.size.value=500,
                             highlight.size.legend.title='depth',
                             highlight.note.col.name = NULL,
                             highlight.note.color = 'blue',
                             highlight.note.size = 3,


                             ordered.x = NULL,
                             order.by.total.vaf=TRUE
){
    library(ggplot2)
    library(gridExtra)

    # make sure factor is converted to string first to avoid factor cluster
    # being treated as number
    df[[cluster.col.name]] = as.character(df[[cluster.col.name]])

    # order variants by decreasing total vafs
    if (!is.null(ordered.x)){
        cluster.orders = seq(1,length(ordered.x))
        names(cluster.orders) = ordered.x
        df = df[order(cluster.orders[df$cluster]),]
    }else if (order.by.total.vaf){
        df$total.vaf = apply(df[,vaf.col.names], 1, sum)
        mean.total.vafs = aggregate(total.vaf ~ cluster, df, mean)
        mean.total.vafs = mean.total.vafs[order(mean.total.vafs$total.vaf,
                                                decreasing=T),]
        rownames(mean.total.vafs) = mean.total.vafs$clusters
        mean.total.vafs.names = mean.total.vafs$cluster
        mean.total.vafs = mean.total.vafs$total.vaf
        names(mean.total.vafs) = mean.total.vafs.names
        df = df[order(mean.total.vafs[df$cluster], decreasing=T),]
    }

    # change cluster id to continous values to enable adjustment of postions
    # of highlighed genes
    cluster.levels = unique(df[[cluster.col.name]])
    df$cluster.label = as.character(df[[cluster.col.name]])
    df[[cluster.col.name]] = as.integer(factor(df[[cluster.col.name]],
                                               levels=cluster.levels))

    clusters = unique(df[[cluster.col.name]])
    x.axis.breaks = c(0, clusters)
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
        sizeName = NULL
        if (!is.null(highlight.size.names)){
            sizeName = highlight.size.names[ii]
        }

        # since violin plot will throw an error if there is zero variance
        # we'll add a very small number to first value of each cluster
        for (cl in clusters){
            df[df[[cluster.col.name]] == cl,][[yName]][1] =
                df[df[[cluster.col.name]] == cl,][[yName]][1] + 0.001
        }


        plotCnt = plotCnt + 1
        p = ggplot(data=df, aes_string(x = cluster.col.name, y = yName,
                                       group=cluster.col.name))
        if (jitter){
            p = p + geom_jitter(height = 0, color=jitter.color, size=jitter.size,
                                alpha=jitter.alpha, shape=jitter.shape)
        }

        if (box && violin){
            p = (p + geom_violin(scale='width', color=boxColor,
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
                                linetype=violin.line.type,
                                size=violin.line.size)
        }else{
            stop('Must specify at least boxplot or violin plot\n')
        }

        if (!is.null(highlight)){
            df.hi = df[highlight,]
            df.hi = randomizeHjust(df.hi, cluster.col.name=cluster.col.name,
                                   vaf.name=yName, hjust=0.75)
            if (!is.null(sizeName)){
                df.hi[[sizeName]] = cutBigValue(df.hi[[sizeName]],
                                                max.highlight.size.value)
            }
            #df.hi$note = paste0(df.hi$gene_name, '\n(',
            #    df.hi$amino_acid_change, ')')
            if (!is.null(highlight.color.col.name)){
                p = p + geom_point(data=df.hi,
                                   aes_string(x = 'newX', y=yName,
                                              size=sizeName,
                                              color=highlight.color.col.name),
                                   shape=1, show_guide=T)
            }else{
                p = p + geom_point(data=df.hi,
                                   aes_string(x = 'newX', y=yName,
                                              size=sizeName),
                                   color=highlight.color, shape=1, show_guide=T)

            }
            if (!is.null(sizeName)){
                size.breaks = c(0, 50, 100, 200, 300, 500)
                size.breaks = size.breaks[size.breaks <=
                                              max.highlight.size.value]
                size.labels = size.breaks
                size.labels[length(size.labels)] =
                    paste0('>', size.labels[length(size.labels)])
                p = (p + scale_size_continuous(name=highlight.size.legend.title,
                                           limits=c(0,max.highlight.size.value),
                                           breaks=size.breaks,
                                           labels=size.labels)
                     + theme(legend.position=c(0.7,0.9))
                )
            }
            p = p + geom_text(data=df.hi,
                              aes_string(x=cluster.col.name, y=yName,
                                         label=highlight.note.col.name,
                                         hjust='hjust'),
                              size=highlight.note.size,
                              color=highlight.note.color)

        }

        p = (
            p + scale_y_continuous(limits = c(0,vaf.limits[plotCnt]))
            + theme_bw(base_size=base_size)
            + theme(panel.border=element_rect(linetype=panel.border.linetype,
                                              size=panel.border.linesize,
                                              color='black'))
            + theme(plot.margin = unit(x = c(1, 1, 1, 1), units = "mm"))
        )
        if (showClusterSize){
            p = p + stat_summary(fun.data = get.n, geom = "text",
                                 position = position_dodge(height = 0,
                                                           width = 0.75),
                                 size = 5, color='blue')
        }
        if (horizontal){
            if (plotCnt > 1){
                p = p + theme(axis.title.x = element_blank())
            }else{
                p = p + scale_x_continuous(breaks = seq(1,nClusters),
                                           labels=clusterSizes)
            }
            p = p + coord_flip()
        }else{
            if (plotCnt < nPlots){
                p = (p + theme(axis.title.x = element_blank())
                     + scale_x_continuous(breaks = x.axis.breaks,
                                          labels=c(cluster.axis.name,cluster.labels))
                )
            }else{
                x.title = cluster.col.name
                if (!is.null(sumCnts)){
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

                    labs = c(cluster.axis.name, cluster.labels)
                }

                p = p + scale_x_continuous(x.title,
                                           breaks = x.axis.breaks,
                                           labels=labs)
            }
        }
        if (!is.null(sample.title.size)){
            p = p + theme(axis.title.y = element_text(size=sample.title.size))
        }
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
    print(e)
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

    if (horizontal){
        multiplot(plotlist=plots, cols=nPlots, horizontal=T)
    }else{
        multiplot(plotlist=plots, cols=1, horizontal=F, e=e)
    }

    return(plots)
}


# testing

boxplot.example <- function(){
    v = read.table('samples/CRC12.new.tsv', header=T, sep='\t', quote='',
                   stringsAsFactors=F)
    v = v[v$cluster != 'c9',]
    v = v[v$cluster != 'c11',]
    v = v[v$cluster != 'c23',]

    colnames(v) = gsub('CRC12_322_', '', colnames(v))
    colnames(v) = gsub('_\\d+.', '.', colnames(v))
    colnames(v) = gsub('C.VAF', 'Primary.VAF', colnames(v))
    colnames(v) = gsub('Li2.VAF', 'Li2.met.VAF', colnames(v))
    colnames(v) = gsub('Li3.VAF', 'Li3.met.VAF', colnames(v))
    colnames(v) = gsub('Li6.VAF', 'Li6.met.VAF', colnames(v))
    colnames(v) = gsub('C_XT1.VAF', 'Primary.xeno.VAF', colnames(v))
    colnames(v) = gsub('Li2_XT1.VAF', 'Li2.met.xeno.VAF', colnames(v))
    colnames(v) = gsub('Li3_XT1.VAF', 'Li3.met.xeno.VAF', colnames(v))
    hi = grepl('GJA8', v$gene_name)
    select=2:8
    vaf.col.names = grep('.VAF', colnames(v), fixed=T, value=T)[select]
    depth.col.names = grep('.depth', colnames(v), fixed=T, value=T)[select]
    pdf('test-out/CRC12.box.pdf', width=7, height=11, useDingbats=FALSE)
    variant.box.plot(v, vaf.col.names=vaf.col.names,
                     variant.class.col.name=NULL,
                     cluster.axis.name='',
                     sample.title.size=12,
                     highlight.size.names=depth.col.names,
                     max.highlight.size.value=200,
                     highlight.note.col.name='gene_name',
                     violin=F, box=T,
                     jitter=T, jitter.alpha=0.75, jitter.color='lightblue',
                     box.alpha=0.1, jitter.size=2,
                     jitter.shape=1,
                     highlight=NULL)
    dev.off()

    # CRC8
    v = read.table('samples/CRC8.tsv', header=T, sep='\t', quote='',
                   stringsAsFactors=F)
    v = v[v$cluster != 4,]
    v = v[v$cluster != 5,]
    v$cluster = paste0('c', v$cluster)

    colnames(v) = gsub('CRC8_237_', '', colnames(v))
    colnames(v) = gsub('_\\d+.', '.', colnames(v))
    colnames(v) = gsub('C.VAF', 'Primary.VAF', colnames(v))
    colnames(v) = gsub('Li2a.VAF', 'Li2a.met.VAF', colnames(v))
    colnames(v) = gsub('Li2b.VAF', 'Li2b.met.VAF', colnames(v))
    colnames(v) = gsub('Li8.VAF', 'Li8.met.VAF', colnames(v))
    select=2:5
    vaf.col.names = grep('.VAF', colnames(v), fixed=T, value=T)[select]
    depth.col.names = grep('.depth', colnames(v), fixed=T, value=T)[select]
    pdf('test-out/CRC8.box.pdf', width=3, height=7, useDingbats=FALSE)
    variant.box.plot(v, vaf.col.names=vaf.col.names,
                     variant.class.col.name=NULL,
                     cluster.axis.name='',
                     sample.title.size=12,
                     highlight.size.names=depth.col.names,
                     max.highlight.size.value=200,
                     highlight.note.col.name='gene_name',
                     violin=F, box=T,
                     jitter=T, jitter.alpha=0.75, jitter.color='lightblue',
                     box.alpha=0.1, jitter.size=2,
                     jitter.shape=1,
                     highlight=NULL)
    dev.off()
}



#boxplot.example()

#stop()





#### END HERE!!!!!!!!!!!!!!!!!!!
#### END HERE!!!!!!!!!!!!!!!!!!!
#### END HERE!!!!!!!!!!!!!!!!!!!
#### END HERE!!!!!!!!!!!!!!!!!!!
#### END HERE!!!!!!!!!!!!!!!!!!!
#### END HERE!!!!!!!!!!!!!!!!!!!
#### END HERE!!!!!!!!!!!!!!!!!!!















# plot values of columns pairwise
# example use cases: plot VAF pairwise
# Prepare the scatter plots
# categoryColNames: use this to color/shape points differently (eg. c('cluster', 'depth'))
# onePage=T => produce one page with all plots
# multiPages = T => produce many plots, each on one page
plotPairwise <- function(data=data.frame(), colNames=c(), categoryColNames=c(), sharedCategoryColName='', onePage=T, multiPages=T, xMin=0, xMax=100, yMin=0, yMax=100, xMinSmall=0, xMaxSmall=70, yMinSmall=0, yMaxSmall=70, outPrefix=''){
    n = length(colNames)
    nPlots = as.integer(n*(n-1)/2)
    smallPlots = list()
    bigPlots = list()
    nCategories = length(categoryColNames)
    if (nCategories != n){
        stop('Number of columns for categories does not match number of column for plot values\n')
    }

    for (i in 1:(n-1)){
        for (j in (i+1):n){
            x = colNames[i]
            y = colNames[j]

            z = data[[x]] + data[[y]]
            nMutations = length(z[z > 0])

            # CI column names
            #xmax = paste(x, '_CI_hi', sep='')
            #xmin = paste(x, '_CI_lo', sep='')
            #ymax = paste(y, '_CI_hi', sep='')
            #ymin = paste(y, '_CI_lo', sep='')
            # big scatter plot
            p = 0
            if (nCategories == 0){
                p = ggplot(data=data, aes_string(x=x, y=y)) + geom_point()
            }else if (nCategories == n){
                p = ggplot(data=data, aes_string(x=x, y=y, shape=categoryColNames[i], color=categoryColNames[j])) + geom_point() + scale_color_brewer()
            }
            p = (
                p
                #+ scale_shape_manual(values=seq(0,25))
                #+ geom_errorbarh(data=data, aes_string(y=y, xmin=xmin, xmax=xmax), size=0.1)
                #+ geom_errorbar(data=data, aes_string(x=x, ymin=ymin, ymax=ymax), size=0.1)
                + theme_bw()
                + theme(panel.border=element_rect(linetype='solid', size=1, color='black'))
                + scale_x_continuous(limits=c(xMin,xMax))
                + scale_y_continuous(limits=c(yMin,yMax))
                + annotate("text", x = xMax*0.8, y = yMax*0.95, label = paste('N=', nMutations, sep=''))
            )
            #print(p)
            bigPlots = c(bigPlots, list(p))
            # small scatter plot
            pSmall = (
                ggplot(data=data, aes_string(x=x, y=y))
                + geom_point()
                #+ scale_shape_manual(values=seq(0,25))
                + theme_bw()
                + theme(panel.border=element_rect(linetype='solid', size=1, color='black'))
                + scale_x_continuous(limits=c(xMinSmall,xMaxSmall))
                + scale_y_continuous(limits=c(yMinSmall,yMaxSmall))
                + annotate("text", x = xMaxSmall*0.8, y = yMaxSmall*0.9, label = paste('N=', nMutations, sep=''))
            )
            #print(pSmall)
            smallPlots = c(smallPlots, list(pSmall))
        }
    }

    # Plot all scatter plots in one page
    nCols = ceiling(sqrt(nPlots))
    nRows = ceiling(nPlots/nCols)
    pdfOutFile = paste(outPrefix, '.scatter.1-page.pdf', sep='')
    pdf(file=pdfOutFile, width=4*nCols, height=3*nRows)
    multiplot(plotlist=smallPlots, cols=nCols, horizontal=T, e=0)
    dev.off()
    system(paste('convert -density 200', pdfOutFile, paste(outPrefix, '.scatter.1-page.png', sep='')))

    # Plot scatter plots, each in one page, all together in a pdf file
    pdf(file=paste(outPrefix, '.scatter.multi-pages.pdf', sep=''), width=7, height=5)
    for (p in bigPlots){
        print(p)
    }
    dev.off()

}


# plot values of column pairs
# example use cases: plot WGS_VAF vs WES_VAF
# Prepare the scatter plots
# categoryColNames: use this to color/shape points differently (eg. c('cluster', 'depth'))
# onePage=T => produce one page with all plots
# multiPages = T => produce many plots, each on one page
plotPairs <- function(data=data.frame(), colNames1=c(), colNames2=c(), categoryColNames1=c(), categoryColNames2=c(), sharedCategoryColName='', onePage=T, multiPages=T, xMin=0, xMax=100, yMin=0, yMax=100, xMinSmall=0, xMaxSmall=70, yMinSmall=0, yMaxSmall=70, outPrefix=''){
    n = length(colNames1)
    n2 = length(colNames2)
    if (n != n2){stop('Number of columns not matched!\n')}
    nPlots = n
    smallPlots = list()
    bigPlots = list()
    nCategories = length(categoryColNames1)
    nCategories2 = length(categoryColNames2)
    if (nCategories != n || nCategories != nCategories2){
        stop('Number of columns for categories does not match, or does not match number of column for plot values\n')
    }

    for (i in 1:(n)){
        x = colNames1[i]
        y = colNames2[i]

        z = data[[x]] + data[[y]]
        nMutations = length(z[z > 0])

        # CI column names
        #xmax = paste(x, '_CI_hi', sep='')
        #xmin = paste(x, '_CI_lo', sep='')
        #ymax = paste(y, '_CI_hi', sep='')
        #ymin = paste(y, '_CI_lo', sep='')
        # big scatter plot
        p = 0
        if (nCategories == 0){
            p = ggplot(data=data, aes_string(x=x, y=y)) + geom_point()
        }else if (nCategories == n){
            p = ggplot(data=data, aes_string(x=x, y=y, shape=categoryColNames1[i], color=categoryColNames2[i])) + geom_point() + scale_color_brewer(type='seq', palette='BuPu')
        }
        p = (
            p
            #+ scale_shape_manual(values=seq(0,25))
            #+ geom_errorbarh(data=t, aes_string(y=y, xmin=xmin, xmax=xmax), size=0.1)
            #+ geom_errorbar(data=t, aes_string(x=x, ymin=ymin, ymax=ymax), size=0.1)
            + theme_bw()
            + theme(panel.border=element_rect(linetype='solid', size=1, color='black'))
            + scale_x_continuous(limits=c(xMin,xMax))
            + scale_y_continuous(limits=c(yMin,yMax))
            + annotate("text", x = xMax*0.8, y = yMax*0.95, label = paste('N=', nMutations, sep=''))
        )
        #print(p)
        bigPlots = c(bigPlots, list(p))
        # small scatter plot
        pSmall = (
            ggplot(data=t, aes_string(x=x, y=y))
            + geom_point()
            #+ scale_shape_manual(values=seq(0,25))
            + theme_bw()
            + theme(panel.border=element_rect(linetype='solid', size=1, color='black'))
            + scale_x_continuous(limits=c(xMinSmall,xMaxSmall))
            + scale_y_continuous(limits=c(yMinSmall,yMaxSmall))
            + annotate("text", x = xMaxSmall*0.8, y = yMaxSmall*0.9, label = paste('N=', nMutations, sep=''))
        )
        #print(pSmall)
        smallPlots = c(smallPlots, list(pSmall))
    }

    # Plot all scatter plots in one page
    nCols = ceiling(sqrt(nPlots))
    nRows = ceiling(nPlots/nCols)
    pdfOutFile = paste(outPrefix, '.scatter.1-page.pdf', sep='')
    pdf(file=pdfOutFile, width=4*nCols, height=3*nRows)
    multiplot(plotlist=smallPlots, cols=nCols, horizontal=T, e=0)
    dev.off()
    system(paste('convert -density 200', pdfOutFile, paste(outPrefix, '.scatter.1-page.png', sep='')))

    # Plot scatter plots, each in one page, all together in a pdf file
    pdf(file=paste(outPrefix, '.scatter.multi-pages.pdf', sep=''), width=7, height=5)
    for (p in bigPlots){
        print(p)
    }
    dev.off()

}




# plot values of column pairs
# example use cases: plot WGS_VAF vs WES_VAF
# Prepare the scatter plots
# categoryColNames: use this to color/shape points differently (eg. c('cluster', 'depth'))
# onePage=T => produce one page with all plots
# multiPages = T => produce many plots, each on one page
# plotPairs2: CI included
plotPairs2 <- function(data=data.frame(), colNames1=c(), colNames2=c(), categoryColNames1=c(), categoryColNames2=c(), sharedCategoryColName='', onePage=T, multiPages=T, xMin=0, xMax=100, yMin=0, yMax=100, xMinSmall=0, xMaxSmall=70, yMinSmall=0, yMaxSmall=70, plotCI=T, outPrefix=''){
    n = length(colNames1)
    n2 = length(colNames2)
    if (n != n2){stop('Number of columns not matched!\n')}
    nPlots = n
    smallPlots = list()
    bigPlots = list()
    nCategories = length(categoryColNames1)
    nCategories2 = length(categoryColNames2)
    if (nCategories != n || nCategories != nCategories2){
        stop('Number of columns for categories does not match, or does not match number of column for plot values\n')
    }

    for (i in 1:(n)){
        x = colNames1[i]
        y = colNames2[i]

        z = data[[x]] + data[[y]]
        nMutations = length(z[z > 0])

        # CI column names
        xmax = paste(x, '_CI_hi', sep='')
        xmin = paste(x, '_CI_lo', sep='')
        ymax = paste(y, '_CI_hi', sep='')
        ymin = paste(y, '_CI_lo', sep='')
        # big scatter plot
        p = 0
        if (nCategories == 0){
            p = ggplot(data=data, aes_string(x=x, y=y)) + geom_point()
        }else if (nCategories == n){
            p = ggplot(data=data, aes_string(x=x, y=y, shape=categoryColNames1[i], color=categoryColNames2[i])) + geom_point() + scale_color_brewer(type='seq', palette='BuPu')
        }
        p = (
            p
            #+ scale_shape_manual(values=seq(0,25))
            + theme_bw()
            + theme(panel.border=element_rect(linetype='solid', size=1, color='black'))
            + scale_x_continuous(limits=c(xMin,xMax))
            + scale_y_continuous(limits=c(yMin,yMax))
            + annotate("text", x = xMax*0.8, y = yMax*0.95, label = paste('N=', nMutations, sep=''))
        )
        if (plotCI){
            p = (
                p + geom_errorbarh(data=t, aes_string(y=y, xmin=xmin, xmax=xmax), size=0.1)
                + geom_errorbar(data=t, aes_string(x=x, ymin=ymin, ymax=ymax), size=0.1)
            )

        }
        #print(p)
        bigPlots = c(bigPlots, list(p))
        # small scatter plot
        pSmall = (
            ggplot(data=t, aes_string(x=x, y=y))
            + geom_point()
            #+ scale_shape_manual(values=seq(0,25))
            + theme_bw()
            + theme(panel.border=element_rect(linetype='solid', size=1, color='black'))
            + scale_x_continuous(limits=c(xMinSmall,xMaxSmall))
            + scale_y_continuous(limits=c(yMinSmall,yMaxSmall))
            + annotate("text", x = xMaxSmall*0.8, y = yMaxSmall*0.9, label = paste('N=', nMutations, sep=''))
        )
        #print(pSmall)
        smallPlots = c(smallPlots, list(pSmall))
    }

    # Plot all scatter plots in one page
    nCols = ceiling(sqrt(nPlots))
    nRows = ceiling(nPlots/nCols)

    # only plot all in 1 page if no error bar is plot
    if (!plotCI){
        pdfOutFile = paste(outPrefix, '.scatter.1-page.pdf', sep='')
        pdf(file=pdfOutFile, width=4*nCols, height=3*nRows)
        multiplot(plotlist=smallPlots, cols=nCols, horizontal=T, e=0)
        dev.off()
        system(paste('convert -density 200', pdfOutFile, paste(outPrefix, '.scatter.1-page.png', sep='')))
    }

    # Plot scatter plots, each in one page, all together in a pdf file
    pdf(file=paste(outPrefix, '.scatter.multi-pages.pdf', sep=''), width=7, height=5)
    for (p in bigPlots){
        print(p)
    }
    dev.off()

}

