# Wrappers for ClonEvol visualization

#' Plot all consensus trees along with the bell plots and sphere of cells
#' showing clonal evolution and admixture of clones for individual samples
#' @description  Plot all consensus trees along with the bell plots and sphere
#' of cells for individual samples
#' @param x Output of infer.clonal.models
#' @param height Height of plot in inch (default = auto)
#' @param width Width of plot in inch (default = 4)
#' @param models Index of models to plot
#' @param trees.per.page Number of models (trees) per page
#' @param tree.branch.width Branch width of trees
#' @param tree.branch.angle Angle of branch of trees
#' @param tree.label Root of tree label (default = "Normal")
#' @panel.widths vector of spaces (widths, size=4) of the panels for tree,
#' bell plot, sphere of cells, and sample title. Default=NULL (auto)
#' tree.branch.rescale Rescale branch length (default = NULL, can be
#' "linear" (equivalent to "none"), "sqrt", or "log2")
#' @param out.pdf.file Output file (ending with .pdf or .png)
#' 
plotTreesBellsCells <- function(x, width=4, height=NULL, panel.widths=NULL,
    models=NULL,
    trees.per.page=1, tree.branch.width=1, tree.branch.angle=40,
    tree.label='Normal', tree.branch.rescale=NULL,
    show.event=TRUE, show.pruned.tree.id=TRUE,
    cell.border.size=0.1, show.frame=TRUE,
    hili.mid=NULL, p2m.clones=NULL, m2m.clones=NULL, xeno.clones=NULL,
    clone.time.step.scale=0.85,
    bell.curve.step=1,
    include.cluster=FALSE, out.pdf.file=NULL){

    w = width; h = height; ncols = trees.per.page
    if (!is.null(tree.branch.rescale)){
        if (tree.branch.rescale == 'linear'){tree.branch.rescale='none'}
        x = convert.consensus.tree.clone.to.branch(x,
            branch.scale=tree.branch.rescale)
    }

    # get sample names and number
    samples = x$params$vaf.col.names
    samples = c(samples[samples != 'P1'], samples[samples == 'P1'])
    samples = rev(samples)
    num.samples = length(samples)
    num.models = x$num.matched.models
    n = length(samples)
    
    # select models to plot
    if (is.null(models)){
        models = 1:num.models
    }else{
        notfound = setdiff(models, 1:num.models)
        if (length(notfound) > 0){
            stop(paste0('ERROR: Models not found: ',
                paste(notfound, collapse=','), '\n'))
        }
    }

    # group models with same pruned trees
    map = x$matched$trimmed.merged.trees.map
    tmp = map[map$tree.idx %in% models,]
    tmp = tmp[order(tmp$pruned.tree.idx),]
    models = tmp$tree.idx
 
    # prepare with/height of plot
    if (is.null(h)){h = n*ncols*0.4}
    if (!is.null(out.pdf.file)){
        pdf(out.pdf.file, width=w, height=h, useDingbats=F, title='')
    }

    # layout the plots on the same page
    mata = c() # mat all
    wwa  = c()
    hha  = c()
    par(mfrow=c(n*ncols,4+include.cluster), mar=c(0,0,0,0))
    for (col in 1:ncols){

        # layout 1st col is tree, 2nd col is multiple bell, each on one row,
        # 3rd and 4rd is something else
        mat = cbind(rep(1, n), seq(2,n+1), seq(n+2, 2*n+1), seq(2*n+2, 3*n+1))
        if (!is.null(mata)){mat = mat + max(mata)}
        mata = rbind(mata, mat)
    }
    hh = rep(1, n*ncols)
    ww = c(7,1.5,1,0.3)
    if (!is.null(panel.widths)){ ww = panel.widths}
    if (include.cluster){
        mata = cbind(rep(1, nrow(mata)),mata+1)
        ww = c(4, ww)
    }
    # reverse layout matrix, so plot is printed from bottom to top
    # and when rotating 90 degree, it becomes left to right
    mata = mata[nrow(mata):1,]
    layout(mata, ww, hh)

    # plot tree/bells/cells for all/selected models
    pruned.tree.id = 0
    box.cols = rep(c('gray', 'lightblue'), 100)
    for (k in models){
        mt = x$matched$merged.trees[[k]]
        prev.pruned.tree.id = pruned.tree.id
        pruned.tree.id = map$pruned.tree.idx[map$tree.idx==k]
        this.tree.label = tree.label
        if (show.pruned.tree.id){
            this.tree.label = paste0('Pruned tree ', pruned.tree.id,
                ' : tree ' , k, ' / ', tree.label)
        }
        cat('Plotting model ', k, '\n')
        #print(mt[, c('lab', 'parent')])
        # manipulate merged tree
        xeno = p2m = p2m.b = m2m = m2m.b = FALSE
        if (is.null(p2m.clones) & is.null(m2m.clones) & is.null(xeno.clones)){
            scol = 'sample.with.nonzero.cell.frac.ci'
            met.founder = grepl('*L', mt[[scol]], fixed=T)
            prim = grepl('P|*P', mt[[scol]], fixed=T)
            p2m = met.founder
        }else{
            p2m = mt$lab %in% p2m.clones
            p2m.b = mt$lab %in% p2m.branches
            m2m = mt$lab %in% m2m.clones
            m2m.b = mt$lab %in% m2m.branches
            xeno = mt$lab %in% xeno.clones
        }
        #mt$branch.border.color=mt$color
        mt$node.border.color = '#525252'
        mt$node.border.width = '0.5'

        if (any(xeno)){
            mt$node.border.color[xeno] = 'green'
            mt$node.border.width[xeno] = 1
        }

        if (any(p2m)){
            mt$node.border.color[p2m] = 'blue'
            mt$node.border.width[p2m] = 1
        }
        if (any(p2m.b)){
            mt$branch.border.color[p2m.b] = 'blue'
            mt$branch.border.width[p2m.b] = 1
            mt$branch.border.linetype[p2m.b] = 'solid'
        }
        if (any(m2m)){
            mt$node.border.color[m2m] = 'red'
            mt$node.border.width[m2m] = 1
        }
        if (any(m2m.b)){
            mt$branch.border.color[m2m] = 'red'
            mt$branch.border.width[m2m] = 2
            mt$branch.border.linetype[m2m] = 'solid'
        }

        # plot tree on 1st col
        mt = normalizeBranchLength(mt, 20)  # normalize length
        plot.tree.clone.as.branch(mt,
                    tree.rotation=90, text.angle=90, angle=tree.branch.angle,
                    branch.width=tree.branch.width, branch.text.size=0.25,
                    node.size=2, node.label.size=0.75, node.text.size=0.5, event.sep.char=',',
                    show.event=show.event, tree.label = this.tree.label)

        # alternate color of surrounding box for different pruned trees
        if (pruned.tree.id != prev.pruned.tree.id){box.cols = box.cols[-1]}
        hili.box.col=box.cols[1]
        if (!is.null(hili.mid) && k == hili.mid){
            hili.box.col='red'
        }
        if (show.frame){ box("figure", lwd=0.5, col=hili.box.col) }

        # plot all bells and cell pops in 2nd col    
        for (s in samples){
            # get model data
            s.match.idx = x$matched$index[[s]][k]
            m = x$models[[s]][[s.match.idx]]
            #cat(s, i, '\n')
            draw.sample.clones(m, x=2, y=0, wid=40, len=9, clone.shape='bell',
                               bell.curve.step=bell.curve.step,
                               clone.time.step.scale=clone.time.step.scale,
                               bell.border.width=0.1, show.clone.label=F,
                               #zero.cell.frac.clone.color='white',
                               zero.cell.frac.clone.border.color='fill',
                               #nonzero.cell.frac.clone.border.color='black',
                               nonzero.cell.frac.clone.border.width=0.1,
                               zero.cell.frac.clone.border.width=0.1,
                               color.border.by.sample.group=F,
                               text.size=0.75,disable.cell.frac=T, show.time.axis=F
            )
            if (show.frame){ box("figure", lwd=0.5, col=hili.box.col) }
        }

        # plot cell population
        for (s in samples){
            # get model data
            s.match.idx = x$matched$index[[s]][k]
            m = x$models[[s]][[s.match.idx]]

            m = m[!m$excluded & !m$is.zero,]
            pa = plot.cell.population(m$free.mean/sum(m$free.mean), m$color, layout='cloud',
                cell.border.size=cell.border.size, cell.border.color='black',
                clone.grouping='horizontal', frame=F,
                num.cells=100)

            plot.new()
            vps = baseViewports()
            pushViewport(vps$figure)
            vp = plotViewport(c(0,0,0,0))
            print(pa, vp=vp)
            if (show.frame) { box("figure", lwd=0.5, col=hili.box.col) }
            popViewport()
        }

        for (s in samples){
            plot(1, type="n", xlab="", ylab="", xlim=c(0, 2), ylim=c(0, 2), axes=F, ann=F);
            text(1,1, s, cex=1, srt=90)
            if (show.frame){ box("figure", lwd=0.5, col=hili.box.col) }
        }

    }

    if (!is.null(out.pdf.file)){ dev.off() }
    cat('\nPlots written to file. Turn the pages 90 degree clockwise for best view.\n')
}





#' Plot all bells and cellular population of consensus trees 
#' showing clonal evolution and admixture of clones for individual samples
#' @description  Plot all consensus trees along with the bell plots and sphere
#' of cells for individual samples
#' @param x Output of infer.clonal.models
#' @param height Height of plot in inch (default = auto)
#' @param width Width of plot in inch (default = auto)
#' @param models Index of models to plot
#' @param trees.per.page Number of models (trees) per page
#' @param tree.branch.width Branch width of trees
#' @param tree.branch.angle Angle of branch of trees
#' @param tree.label Root of tree label (default = "Normal")
#' tree.branch.rescale Rescale branch length (default = NULL, can be
#' "linear" (equivalent to "none"), "sqrt", or "log2")
#' @param out.pdf.file Output file (ending with .pdf or .png)
#' 
plotBellsCells <- function(x, width=NULL, height=NULL, models=NULL,
    trees.per.page=1, cell.border.size=0.1, show.frame=TRUE,
    hili.mid=NULL, p2m.clones=NULL, m2m.clones=NULL, xeno.clones=NULL,
    clone.time.step.scale=0.85, bell.curve.step=1,
    out.pdf.file=NULL){

    w = width; h = height; ncols = trees.per.page

    # get sample names and number
    samples = x$params$vaf.col.names
    samples = c(samples[samples != 'P1'], samples[samples == 'P1'])
    samples = rev(samples)
    num.samples = length(samples)
    num.models = x$num.matched.models
    n = length(samples)
    
    # select models to plot
    if (is.null(models)){
        models = 1:num.models
    }else{
        notfound = setdiff(models, 1:num.models)
        if (length(notfound) > 0){
            stop(paste0('ERROR: Models not found: ',
                paste(notfound, collapse=','), '\n'))
        }
    }

    # group models with same pruned trees
    map = x$matched$trimmed.merged.trees.map
    tmp = map[map$tree.idx %in% models,]
    tmp = tmp[order(tmp$pruned.tree.idx),]
    models = tmp$tree.idx
 
    # prepare with/height of plot
    if (is.null(h)){h = n*0.4}
    if (is.null(w)){w = 1*ncols}
    if (!is.null(out.pdf.file)){
        pdf(out.pdf.file, width=w, height=h, useDingbats=F, title='')
    }

    # layout the plots on the same page
    mata = c() # mat all
    wwa  = c()
    hha  = c()
    par(mfrow=c(n*ncols,3), mar=c(0,0,0,0))
    for (col in 1:ncols){
        # layout 1st col is bells,  2nd col is cells, 3rd col is sample name
        mat = matrix(1:(3*n), nrow=n)
        if (!is.null(mata)){mat = mat + max(mata)}
        mata = cbind(mata, mat)
    }
    hh = rep(1, n*ncols)
    ww = c(1.5,1,0.3)
    # reverse layout matrix, so plot is printed from bottom to top
    # and when rotating 90 degree, it becomes left to right
    mata = mata[nrow(mata):1,]
    layout(mata, ww, hh)

    # plot tree/bells/cells for all/selected models
    pruned.tree.id = 0
    box.cols = rep(c('gray', 'lightblue'), 100)
    for (k in models){
        mt = x$matched$merged.trees[[k]]
        prev.pruned.tree.id = pruned.tree.id
        pruned.tree.id = map$pruned.tree.idx[map$tree.idx==k]
        #this.tree.label = tree.label
        #if (show.pruned.tree.id){
        #    this.tree.label = paste0('Pruned tree ', pruned.tree.id,
        #        ' : tree ' , k, ' / ', tree.label)
        #}
        cat('Plotting model ', k, '\n')
        #print(mt[, c('lab', 'parent')])
        # manipulate merged tree
        xeno = p2m = p2m.b = m2m = m2m.b = FALSE
        if (is.null(p2m.clones) & is.null(m2m.clones) & is.null(xeno.clones)){
            scol = 'sample.with.nonzero.cell.frac.ci'
            met.founder = grepl('*L', mt[[scol]], fixed=T)
            prim = grepl('P|*P', mt[[scol]], fixed=T)
            p2m = met.founder
        }else{
            p2m = mt$lab %in% p2m.clones
            p2m.b = mt$lab %in% p2m.branches
            m2m = mt$lab %in% m2m.clones
            m2m.b = mt$lab %in% m2m.branches
            xeno = mt$lab %in% xeno.clones
        }
        #mt$branch.border.color=mt$color
        mt$node.border.color = '#525252'
        mt$node.border.width = '0.5'

        if (any(xeno)){
            mt$node.border.color[xeno] = 'green'
            mt$node.border.width[xeno] = 1
        }

        if (any(p2m)){
            mt$node.border.color[p2m] = 'blue'
            mt$node.border.width[p2m] = 1
        }
        if (any(p2m.b)){
            mt$branch.border.color[p2m.b] = 'blue'
            mt$branch.border.width[p2m.b] = 1
            mt$branch.border.linetype[p2m.b] = 'solid'
        }
        if (any(m2m)){
            mt$node.border.color[m2m] = 'red'
            mt$node.border.width[m2m] = 1
        }
        if (any(m2m.b)){
            mt$branch.border.color[m2m] = 'red'
            mt$branch.border.width[m2m] = 2
            mt$branch.border.linetype[m2m] = 'solid'
        }

        # alternate color of surrounding box for different pruned trees
        if (pruned.tree.id != prev.pruned.tree.id){box.cols = box.cols[-1]}
        hili.box.col=box.cols[1]
        if (!is.null(hili.mid) && k == hili.mid){
            hili.box.col='red'
        }
        if (show.frame){ box("figure", lwd=0.5, col=hili.box.col) }

        # plot all bells and cell pops in 2nd col    
        for (s in samples){
            # get model data
            s.match.idx = x$matched$index[[s]][k]
            m = x$models[[s]][[s.match.idx]]
            #cat(s, i, '\n')
            draw.sample.clones(m, x=2, y=0, wid=40, len=9, clone.shape='bell',
                               bell.curve.step=bell.curve.step,
                               clone.time.step.scale=clone.time.step.scale,
                               bell.border.width=0.1, show.clone.label=F,
                               #zero.cell.frac.clone.color='white',
                               zero.cell.frac.clone.border.color='fill',
                               #nonzero.cell.frac.clone.border.color='black',
                               nonzero.cell.frac.clone.border.width=0.1,
                               zero.cell.frac.clone.border.width=0.1,
                               color.border.by.sample.group=F,
                               text.size=0.75,disable.cell.frac=T, show.time.axis=F
            )
            if (show.frame){ box("figure", lwd=0.5, col=hili.box.col) }
        }

        # plot cell population
        for (s in samples){
            # get model data
            s.match.idx = x$matched$index[[s]][k]
            m = x$models[[s]][[s.match.idx]]

            m = m[!m$excluded & !m$is.zero,]
            pa = plot.cell.population(m$free.mean/sum(m$free.mean), m$color, layout='cloud',
                cell.border.size=cell.border.size, cell.border.color='black',
                clone.grouping='horizontal', frame=F,
                num.cells=100)

            plot.new()
            vps = baseViewports()
            pushViewport(vps$figure)
            vp = plotViewport(c(0,0,0,0))
            print(pa, vp=vp)
            if (show.frame) { box("figure", lwd=0.5, col=hili.box.col) }
            popViewport()
        }

        for (s in samples){
            plot(1, type="n", xlab="", ylab="", xlim=c(0, 2), ylim=c(0, 2), axes=F, ann=F);
            text(1,1, s, cex=1, srt=90)
            if (show.frame){ box("figure", lwd=0.5, col=hili.box.col) }
        }

    }

    if (!is.null(out.pdf.file)){ dev.off() }
    cat('\nPlots written to file. Turn the pages 90 degree clockwise for best view.\n')
}



