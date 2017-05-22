#' Plot all trees
plot.all.trees.clone.as.branch <- function(x, ...){
    for (i in 1:length(x$matched$merged.trees)){
        tg = x$matched$merged.trees[[i]]
        plot.tree.clone.as.branch(tg, ...)
    }
}


#' Save clonevol output to file
#' @param x: Output of infer.clonal.models
#' @param out.prefix: Output prefix
#'
save.clonevol.results <- function(x, out.prefix){
    save.table(x$variants, file=paste0(out.prefix, '.variants.tsv'))
    save.table(x$params, file)
    if (x$num.matched.models == 0){
        cat('No model to save.\n')
    }else{
        for (i in 1:x$num.matched.models){
            mt = x$matched$merged.trees[[i]]
            colnames(mt) = gsub('^lab$', 'clone', colnames(mt))
            save.table(mt, file=paste0(out.prefix, '.tree-', i, '.tsv'))
        }
        cat(x$num.matched.models , 'model(s) saved!')
    }
}


import.tree <- function(tree.file, variant.file){
    # read variants
    variants = read.table(variant.file, header=T, stringsAsFactors=T, sep='\t',
                          na.strings = c('', 'NA', '<NA>'), comment.char='')
    # read tree, and prepare clonevol merged tree data frame
    mt = read.table(tree.file, header=T, stringsAsFactors=F, sep='\t',
                    na.strings = c('', 'NA', '<NA>'), comment.char='')
    # check required columns
    required.cols = c('clone', 'parent', 'sample.with.nonzero.cell.frac.ci')
    if (!all(required.cols %in% colnames(mt))){
        cat('ERROR: The following required columns are missing:\n')
        cat(paste(required.cols[!(required.cols %in% colnames(mt))],
                  collapse='\n'))
    }
    # clonevol uses lab as col label for clone
    colnames(mt) = gsub('^clone$', 'lab', colnames(mt))
        # fill missing columns:
    if (!('color' %in% colnames(mt))){
        cat('Generating colors for clone using ClonEvol pallette...')
        mt$color = get.clonevol.colors(nrow(mt))
    }
    if (!('sample.group' %in% colnames(mt))){
        mt$sample.group = 'group1'
        mt$sample.group.color = 'black'
    }

    mt$excluded = FALSE
    y = list()
    y$num.matched.models = 1
    y$variants = variants
    y$matched = list()
    y$matched$merged.trees = list()
    y$matched$merged.trees[[1]] = mt
    y$matched$trimmed.merged.trees = list()
    y$matched$trimmed.merged.trees[[1]] = mt
    return(y)
}
