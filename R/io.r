#' Plot all branch-based consensus clonal evolution trees
#' @description Plot all consensus trees using the branch-based visualizations
#' @param x Output of infer.clonal.models function
#' @param ... Other parammeters of plot.tree.clone.as.branch function
#' @seealso plot.tree.clone.as.branch
#' @export plot.all.trees.clone.as.branch
#' @examples
#' data(aml1)
#' y = aml1
#' \dontrun{
#' plot.all.trees.clone.as.branch(y, branch.width = 0.5,
#'      node.size = 1, node.label.size = 0.5)
#' }
#'
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
    write.table(x$variants, file=paste0(out.prefix, '.variants.tsv'))
    write.table(x$params, file)
    if (x$num.matched.models == 0){
        cat('No model to save.\n')
    }else{
        for (i in 1:x$num.matched.models){
            mt = x$matched$merged.trees[[i]]
            colnames(mt) = gsub('^lab$', 'clone', colnames(mt))
            write.table(mt, file=paste0(out.prefix, '.tree-', i, '.tsv'))
        }
        cat(x$num.matched.models , 'model(s) saved!')
    }
}


#' Import trees to clonevol
#' @description Import trees (eg. those predicted by other methods) to clonevol
#' thus we can visualize them using clonevol visualization features
#' @param tree.file A TSV file that contains the tree with minimum 3 labeled columns:
#' clone, parent, sample.with.nonzero.cell.frac.ci. It can also have additional
#' useful columns such as color, and events
#' @param variant.file A TSV file containing the variants, with cluster column,
#' vaf.col.names, etc.
#' @return The clonal evolution tree similar to the output of infer.clonal.models
#' function that can be used in multiple plotting functions as if it is the output
#' of infer.clonal.models.
#' @export import.tree
#' @examples
#' \dontrun{
#' y = import.tree('trees.tsv', 'variants.tsv')
#' }
#'
import.tree <- function(tree.file, variant.file){
    # read variants
    variants = read.table(variant.file, header=TRUE, stringsAsFactors=TRUE, sep='\t',
                          na.strings = c('', 'NA', '<NA>'), comment.char='')
    # read tree, and prepare clonevol merged tree data frame
    mt = read.table(tree.file, header=TRUE, stringsAsFactors=FALSE, sep='\t',
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
