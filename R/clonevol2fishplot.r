#' Generate fishplot ready data from clonevol clonal evolution models
#' @description Generate fishplot ready data from clonevol clonal evolution models
#'
#' @param results: output from infer.clonal.models
#' @param rescale: rescale VAF such that no sum rule is violated (eg. parent
#' should have VAF >= sum of its chilren's VAF
#' @param samples: the names of the samples to be used in fishplot (this should
#' be a subset of vaf.col.names parameter provided to infer.clonal.models); The
#' order of samples provided will be prerseved in the fishplot
#' @export
#' @examples
#'   samples = c('P', 'R')
#'   v = data.frame(cluster=c(1,1,2,2),P=c(50,40,0,0),R=c(50,60,20,20))
#'   y = infer.clonal.models(variants=v, vaf.col.names=samples,
#'      founding.cluster=1)
#'   f = generateFishplotInputs(results=y)
#'\dontrun{
#'   fishes = createFishPlotObjects(f)
#'   pdf(tempfile(), width=8, height=5)
#'   for (i in 1:length(fishes)){
#'     fish = layoutClones(fishes[[i]])
#'     fish = setCol(fish,f$clonevol.clone.colors)
#'     fishPlot(fish,shape="spline", title.btm="Patient", cex.title=0.5,
#'            vlines=seq(1, length(samples)), vlab=samples, pad.left=0.5)
#'
#'   }
#'   dev.off()
#'}
#'
generateFishplotInputs <- function(results, rescale=T, samples=NULL){
  #no results, punt
  if (is.null(results$matched)){return(NULL)}

  #if specific samples aren't given, use them all
  if (is.null(samples)){
    samples = names(results$models)
  }

  #store the clonevol results in a list
  res = list(samples=samples, clonevol.clone.names=NULL, clonevol.clone.colors=NULL,
    timepoints=seq(1, length(samples)), num.models = nrow(results$matched$index),
    parents=list(), cell.fractions=list(), all=list())


  clonevol.clone.names = NULL
  clone.nums = NULL
  clonevol.clone.colors = NULL

  #create the needed inputs to fishplot
  for (i in 1:nrow(results$matched$index)){
    vv = NULL
    for (s in samples){
      v = results$models[[s]][[results$matched$index[i, s]]]
      if (rescale){v = rescale.vaf(v)}
      v = v[, c('lab', 'vaf', 'parent', 'color')]

      ## scale vaf and make cell.frac
      max.vaf = max(v$vaf)
      scale = 0.5/max.vaf*2*100
      v$vaf = v$vaf*scale
      v$vaf[v$vaf > 100] = 100# safeguard against rounding error making some vaf slightly > 100

      colnames(v) = c('clone', s , 'parent', 'color')
      v = v[!is.na(v$parent) & v$clone != '0',]
      if (is.null(vv)){vv = v}else{vv = merge(vv, v, all=T)}
    }
    for (s in samples){
      vv[is.na(vv[[s]]),s] = 0
    }
    vv = vv[order(as.integer(vv$clone)),]
    vv$parent[vv$parent == '-1'] = 0
    rownames(vv) = vv$clone

    ## fishplot requires clones to be named in sequential order. Do that, but
    ## store the clonevol-generated names and colors for pass-through
    if (is.null(clone.nums)){
      clone.nums = c(0, seq(1, nrow(vv)))
      names(clone.nums) = c(0, vv$clone)

      clonevol.clone.names = names(clone.nums)
      names(clonevol.clone.names) = as.character(clone.nums)
      res$clonevol.clone.names = clonevol.clone.names[-1]

      clonevol.clone.colors = c( 'white', vv$color)
      names(clonevol.clone.colors) = as.character(clone.nums)
      res$clonevol.clone.colors = clonevol.clone.colors[-1]

    }
    vv$clone = clone.nums[vv$clone]
    vv$parent = clone.nums[vv$parent]

    par = vv$parent
    frac = vv[, samples]
    res$parents[[i]] = par
    res$cell.fractions[[i]] = as.matrix(frac)
    res$all[[i]] = vv
  }
  return(res)
}

#' Create a list of fishplot objects that can then be called
#' @description Create a list of fishplot objects that can then be called
#' by layoutClones, then fishPlot
createFishPlotObjects <- function(results){
  #library(fishplot)

  fishes = list()
  for(i in 1:results$num.models){
    print(i)
    fishes[[i]]=createFishObject(results$cell.fractions[[i]],
      parents=as.integer(results$parents[[i]]),
      timepoints=results$timepoints,
      clone.labels=results$clonevol.clone.names[1:length(results$clonevol.clone.names)],
      fix.missing.clones=TRUE)
  }
  return(fishes)
}

