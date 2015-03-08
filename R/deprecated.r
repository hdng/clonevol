#' (Deprecated)
#' This function does not work well for >=3 samples
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
#' @param model: cancer evolution model to use, c('monoclonal', 'polyclonal').
#' monoclonal model assumes the orginal tumor (eg. primary tumor) arises from
#' a single normal cell; polyclonal model assumes the original tumor can arise
#' from multiple cells (ie. multiple founding clones). In the polyclonal model,
#' the total VAF of the separate founding clones must not exceed 0.5
#'
#'
infer.clonal.models.2 <- function(c, model='monoclonal'){
  # check format of input, find samples
  cnames = colnames(c)
  if (!('cluster' %in% cnames && length(cnames) >= 2)){
    stop('ERROR: Wrong dat format. No cluster column and/or no sample\n')
  }
  samples = setdiff(cnames, 'cluster')
  nSamples = length(samples)
  cat('Primary tumor sample: ', samples[1], '\n')
  if (nSamples >= 2){
    for (i in 2:nSamples){
      cat('Met/relapse sample ', i-1, ': ', samples[i], '\n', sep='')
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
  cat('Using ', model, ' model\n', sep='')

  # prepare cluster data and infer clonal models for individual samples
  vv = list()
  for (s in samples){
    cat(s, ':')
    v = make.clonal.data.frame(c[[s]], c$cluster, add.normal)
    models = enumerate.clones(v)
    cat(length(models), 'clonal architecture model(s) found\n')
    if (length(models) == 0){
      stop(paste('ERROR: No clonal models for sample:', s,
                 '\nCheck data or remove this sample, then re-run.\n'))
    }else{
      vv[[s]] = models
    }
  }

  # infer clonal evolution models, given all evolve from the 1st sample
  matched = NULL
  if (nSamples == 1 && length(vv[[1]]) > 0){
    matched = data.frame(x=1:length(vv[[1]]))
    colnames(matched) = c(samples[1])
  }
  if (nSamples >= 2){
    prims = vv[[1]]
    for (prim.model in 1:length(prims)){
      matched.models = rep(0, nSamples)
      names(matched.models) = samples
      matched.models[1] = prim.model
      prim = prims[[prim.model]]
      for (met.idx in 2:length(vv)){
        mets = vv[[met.idx]]
        for (met.model in 1:length(mets)){
          matched.models[met.idx] = 0
          met = mets[[met.model]]
          if (match.sample.clones(prim, met)){
            matched.models[met.idx] = met.model
            # if all met samples are compatible with primary
            if (all(matched.models > 0)){
              if (is.null(matched)){
                matched = as.data.frame(as.list(matched.models))
              }else{
                matched = rbind(matched,
                                as.data.frame(as.list(matched.models)))
              }
            }
          }
        }
      }
    }
    if (!is.null(matched)){
      rownames(matched) = seq(1,nrow(matched))
    }
  }
  cat(paste0('Found ', ifelse(is.null(matched), 0, nrow(matched)),
             ' compatible evolutional models\n'))
  return (list(models=vv, matched=matched))

}



#' (Deprecated)
#' Find matched models between samples
#' infer clonal evolution models, given all evolve from the 1st sample
#' TODO: This function does not care about met met matching. Obsolete!
find.matched.models.1 <- function(vv, samples){
  nSamples = length(samples)
  matched = NULL
  find.next.match <- function(prim.idx, prim.model.idx,
                              met.idx, met.model.idx, matched.models){
    if (met.idx > nSamples){
      if (all(matched.models > 0)){
        matched <<- rbind(matched, matched.models)
      }
    }else{
      if (match.sample.clones(vv[[prim.idx]][[prim.model.idx]],
                              vv[[met.idx]][[met.model.idx]])){
        matched.models[[met.idx]] = met.model.idx
        find.next.match(prim.idx, prim.model.idx, met.idx+1, 1, matched.models)
      }else{
        if (length(vv[[met.idx]]) > met.model.idx){
          find.next.match(prim.idx, prim.model.idx,
                          met.idx, met.model.idx+1, matched.models)
        }
        #find.next.match(prim.idx)
      }
    }
  }
  for (prim.model in 1:length(vv[[1]])){
    # cat('prim.model =', prim.model, '\n')
    find.next.match(1, prim.model, 2, 1, c(prim.model,rep(0,nSamples-1)))
  }
  return(matched)
}


#' (Deprecated)
#' Enumerate all possible clonal structures for a single sample using
#' mono-clonal model (ie. the tumor originated from a single cell)
#'
#' @param v: a data frame that contains clonal structure pre-identified from
#' clonality analysis (most often by clustering variants using their variant
#' allele frequencies). This data frame has two columns:
#'    cluster: a number indicating the cluster of variants
#'    vaf: variant allele frequency estimated for the cluster
#'
#' Example data frame:\cr
#' cluster vaf\cr
#' 1       0.5\cr
#' 2       0.3\cr
#' 3       0.1\cr
enumerate.clones.monoclonal <- function(v){
  vv = list()
  findParent <- function(v, i){
    if (i > nrow(v)){
      #debug
      #print(v)
      vv <<- c(vv, list(v))
    }else{
      vaf = v[i,]$vaf
      if (vaf == 0){
        vx = v
        vx$parent[i] = NA
        findParent(vx, i+1)
      }else{
        #for (j in 1:(i-1)){
        for (j in 1:nrow(v)){
          if (i != j && v$free[j] >= vaf && v$parent[j] != v$lab[i]){
            vx = v
            vx$free[j] = vx$free[j] - vaf
            vx$occupied[j] = vx$occupied[j] + vaf
            vx$parent[i] = vx[j,]$lab
            # debug
            #cat(i, '<-', j, 'vaf=', vaf, '\n')
            findParent(vx, i+1)
          }
        }
      }
    }
  }
  findParent(v,2)
  return(vv)
}

# (Incompleted code)
# Test if total means of x clusters greater than mean of y cluster
# c has two columns (1st is cluster id, 2nd is vaf)
fittable <- function(c, x, y){
  variants = crc12.variants
  variants[variants$cluster %in% c(5,6,7),]$cluster = 5
  variants[variants$cluster > 5,]$cluster = variants[variants$cluster > 5,]$cluster - 2
  vaf.col.names = grep('.WGS_VAF', colnames(variants), value=T, fixed=T)
  vaf.col.names = vaf.col.names[!grepl('normal', vaf.col.names)]

  cc = variants[, c('cluster', vaf.col.names)]
  c = cc[,c(1,)]

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
#' @param model: cancer evolution model to use, c('monoclonal', 'polyclonal').
#' monoclonal model assumes the orginal tumor (eg. primary tumor) arises from
#' a single normal cell; polyclonal model assumes the original tumor can arise
#' from multiple cells (ie. multiple founding clones). In the polyclonal model,
#' the total VAF of the separate founding clones must not exceed 0.5
#'
#'
infer.clonal.models.absolute <- function(c, model='monoclonal'){
  # check format of input, find samples
  cnames = colnames(c)
  if (!('cluster' %in% cnames && length(cnames) >= 2)){
    stop('ERROR: Wrong dat format. No cluster column and/or no sample\n')
  }
  samples = setdiff(cnames, 'cluster')
  nSamples = length(samples)
  cat('Primary tumor sample: ', samples[1], '\n')
  if (nSamples >= 2){
    for (i in 2:nSamples){
      cat('Met/relapse sample ', i-1, ': ', samples[i], '\n', sep='')
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
  cat('Using ', model, ' model\n', sep='')

  # prepare cluster data and infer clonal models for individual samples
  vv = list()
  for (s in samples){
    cat(s, ': ')
    v = make.clonal.data.frame(c[[s]], c$cluster, add.normal)

    models = enumerate.clones.absolute(v)
    cat(length(models), 'clonal architecture model(s) found\n')
    if (length(models) == 0){
      print(v)
      stop(paste('ERROR: No clonal models for sample:', s,
                 '\nCheck data or remove this sample, then re-run.\n'))
    }else{
      vv[[s]] = models
    }
  }

  # infer clonal evolution models, given all evolve from the 1st sample
  matched = NULL
  scores = NULL
  if (nSamples == 1 && length(vv[[1]]) > 0){
    num.models = length(vv[[1]])
    matched = data.frame(x=1:num.models)
    colnames(matched) = c(samples[1])
  }
  if (nSamples >= 2){
    z = find.matched.models(vv, samples)
    matched = z$matched.models
    scores = z$scores
    if (!is.null(matched)){
      rownames(matched) = seq(1,nrow(matched))
      colnames(matched) = samples
      matched = as.data.frame(matched)
    }

  }
  cat(paste0('Found ', ifelse(is.null(matched), 0, nrow(matched)),
             ' compatible evolutional models\n'))
  return (list(models=vv, matched=matched, scores=scores))

}


#' Enumerate all possible clonal structures for a single sample, not using
#' any test (called absolute method)
#'
#' @description Enumerate all possible clonal structures for a single sample
#' using monoclonal (ie. the primary tumor is originated from a single
#' cancer cell) or polyclonal model (ie. the primary tumor can originate from
#' multi cancer cells)
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
enumerate.clones.absolute <- function(v){
  vv = list() # to hold list of output clonal models
  findParent <- function(v, i){
    if (i > nrow(v)){
      #debug
      #print(v)
      vv <<- c(vv, list(v))
    }else{
      vaf = v[i,]$vaf
      if (!is.na(v[i,]$parent && v[i,]$parent == -1)){# root
        vx = v
        findParent(vx, i+1)
      }else if (vaf == 0){
        vx = v
        vx$parent[i] = NA
        findParent(vx, i+1)
        #}else if (is.na(v$parent[i]) || v$parent[i] != -1){ # not the root, no parent yet, look for parent
      }else{
        #for (j in 1:(i-1)){
        for (j in 1:nrow(v)){
          # asign cluster in row j as parent of cluster in row i
          if (i != j && v$free[j] >= vaf &&
                ifelse(is.na(v$parent[j]), TRUE, v$parent[j] != v$lab[i])){
            vx = v
            vx$free[j] = vx$free[j] - vaf
            vx$occupied[j] = vx$occupied[j] + vaf
            vx$parent[i] = vx[j,]$lab
            # debug
            #cat(i, '<-', j, 'vaf=', vaf, '\n')
            findParent(vx, i+1)
          }
        }
      }
    }
  }

  # if normal sample (0) is included, the normal sample
  # will be root (polyclonal model), otherwise
  if (v[1,]$lab == 0){
    v[1,]$parent = -1
    findParent(v, 2)
  }else{
    max.vaf = max(v$vaf)
    roots = rownames(v)[v$vaf == max.vaf]
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


