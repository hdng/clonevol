#' Generate and calculate bootstrap means for all clusters
#'
#' @param variants: data frame of variants with cluster assignments and VAF
#' of samples. Columns are c('cluster', 'sample1', 'sample2', ....)
#' @param zero.sample: The sample of zero vaf (to use to compare with other
#' clusters to determine if the cluster should be considered zero VAF, and
#' not included in the models)
#' @param model: specifies the statistical model used in parametric bootstrap resampling (can be normal, beta, or beta-binomial)

#Last update: Steven Mason Foltz 2015-03-10 13:45
#Original by Ha X. Dang
#SMF added parametric bootstrap functionality

generate.boot.parametric <- function(variants, cluster.col.name='cluster', vaf.col.names=NULL, vaf.in.percent=TRUE, num.boots=1000, zero.sample=NULL, model=NULL){
  cat('Generating boostrap samples...\n')
  boot.means = NULL
  if (is.null(vaf.col.names)){
    vaf.col.names = setdiff(colnames(variants), cluster.col.name)
  }
  
  ### SMF addition begin ###
  if(is.null(model) | !(model %in% c("normal","beta","beta-binomial"))){
    print("User must specify statistical model for parametric bootstrap resampling. Model can be 'normal', 'beta', or 'beta-binomial'.")
    stop()
  }
  ### SMF addition end #####

  #if no cluster or no sample provided, return NULL
  clusters = unique(variants[[cluster.col.name]])
  num.clusters = length(clusters)
  num.samples = length(vaf.col.names)
  if (num.samples == 0 || num.clusters == 0){return(NULL)}
  if (vaf.in.percent){
    variants[,vaf.col.names] = variants[,vaf.col.names]/100.00
  }
  #make separate data frame for each cluster
  v = list()
  for (cl in clusters){
    v1 = variants[variants[[cluster.col.name]]==cl, vaf.col.names]
    v[[as.character(cl)]] = v1
  }

  # debug
  # vv <<- v

  # generate bootstrap samples for each cluster, each sample
  num.variants.per.cluster = table(variants[[cluster.col.name]])
  #print(num.variants.per.cluster)

  boot.means = list()
  clusters = as.character(clusters)
  zeros = c()
  for (vaf.col.name in vaf.col.names){
    #cat('Booting sample: ', vaf.col.name, '\n')
    sample.boot.means = matrix(NA, nrow=num.boots, ncol=num.clusters)
    colnames(sample.boot.means) = clusters
    rownames(sample.boot.means) = seq(1, num.boots)
    for (cl in clusters){
      # learn zero samples from data, if a cluster has median VAF = 0, consider
      # it as a sample generated from true VAF = 0
      vafs = v[[cl]][[vaf.col.name]]
      if (median(vafs)==0){zeros = c(zeros, vafs)}
      
      boot.size = num.variants.per.cluster[cl]

      ### SMF addition begin ###
      #find the mean and standard deviation of the cluster
      #mean and sd are used as parameters in bootstrapping
      this.mean <- mean(vafs)
      this.sd <- sd(vafs)
      
      if(model == "normal"){
        for (b in 1:num.boots){
          #use mean and standard deviation as normal MLEs
          s.mean <- mean(rnorm(n=boot.size, mean=this.mean, sd=this.sd))
          sample.boot.means[b,cl] <- s.mean
        }
        next
      }

      if(model == "beta"){
        #use mean and sd to calculate alpha and beta (method of moments)
        m <- this.mean; var <- this.sd^2
        alpha <- m*((m-m*m)/var-1); beta <- (1-m)*((m-m*m)/var-1)
        for(b in 1:num.boots){
          s.mean <- mean(rbeta(n=boot.size, shape1=alpha, shape2=beta))
          sample.boot.means[b,cl] <- s.mean
        }
        next
      }

      if(model == "beta-binomial"){
        #uses binomial distribution with probability of success drawn from a beta distribution
        #use mean and sd to calculate alpha and beta (method of moments)
        m <- this.mean; var <- this.sd^2
        alpha <- m*((m-m*m)/var-1); beta <- (1-m)*((m-m*m)/var-1)
        for(b in 1:num.boots){
          s.mean <- mean(rbinom(n=boot.size, size=100, prob=rbeta(n=1, shape1=alpha, shape2=beta)))/100
          sample.boot.means[b,cl] <- s.mean
        }
        next
      }

      ### SMF addition end #####

      ### SMF comment begin ###
      #cat('Booting cluster: ', cl, 'boot.size=', boot.size, '\n')
      #for (b in 1:num.boots){
      #  s.mean = mean(sample(v[[cl]][[vaf.col.name]], boot.size, replace=T))
      #  sample.boot.means[b, cl] = s.mean
      #}
      ### SMF comment end #####
    }
    boot.means[[vaf.col.name]] = sample.boot.means
  }

  #generate bootstrap means for zero sample
  if (is.null(zero.sample)){
    zero.sample = zeros
  }
  if (length(zero.sample) > 0){
    zero.sample.boot.means = rep(NA, num.boots)
    zero.sample.boot.size = length(zero.sample)
    for (b in 1:num.boots){
      s.mean = mean(sample(zero.sample, zero.sample.boot.size, replace=T))
      zero.sample.boot.means[b] = s.mean
    }
    boot.means$zero.means = zero.sample.boot.means
  }
  return(boot.means)
}
