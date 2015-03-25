#' Generate and calculate bootstrap means for all clusters
#'
#' @param variants: data frame of variants with cluster assignments and VAF
#' of samples. Columns are c('cluster', 'sample1vaf', 'sample2vaf', ...,
#' 'sample1count', 'sample2count', ....). Corresponding VAF and count columns
#' should appear in the same relative order unless both are specified. 
#' @param cluster.col.name: name of column containing cluster information.
#' @param vaf.col.names: names of columns containing VAFs for each sample. If
#' NULL, columns ending in vaf will be used. Default NULL.
#' @param vaf.col.name.suffix: when VAF columns are not specified, this pattern
#' is matched at the end of column names to define VAF columns. Default vaf.
#' @param count.col.names: names of columns containing count/depth for each
#' sample. If NULL, columns ending in count will be used. Not needed if
#' weighted parameter is FALSE. Default count.
#' @param count.col.name.suffix: when count columns are not specified, this
#' pattern is matched at the end of column names to define count columns.
#' Default NULL.
#' @param vaf.in.percent: If TRUE, VAFs will be converted to proportions
#' between 0 and 1. Default TRUE.
#' @param num.boots: Number of times to resample. Default 1000.
#' @param model: specifies the statistical model used in bootstrap resampling. 
#' Model can be normal, normal-truncated, beta, binomial, beta-binomial, or
#' non-parametric.
#' @param weighted: If TRUE, weights points proportionally to read count.
#' Default TRUE.
#' @param zero.sample: The sample of zero vaf (to use to compare with other
#' clusters to determine if the cluster should be considered zero VAF, and
#' not included in the models)  

#Last update: Steven Mason Foltz 2015-03-23
#Original by Ha X. Dang
#SMF added weighted parametric bootstrap functionality

generate.boot <- function(variants,
                          cluster.col.name='cluster',
                          vaf.col.names=NULL,
                          vaf.col.name.suffix='vaf',
                          count.col.names=NULL,
                          count.col.name.suffix='count',
                          vaf.in.percent=TRUE,
                          num.boots=1000,
                          model=NULL,
                          weighted=TRUE,
                          zero.sample=NULL){

    #check that the model is not NULL
    if(is.null(model) | !(model %in% c("normal", "normal-truncated", "beta",
        "binomial", "beta-binomial", "non-parametric"))){
        stop("User must specify statistical model for parametric bootstrap resampling. Model can be 'normal', 'normal-truncated', 'beta', 'binomial', 'beta-binomial', or 'non-parametric'.")
    }

    #check that the weighted parameter is logical
    if(!is.logical(weighted)){
        stop("'weighted' parameter must be TRUE or FALSE. Default is TRUE.")
    }

    #check VAF and count column specifications
    other.col.names = setdiff(colnames(variants), cluster.col.name)
    if(is.null(vaf.col.names)){
        vaf.col.names = other.col.names[grepl(paste0(vaf.col.name.suffix,"$"),
        other.col.names,ignore.case=T)]
    }
    if(is.null(count.col.names) & weighted){
        count.col.names = other.col.names[grepl(paste0(count.col.name.suffix,
        "$"),other.col.names,ignore.case=T)]
    }
    if(weighted & length(vaf.col.names)!=length(count.col.names)){
        stop("Number of VAF columns differs from number of count columns. If column names for VAF and count are not specified, column names ending in 'vaf' and 'count' (not case sensitive) are used by default.")
    }
 
    #if no cluster or no sample provided, return NULL
    clusters = unique(variants[[cluster.col.name]])
    num.clusters = length(clusters)
    num.samples = length(vaf.col.names)
    if (num.samples == 0 || num.clusters == 0){return(NULL)}
    if (vaf.in.percent){
        variants[,vaf.col.names] = variants[,vaf.col.names]/100.00
    }

    cat('Generating boostrap samples...\n')
    boot.means = NULL
        
    #make separate data frame for each cluster
    v = list()
    for (cl in clusters){
        v1 = variants[variants[[cluster.col.name]]==cl, c(vaf.col.names,
        count.col.names)]
        v[[as.character(cl)]] = v1
    }

    # generate bootstrap samples for each cluster, each sample
    num.variants.per.cluster = table(variants[[cluster.col.name]])
    #print(num.variants.per.cluster)

    boot.means = list()
    clusters = as.character(clusters)
    zeros = c()
    for (col.name in 1:length(vaf.col.names)){
        vaf.col.name = vaf.col.names[col.name]
        count.col.name = count.col.names[col.name]
        #cat('Booting sample: ', vaf.col.name, '\n')
        sample.boot.means = matrix(NA, nrow=num.boots, ncol=num.clusters)
        colnames(sample.boot.means) = clusters
        rownames(sample.boot.means) = seq(1, num.boots)

        for (cl in clusters){
            boot.size = num.variants.per.cluster[cl]
            vafs = v[[cl]][[vaf.col.name]]

            # learn zero samples from data,
            # if a cluster has median VAF = 0, consider
            # it as a sample generated from true VAF = 0
            if (median(vafs)==0){zeros = c(zeros, vafs)}
          
            #find the mean and standard deviation of the cluster
            #mean and sd are used as parameters in bootstrapping
            if(weighted){ #weighted sum and sd
                counts = v[[cl]][[count.col.name]]
                #counts = log(counts+1)
                this.mean = (1/sum(counts))*sum(vafs*counts)
                this.sd = sqrt((sum(counts)/(sum(counts)^2-sum(counts^2)))*
                sum(counts*(vafs-this.mean)^2))
            }
            else{ #not weighted
                counts = rep(1,boot.size)
                this.mean = mean(vafs)
                this.sd = sd(vafs)
            }
            
            #uses normal - could produce values below 0 or above 1 (bad)
            if(model == "normal"){
                for (b in 1:num.boots){
                    #use mean and standard deviation as normal MLEs
                    s.mean = mean(rnorm(n=boot.size,mean=this.mean,sd=this.sd))
                    sample.boot.means[b,cl] = s.mean
                }
            }

            #uses zero-one truncated Normal distribution
            else if(model == "normal-truncated"){
                library(truncnorm) #use truncnorm library
                for (b2 in 1:num.boots){ #b2 since b in rtruncnorm()
                    #use mean and standard deviation as normal MLEs
                    s.mean = mean(rtruncnorm(n=boot.size,a=0,b=1,mean=this.mean,
                        sd=this.sd))
                    sample.boot.means[b2,cl] = s.mean
                }
            }

            else if(model == "beta"){
                #use mean and sd to calculate alpha and beta (method of moments)
                m = this.mean; var = this.sd^2
                alpha = m*((m-m*m)/var-1); beta = (1-m)*((m-m*m)/var-1)
                for(b in 1:num.boots){
                    s.mean = mean(rbeta(n=boot.size, shape1=alpha, shape2=beta))
                    sample.boot.means[b,cl] = s.mean
                }
            }

            else if(model == "binomial"){
                #use mean to define probability of success in 100 trials
                for(b in 1:num.boots){
                    s.mean = mean(rbinom(n=boot.size,size=100,prob=this.mean))
                    sample.boot.means[b,cl] = s.mean/100
                }
            }

            else if(model == "beta-binomial"){
                #binomial with probability drawn from a beta distribution
                #use mean and sd to calculate alpha and beta (method of moments)
                m = this.mean; var = this.sd^2
                alpha = m*((m-m*m)/var-1); beta = (1-m)*((m-m*m)/var-1)
                for(b in 1:num.boots){
                    beta.probs = rbeta(n=boot.size,shape1=alpha,shape2=beta)
                    s.mean = mean(rbinom(n=boot.size,size=100,prob=beta.probs))
                    sample.boot.means[b,cl] = s.mean/100
                }
            }

            else { #if(model == "non-parametric"){
                #cat('Booting cluster: ', cl, 'boot.size=', boot.size, '\n')
                for (b in 1:num.boots){
                    s.mean = mean(suppressWarnings(
                    sample(v[[cl]][[vaf.col.name]], boot.size, replace=T,
                    prob=counts)))
                    sample.boot.means[b, cl] = s.mean
                }
            }
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
