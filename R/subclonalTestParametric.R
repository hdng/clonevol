#variants is a data.frame with columns $cluster $primary $relapse
#and one row for each variant

#assumes normal parameterization but should be updated to include beta, beta-binomial, gamma
#perhaps an option to mathemagically pick the best model for each cluster
#also want it to have bootstrap=TRUE or FALSE option, where FALSE just uses the theoretical distribution of the bootstrap

subclonalTestParametric <- function(variants, clusterCol='cluster', vafCol='primary', parentCluster, subClusters, model='normal', alpha=.05, iterations=1000){
  #Takes variants file and performs parametric bootstrap resampling to find distribution of cluster means and variances

  #Variants file is a data.frame with column names $cluster, $primary, and $relapse and one row for each variant

  #Model can be 'normal', 'zero-normal', 'beta', or 'beta-binomial'
  #Zero-normal is the zero truncated normal distribution.
  if(!(model %in% c("normal","zero-normal","beta","beta-binomial"))){
    #checks to make sure the input variable 'model' is an option
    print("Model type not currently an option. Must be 'normal', 'zero-normal', 'beta', or 'beta-binomial'.")
    stop()
  }
  
  #Uses bootstrap of clusters to compute value of D = mean(parent.cluster)-sum(mean(sub.clusters))
  calculateD <- function(means,vars,ns,nclusters){
    clusterMeans <- vector()
    clusterMeans[1] <- mean(rnorm(n=ns[1], mean=means[1], sd=sds[1]))
    for(i in 2:nclusters){
      clusterMeans[i] <- mean(rnorm(n=ns[i], mean=means[i], sd=sds[i]))
      }
    return(clusterMeans[1]-sum(clusterMeans[2:nclusters]))
  }

  #import data
  variants <- as.data.frame(read.table(variants,header=T))
  #force data types
  variants[,1] <- as.character(variants[,clusterCol])

  #find mean, variance, and number of variants of each cluster
  #these are used to calculate MLEs
  means <- vector(); sds <- vector(); ns <- vector(); count <- 0	
  allClusters <- c(parentCluster, subClusters); nclusters <- length(allClusters)
  for(clusterID in allClusters){
    count <- count + 1
    means[count] <- mean(variants[variants[,clusterCol] == clusterID, vafCol])
    sds[count] <- sd(variants[variants[,clusterCol] == clusterID, vafCol])
    ns[count] <- length(variants[variants[,clusterCol] == clusterID, vafCol])
  }

  #find maximum likelihood estimators (MLEs), depending on the model option
  #MLEs are use as the parameters for the parameterized bootstrap
  if(model %in% c("normal","zero-normal")){
    #use same MLEs for normal and zero-normal models
    normalMLEsMeans <- means
    normalMLEsSDs <- sds
  }

  if(model == "beta"){
    #use method of moments to estimate MLEs
    #for sample mean m and sample variance v, if v < m*(1-m):
    #alphahat = m*((m-m^2)/v - 1)
    #betahat = (1-m)*((m-m^2)/v - 1)
    if(all(v < m*(1-m))){
      betaMLEsAlphas <- means*((means - means*means)/(sds*sds) - 1)
      betaMLEsBetas <- (1 - means)*((means - means*means)/(sds*sds) - 1)
    }
  }

  if(model == "beta-binomial"){
    print("Beta-binomial not yet implemented. Working on it...")
    stop()
  }
       
  #calculate large number of D values by resampling from each cluster
  t <- vector(); t[1] <- proc.time()[3]
  D <- vector()
  for(i in 1:iterations){
    D[i] <- calculateD(means,sds,ns,nclusters)
  }
     
  #what is the probability that D<0?
  p <- mean(D < 0)
  #what is the 1-alpha confidence interval of D?
  lb <- quantile(D, probs=alpha/2) #lowerBound
  ub <- quantile(D, probs=1-alpha/2) #upperBound
  
  print(paste0("Confidence interval of D = (",lb,",",ub,")"))
  print(paste0("Probability D<0 is ",p))
  print(paste0("D has bootstrap mean ",mean(D)," and sd ",sd(D)))
  
  t[2] <- proc.time()[3]

  #theoretical distribution of D
  print(paste0("D has theoretical distribution Normal(",means[1]-sum(means[2:nclusters]),",",sqrt(sum(sds*sds/ns)),")"))
  print(paste0("Theoretical P(D<0)=",pnorm(0,means[1]-sum(means[2:nclusters]),sqrt(sum(sds*sds/ns))),"."))
  print(paste0("Theoretical 95% CI around mean(D) is (",qnorm(.025,means[1]-sum(means[2:nclusters]),sqrt(sum(sds*sds/ns))),",",qnorm(.975,means[1]-sum(means[2:nclusters]),sqrt(sum(sds*sds/ns))),")."))

  t[3] <- proc.time()[3]
  print(c(t[2]-t[1],t[3]-t[2]))
        
  output <- list(probability=p, lowerBound=lb, upperBound=ub)
  return(output)
}

subclonalTestParametric("AML31.variants", clusterCol='cluster', vafCol='primary', parentCluster='1', subClusters=c('2','4'), model='normal', alpha=.05, iterations=1000)
