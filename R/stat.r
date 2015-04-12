# bootstrap resampling and test

resample.test.greater <- function(x, y, nBoots=10000){
    n.x = length(x)
    n.y = length(y)
    x.bmean = rep(NA, nBoots)
    y.bmean = rep(NA, nBoots)
    t = rep(NA, nBoots)
    for (i in 1:nBoots){
        x.bmean[i] = mean(sample(x, n.x, replace=T))
        y.bmean[i] = mean(sample(y, n.y, replace=T))
        t[i] = x.bmean[i] - y.bmean[i]
    }
    p = (sum(t < 0)+1)/(length(t)+1)
    return(p)
}

resample.test.pooled <- function(x, y, nBoots=10000){
    mean.diff = mean(x) - mean(y)
    z = c(x,y)
    n.x = length(x)
    n.y = length(y)
    labs = c(rep('x',n.x), rep('y', n.y))
    t = rep(NA, nBoots)
    for (i in 1:nBoots){
        labs.i = sample(labs)
        x.mean = mean(z[which(labs.i=='x')])
        y.mean = mean(z[which(labs.i=='y')])
        t[i] = x.mean - y.mean
    }
    p = (sum(t > mean.diff)+1)/(length(t)+1)
    return(p)
}


xxx <- function(){

    v = crc12.variants
    vaf.col.names = grep('WGS_VAF', colnames(v), value=T)
    v = v[, c('cluster', vaf.col.names)]
    clone.vafs = estimate.clone.vaf(v, 'cluster', vaf.col.names)

    v = v[,c('cluster','met7.WGS_VAF')]

    x = v[v$cluster==1,2]
    y = v[v$cluster==9,2]
    resample.test.greater(x, y, nBoots=10000)
    resample.test.pooled(x, y, nBoots=10000)
    t.test(x, y, alternative='greater')

    v = aml31.variants
    x = v[v$cluster==6,2]
    y = rep(0, 100)
    resample.test(x, y, nBoots=10000)
    resample.test.pooled(x, y, nBoots=10000)
    t.test(x, y, alternative='greater')



}



#' Generate and calculate bootstrap means for all clusters
#' Depricated!
#'
#' @param variants: data frame of variants with cluster assignments and VAF
#' of samples. Columns are c('cluster', 'sample1', 'sample2', ....)
#' @param zero.sample: The sample of zero vaf (to use to compare with other
#' clusters to determine if the cluster should be considered zero VAF, and
#' not included in the models)
generate.boot.nonparametric <- function(variants, cluster.col.name='cluster',
                          vaf.col.names=NULL, vaf.in.percent=TRUE,
                          num.boots=1000,
                          zero.sample=NULL){
    cat('Generating boostrap samples...\n')
    boot.means = NULL
    if (is.null(vaf.col.names)){
        vaf.col.names = setdiff(colnames(variants), cluster.col.name)
    }

    # if no cluster or no sample provided, return NULL
    clusters = unique(variants[[cluster.col.name]])
    num.clusters = length(clusters)
    num.samples = length(vaf.col.names)
    if (num.samples == 0 || num.clusters == 0){return(NULL)}
    if (vaf.in.percent){
        variants[,vaf.col.names] = variants[,vaf.col.names]/100.00
    }
    # make separate data frame for each cluster
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
            #cat('Booting cluster: ', cl, 'boot.size=', boot.size, '\n')
            for (b in 1:num.boots){
                s.mean = mean(sample(v[[cl]][[vaf.col.name]], boot.size, replace=T))
                sample.boot.means[b, cl] = s.mean
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

#' Calculate CI for a cluster in a sample
boot.vaf.ci <- function(boot, sample, cluster, alpha=0.05){
    vaf.means = boot[[sample]][,as.character(cluster)]
    mean.vaf = mean(vaf.means)
    upper.vaf = quantile(vaf.means, 1-alpha/2)
    lower.vaf = quantile(vaf.means, alpha/2)
    upper.vaf.fmt = sprintf('%0.2f%%', upper.vaf)
    lower.vaf.fmt = sprintf('%0.2f', lower.vaf)
    vaf.ci.str = paste0(lower.vaf.fmt, ' - ', upper.vaf.fmt)
    return(list(vaf.mean=mean.vaf,
                vaf.lower=lower.vaf, vaf.upper=upper.vaf,
                vaf.ci.str=vaf.ci.str))
}

#' Subclonal test H0: mean_X1 + mean_X2 <= mean_X
#'
#' @param boot: output of generate.boot(), if NULL, test using VAFs
#' @param vaf.col.name: name of the VAF column
#' @param parent.cluster: the cluster of mutations representing the parent
#' clone
#' @param sub.clusters: vector of sub clusters to test if they can be
#' all directly arise from the parent clone; if NULL, the parent.cluster
#' will be tested against 0 (this will help determine if this (parent)clone
#' represents with high enough cell frac to report)
#' @param cdf: Clonal data frame with VAF
#'
#' @description Return a list of p-value, confidence intervals of X-(X1+X2),etc.
#' if p-value is small, reject H0, otherwise, not enough evidence to reject H0
#'
subclonal.test <- function(vaf.col.name, parent.cluster, sub.clusters=NULL,
                           boot=NULL, cdf=NULL, min.cluster.vaf=0, alpha=0.05){
    # debug
    #cat('subclonal.test: sample=', vaf.col.name, 'parent.cluster=', parent.cluster,
    #    'sub.clusters=', paste(sub.clusters, collapse=','),'\n')

    if (is.null(boot) || length(boot) == 0){
        # test using absolute values of VAF
        # debug
        #cat('No booting! absolute value comparison.\n')
        if (is.null(sub.clusters)){
            mean.free.vaf = cdf$vaf[cdf$lab==parent.cluster]
            p = ifelse(mean.free.vaf > min.cluster.vaf, 1, 0)
        }else{
            mean.free.vaf = (cdf$vaf[cdf$lab==parent.cluster] -
                sum(cdf$vaf[cdf$lab %in% sub.clusters]))
            p = ifelse(mean.free.vaf >= 0, 1, 0)
        }
        #debug
        #if (parent.cluster == 'c1' && (all(c('c1a', 'c1aC1') %in% sub.clusters))){
        #  cat('mean.free.vaf=', mean.free.vaf, 'p=', p, '\n')
        #stop()
        #}
        #zz <<- mean.free.vaf
        free.vaf.mean=mean.free.vaf
        upper.free.vaf = NA
        lower.free.vaf = NA
        upper.free.vaf.fmt = NA
        lower.free.vaf.fmt = NA
        confident.level = NA
        confident.level.non.negative = NA
        free.vaf.ci.str = NA

    }else{
        num.boots = nrow(boot[[1]])
        if (num.boots == 0){return(NULL)}
        if (is.null(sub.clusters)){
            #free.vaf = boot[[vaf.col.name]][,parent.cluster]
            # -boot$zero.means            
            free.vaf = boot[[vaf.col.name]][,parent.cluster]
        }else{
            free.vaf = apply(boot[[vaf.col.name]], 1,
                             function(row) (row[parent.cluster] -
                                                sum(row[sub.clusters])))
        }
        zz <<- free.vaf
        zero.vaf = ifelse(is.null(sub.clusters), min.cluster.vaf, 0)
        p = sum(free.vaf > zero.vaf)/length(free.vaf)
        mean.free.vaf = mean(free.vaf)
        upper.free.vaf = quantile(free.vaf, 1-alpha/2)
        lower.free.vaf = quantile(free.vaf, alpha/2)
        upper.free.vaf.fmt = sprintf('%0.2f%%', upper.free.vaf)
        lower.free.vaf.fmt = sprintf('%0.2f', lower.free.vaf)
        confident.level = 1 - alpha
        confident.level.non.negative = ifelse(lower.free.vaf >= 0,
            confident.level, sum(free.vaf >= 0 &
            free.vaf <= upper.free.vaf)/length(free.vaf)
        )
        free.vaf.ci.str = ifelse(lower.free.vaf > 0,
                                 paste0(lower.free.vaf.fmt, ' - ',
                                        upper.free.vaf.fmt),
                                 paste0('0 - ', upper.free.vaf.fmt))
    }
    #debug
    # cat('p-value =', p, '\n')
    #cat('CI =', lower.free.vaf, '-', upper.free.vaf, 'free.mean=', mean.free.vaf, '\n')
    
    return(list(free.vaf.ci=free.vaf.ci.str,
                free.vaf.mean=mean.free.vaf,
                free.vaf.lower=lower.free.vaf,
                free.vaf.upper=upper.free.vaf,
                free.vaf.confident.level=confident.level,
                free.vaf.confident.level.non.negative=confident.level.non.negative,
                p.value=p))
}



testttt <- function(){
    variants = crc8.variants
    vaf.col.names = grep('.VAF', colnames(variants), value=T, fixed=T)
    vaf.col.names = vaf.col.names[!grepl('PBMC', vaf.col.names)]
    clone.vafs = estimate.clone.vaf(variants, 'cluster', vaf.col.names)

    sample = 'C.VAF'
    v = make.clonal.data.frame(vafs=clone.vafs[[sample]],
                               labels=clone.vafs$cluster,
                               #founding.label='1'
                               )
    boot = generate.boot(variants, vaf.col.names=vaf.col.names, num.boots=1000)
    source('R/clonevol.r');

    xb = enumerate.clones(v, sample, variants, boot=boot, founding.cluster='1');
    draw.sample.clones.all(xb, paste0('test-out/CRC8/', sample))

    x = enumerate.clones.absolute(v)
    draw.sample.clones.all(x, 'test-out/CRC8-absolute')
    draw.sample.clones(rescale.vaf(xb[[1]]), cell.frac.ci=T)


    variants = crc12.variants
    #variants[variants$cluster %in% c(5,6,7),]$cluster = 5
    #variants[variants$cluster > 5,]$cluster = variants[variants$cluster > 5,]$cluster - 2
    vaf.col.names = grep('WGS_VAF', colnames(variants), value=T)
    vaf.col.names = vaf.col.names[!grepl('normal', vaf.col.names)]
    clone.vafs = estimate.clone.vaf(variants, 'cluster', vaf.col.names)
    x = infer.clonal.models(variants=variants, vaf.col.names=vaf.col.names,
                            subclonal.test='bootstrap', num.boots=1000,
                            founding.cluster=1, min.cluster.vaf=0)
    plot.clonal.models(x$models,
                       out.dir='test-out/bootstrap-test',
                       matched=x$matched,
                       out.format='png', overwrite.output=T,
                       scale.monoclonal.cell.frac=TRUE,
                       cell.frac.ci=T,
                       tree.node.shape='circle')
    for (s in vaf.col.names){
        draw.sample.clones.all(x$models[[s]], paste0('test-out/bootstrap-test/', s))
    }


    variants = crc12.variants.new
    variants.xeno = variants[, grepl('_XT1', colnames(variants))]
    xeno.vaf.col.names = grep('_XT1.VAF', colnames(variants), value=T, fixed=T)

    #exlude xeno
    variants = variants[, !grepl('_XT1', colnames(variants))]

    #variants[variants$cluster == 'c9',]$Li3_XT1.VAF = variants[variants$cluster == 'c9',]$Li3_XT1.VAF - 18.7

    # remove some small clusters that conflict!
    variants = variants[!(variants$cluster %in% c('c9', 'c11','c23')),]
    model = 'non-parametric'
    out.dir = paste0('test-out/CRC12-new-', model)
    #variants[variants$cluster %in% c(5,6,7),]$cluster = 5
    #variants[variants$cluster > 5,]$cluster = variants[variants$cluster > 5,]$cluster - 2
    vaf.col.names = grep('.VAF', colnames(variants), value=T, fixed=T)

    # remove normal sample
    vaf.col.names = vaf.col.names[!grepl('PBMC', vaf.col.names)]

    clone.vafs = estimate.clone.vaf(variants, 'cluster', vaf.col.names)

    x = infer.clonal.models(variants=variants, vaf.col.names=vaf.col.names,
                            subclonal.test='bootstrap',
                            subclonal.test.model=model,
                            num.boots=1000,
                            founding.cluster='c1', min.cluster.vaf=0.025,
                            p.value.cutoff=0.1)
    # xeno model merging
    #matched = x$matched
    #matched$index$C_XT1.VAF = 1
    #matched$index$Li3_XT1.VAF = 1
    #matched$index$Li2_XT1.VAF = 1
    #colnames(x$matched$index) = vaf.col.names

    plot.clonal.models(x$models,
                       out.dir=out.dir,
                       matched=x$matched,
                       variants=variants,
                       box.plot=T,
                       out.format='pdf', overwrite.output=T,
                       scale.monoclonal.cell.frac=TRUE,
                       cell.frac.ci=F,
                       tree.node.shape='circle',
                       tree.node.size=35,
                       max.num.models.to.plot=10,
                       width=15, height=15,
                       )

    dir.create(out.dir)
    for (s in vaf.col.names){
        draw.sample.clones.all(x$models[[s]], paste0(out.dir, '/', s, '.tree'),
                               object.to.plot='tree')
    }

    variants = crc8.variants
    out.dir = 'test-out/CRC8'
    variants = variants[!(variants$cluster %in% c(4,5)),]
    vaf.col.names = grep('.VAF', colnames(variants), value=T, fixed=T)
    vaf.col.names = vaf.col.names[!grepl('PBMC', vaf.col.names)]
    clone.vafs = estimate.clone.vaf(variants, 'cluster', vaf.col.names)
    x = infer.clonal.models(variants=variants, vaf.col.names=vaf.col.names,
                            subclonal.test='bootstrap', num.boots=1000,
                            founding.cluster='1', min.cluster.vaf=0.025,
                            p.value.cutoff=0.25)
    plot.clonal.models(x$models,
                       out.dir=out.dir,
                       matched=x$matched,
                       variants=variants,
                       box.plot=T,
                       out.format='png', overwrite.output=T,
                       scale.monoclonal.cell.frac=TRUE,
                       cell.frac.ci=F,
                       tree.node.shape='circle',
                       tree.node.size=55,
                       max.num.models.to.plot=1,
                       width=11, height=10)
    dir.create(out.dir)
    for (s in vaf.col.names){
        draw.sample.clones.all(x$models[[s]], paste0(out.dir, '/', s, '.tree'),
                               object.to.plot='tree')
    }


    variants = aml31.variants
    vaf.col.names = c('Tumor', 'Relapse')
    clone.vafs = estimate.clone.vaf(variants, 'cluster', vaf.col.names)
    x = infer.clonal.models(variants=variants, vaf.col.names=vaf.col.names,
                            subclonal.test='bootstrap', num.boots=1000,
                            founding.cluster=1)
    plot.clonal.models(x$models, variants=variants,
                       out.dir='test-out/bootstrap-test-aml31',
                       matched=x$matched,
                       out.format='pdf', overwrite.output=T,
                       tree.node.shape='circle',
                       scale.monoclonal.cell.frac=TRUE,
                       cell.frac.ci=T,
                       box.plot=T,
                       width=10)
    for (s in vaf.col.names){
        draw.sample.clones.all(x$models[[s]],
                               paste0('test-out/bootstrap-test-aml31/', s))
    }

    variants = aml1.variants
    vaf.col.names = c('Tumor', 'Relapse')
    clone.vafs = estimate.clone.vaf(variants, 'cluster',
                                    vaf.col.names, vaf.in.percent=F)
    x = infer.clonal.models(variants=variants, vaf.col.names=vaf.col.names,
                            vaf.in.percent=F,
                            subclonal.test='bootstrap', num.boots=10000,
                            min.cluster.vaf=0.01,
                            founding.cluster=1)
    plot.clonal.models(x$models,
                       out.dir='test-out/bootstrap-test-aml1',
                       matched=x$matched,
                       out.format='png', overwrite.output=T,
                       tree.node.shape='square',
                       scale.monoclonal.cell.frac=TRUE,
                       width=10,
                       cell.frac.ci=T)
    for (s in vaf.col.names){
        draw.sample.clones.all(x$models[[s]],
                               paste0('test-out/bootstrap-test-aml1/', s))
    }

    c = msclc3588.ccf
    c[,2:4] = c[,2:4]/2
    clones = c
    clones = clones[order(clones$CCF1, decreasing=T),]
    clones$cluster = as.character(clones$cluster)
    #clones = clones[,1:3]
    colnames(clones) = c('cluster', 'Primary', 'Met1', 'Met2')
    x = infer.clonal.models(c=clones)
    plot.clonal.models(x$models,
                       out.dir='test-out/MSCLC3588',
                       matched=x$matched,
                       out.format='png', overwrite.output=T,
                       tree.node.shape='square',
                       tree.node.size = 40,
                       scale.monoclonal.cell.frac=TRUE,
                       cell.frac.ci=F)



    c = msclc3151.ccf.touched
    clones = c
    clones$cluster = as.character(clones$cluster)
    #clones = clones[,1:3]
    colnames(clones) = c('cluster', 'Primary', 'Lymph_Met', 'Liver_Met')
    x = infer.clonal.models(c=clones)
    plot.clonal.models(x$models,
                       out.dir='test-out/MSCLC3151',
                       matched=x$matched,
                       out.format='png', overwrite.output=T,
                       tree.node.shape='square',
                       tree.node.size = 40,
                       scale.monoclonal.cell.frac=TRUE,
                       cell.frac.ci=F)

    c = msclc984.ccf
    clones = c
    clones$cluster = as.character(clones$cluster)
    colnames(clones) = c('cluster', 'Primary', 'Lymph_Met')
    x = infer.clonal.models(c=clones)
    plot.clonal.models(x$models,
                       out.dir='test-out/MSCLC984',
                       matched=x$matched,
                       out.format='png', overwrite.output=T,
                       tree.node.shape='square',
                       tree.node.size = 40,
                       scale.monoclonal.cell.frac=TRUE,
                       cell.frac.ci=F)

    c = msclc984.ccf
    ggplot(c, aes(x=CCF1, y=CCF2)) + geom_text(aes(label=cluster))



    #### CRC12 Mar. 30, 2015
    variants = read.table('samples/CRC12-finer.tsv', header=T, sep='\t', quote='', stringsAsFactors=F)
    #variants$cluster = paste0('c', variants$cluster)
    variants.xeno = variants[, grepl('_XT1', colnames(variants))]
    xeno.vaf.col.names = grep('_XT1.VAF', colnames(variants), value=T, fixed=T)

    #exlude xeno
    variants = variants[, !grepl('_XT1', colnames(variants))]

    # remove cluster 2
    variants = variants[variants$cluster != 2,]

    model = 'normal'
    out.dir = paste0('test-out/CRC12-March30-', model)
    #variants[variants$cluster %in% c(5,6,7),]$cluster = 5
    #variants[variants$cluster > 5,]$cluster = variants[variants$cluster > 5,]$cluster - 2
    vaf.col.names = grep('.VAF', colnames(variants), value=T, fixed=T)

    # remove normal sample
    vaf.col.names = vaf.col.names[!grepl('PBMC', vaf.col.names)]

    clone.vafs = estimate.clone.vaf(variants, 'cluster', vaf.col.names)

    x = infer.clonal.models(variants=variants, vaf.col.names=vaf.col.names,
                            subclonal.test='bootstrap',
                            subclonal.test.model=model,
                            num.boots=5000,
                            founding.cluster='1', min.cluster.vaf=0.01,
                            p.value.cutoff=0.01)
    # xeno model merging
    #matched = x$matched
    #matched$index$C_XT1.VAF = 1
    #matched$index$Li3_XT1.VAF = 1
    #matched$index$Li2_XT1.VAF = 1
    #colnames(x$matched$index) = vaf.col.names

    plot.clonal.models(x$models,
                       out.dir=out.dir,
                       matched=x$matched,
                       variants=variants,
                       box.plot=T,
                       out.format='pdf', overwrite.output=T,
                       scale.monoclonal.cell.frac=TRUE,
                       cell.frac.ci=F,
                       tree.node.shape='circle',
                       tree.node.size=35,
                       max.num.models.to.plot=10,
                       width=15, height=15,
    )

    for (s in vaf.col.names){
        draw.sample.clones.all(x$models[[s]], paste0(out.dir, '/', s, '.polygon'),
                               object.to.plot='polygon')
    }

    pdf(paste0(out.dir, '/variants.box.pdf'), width=7, height=15)
    variant.box.plot(variants, vaf.col.names = vaf.col.names, sample.title.size=10)
    dev.off()


    #### CRC8 Mar. 30, 2015
    variants = read.table('samples/CRC8-finer.tsv', header=T, sep='\t', quote='', stringsAsFactors=F)

    # remove cluster 1,5
    variants = variants[!(variants$cluster %in% c(1,5)),]

    model = 'non-parametric'
    out.dir = paste0('test-out/CRC8-March30-', model)
    #variants[variants$cluster %in% c(5,6,7),]$cluster = 5
    #variants[variants$cluster > 5,]$cluster = variants[variants$cluster > 5,]$cluster - 2
    vaf.col.names = grep('.VAF', colnames(variants), value=T, fixed=T)
    vaf.col.names = vaf.col.names[!grepl('PBMC', vaf.col.names)]

    clone.vafs = estimate.clone.vaf(variants, 'cluster', vaf.col.names)
    sample = 'CRC8_237_C_20120525.VAF'
    v = make.clonal.data.frame(vafs=clone.vafs[[sample]],
                               labels=clone.vafs$cluster,
                               founding.label='1'
                            )
    boot = generate.boot(variants, vaf.col.names=vaf.col.names, num.boots=1000)
    source('R/clonevol.r')
    xb = enumerate.clones(v, sample, variants, boot=boot, founding.cluster='2')
    draw.sample.clones.all(xb, paste0('test-out/CRC8/', sample))


    e = enumerate.clones(v = clone.vafs, sample = '', variants = variants, founding.cluster = 1)

    x = infer.clonal.models(variants=variants, vaf.col.names=vaf.col.names,
                            subclonal.test='bootstrap',
                            subclonal.test.model=model,
                            num.boots=1000,
                            founding.cluster='2', min.cluster.vaf=0.025,
                            p.value.cutoff=0.01)
    # xeno model merging
    #matched = x$matched
    #matched$index$C_XT1.VAF = 1
    #matched$index$Li3_XT1.VAF = 1
    #matched$index$Li2_XT1.VAF = 1
    #colnames(x$matched$index) = vaf.col.names

    plot.clonal.models(x$models,
                       out.dir=out.dir,
                       matched=x$matched,
                       variants=variants,
                       box.plot=T,
                       out.format='pdf', overwrite.output=T,
                       scale.monoclonal.cell.frac=TRUE,
                       cell.frac.ci=T,
                       tree.node.shape='circle',
                       tree.node.size=35,
                       max.num.models.to.plot=10,
                       width=15, height=15,
    )

    for (s in vaf.col.names){
        draw.sample.clones.all(x$models[[s]], paste0(out.dir, '/', s, '.polygon'),
                               object.to.plot='polygon')
    }

    pdf(paste0(out.dir, '/variants.box.pdf'), width=7, height=15)
    variant.box.plot(variants, vaf.col.names = vaf.col.names, sample.title.size=10)
    dev.off()

}
