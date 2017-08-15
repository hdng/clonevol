#!/usr/bin/env Rscript

USAGE = 'thisfile <subpopulations.truth.tsv> <parental-constrainted:T/F> <depth> <num.trees> <num.replicates> <error.rate> <cn.rate> <rand.seed> <out.dir> <prefix>\n'

# Given M clones and their subpopulation fraction (cancer cell fraction) in N samples
# this script generate random trees and simulate read counts for variants in each clone
# considering sequencing error rate and copy number variations

library(clonevol)

#' Generate n variants from true.vafs across samples
#' @param clone: An integer number indicating clone ID
#' @param n: Number of variants
#' @param sample.names: vector of sample name strings
#' @param true.vafs: vector of true VAFs (0-1)
#' @param mean.depth: mean depth of sequencing
#' @param total.reads: total number of reads to draw from
#' @param err.rate: sequencing error rate (0-1)
#' @param cn.rate: a rate (0-1) at which a variant
#' has undetectable copy number event (diff. from 2) that
#' affects VAF estimate.
# 
generate.variants <- function(clone, n, sample.names, true.vafs, mean.depth=100,
                                total.reads=100000, err.rate=0.01, cn.rate=0.01){
    v = data.frame(clone=rep(clone,n), stringsAsFactors=F)
    # generate reads for each sample
    for (i in 1:length(sample.names)){
        s = sample.names[i]
        vaf = true.vafs[i]
        # generate depths for each variant following Poisson dist.
        depths = rpois(n, mean.depth)
        sim.vafs = c()
        sim.ref.cnts = c()
        sim.var.cnts = c()
        # decide (randomly) what reads will be cn-altered and
        # assigned variant allele specific cn between 0-2
        cn = rep(1, n)
        if (cn.rate > 0){
            num.cn.variants = ceiling(cn.rate*n)
            cn[sample(1:n, num.cn.variants, replace=F)] = 1 + runif(num.cn.variants,-1,1)
        }

        # draw reads for each variant
        for (j in 1:n){
            # set up underlying true var/ref reads in the sample
            reads = rep(0, total.reads)
            # calc num of ref reads given the cn of variant allele, and cn-neutral vaf
            total.num.ref.reads = floor((1-vaf)*total.reads/((cn[j]-1)*vaf+1))
            total.num.var.reads = total.reads - total.num.ref.reads
            reads[sample(1:total.reads, total.num.var.reads , replace=F)] = 1

            # sampling var/ref reads (this could also be done using existing generator
            # such as rbinom, however, let's just simulate via (slower) random sampling
            # as it mimics the sequencing process and introduces more variablility
            depth = depths[j]
            sim.reads = reads[sample(1:total.reads, depth, replace=F)]
            tmp = data.frame(sim.read=sim.reads)
            sim.var.cnt = sum(sim.reads)
            # scale # of variant reads according to cn.
            # this should result in the same var read counts as if it is drawn from
            # cn-scaled reads in the sample, but scaling here is more convenient
            # sim.var.cnt = round(sim.var.cnt*cn[j]/2)
            sim.ref.cnt = depth - sim.var.cnt
            # throw in sequencing errors by changing the nucleotide of some reads to
            # one of the three other nucleotides.
            num.err.reads = ceiling(err.rate*depth)
            err.reads = rep(F, depth)
            err.reads[sample(1:depth, num.err.reads, replace=F)] = T
            ref.err.reads.idx = sim.reads==0 & err.reads
            var.err.reads.idx = sim.reads==1 & err.reads
            num.ref.err.reads = sum(ref.err.reads.idx)
            num.var.err.reads = sum(var.err.reads.idx)
            # assuming ref base = 0, variant base = 1, other two bases = 2, 3.
            if (num.ref.err.reads > 0){
                sim.reads[ref.err.reads.idx] = sample(c(1,2,3), num.ref.err.reads, replace=T)
            }
            if (num.var.err.reads > 0){
                sim.reads[var.err.reads.idx] = sample(c(0,2,3), num.var.err.reads, replace=T)
            }
            tmp = cbind(tmp, sim.read.w.err=sim.reads)
            tmp$err = tmp$sim.read != tmp$sim.read.w.err
            # remove error reads that are not either var or ref reads
            sim.reads = sim.reads[sim.reads < 2]
            actual.depth = length(sim.reads)
            sim.var.cnt = sum(sim.reads)
            sim.ref.cnt = actual.depth - sim.var.cnt
            sim.vaf = sim.var.cnt/actual.depth*100
            sim.vafs = c(sim.vafs, sim.vaf)
            sim.ref.cnts = c(sim.ref.cnts, sim.ref.cnt)
            sim.var.cnts = c(sim.var.cnts, sim.var.cnt)
        }
        v[[paste0(s, '.vaf')]] = sim.vafs
        v[[paste0(s, '.ref.cnt')]] = sim.ref.cnts
        v[[paste0(s, '.var.cnt')]] = sim.var.cnts
    }
    return(v)
}

# determine mean/median of a vector of vaf
center.vaf <- function(x, method='mean'){
    if (method == 'mean'){
        return(mean(x))
    }else if (method == 'median'){
        return(median(x))
    }
}


# check if sum rule is strictly violated (ie. sum of children mean vafs > parent mean vaf)
# for a simulation
check.sum.violation <- function(v, x, vaf.cols){
    samples = vaf.cols
    clones = x$clone
    v = estimate.clone.vaf(v, cluster.col.name='clone', vaf.col.names=vaf.cols, method='mean')
    stats = NULL
    for (s in samples){
        # traverse from leaves to root
        #assume that parent clone is always numbered smaller than children clones
        cluster.vafs = v[[s]]
        violated.clones = NULL
        status = 'agree'
        for (i in length(clones):1){
            cl = clones[i]
            children = x$clone[x$parent == cl]
            if (length(children) > 0){#parental clones, addup all children
                if (cluster.vafs[cl] < sum(cluster.vafs[children])){
                    status = 'contr'
                    violated.clones = c(violated.clones, cl)
                }
            }
        }
        this.stat = data.frame(sample=s, status=status,
            violated.clones=paste(violated.clones, collapse=','), stringsAsFactors=F)
        if (is.null(stats)){stats = this.stat}else{stats = rbind(stats, this.stat)}
    }
    return(stats)
}

#' convert ground-truth clone ccf to variant vaf
#' @description Given a clonal evolution tree and the observed ccf of the clones
#' in individual samples, calculate vaf of the marker variants
#' in individual samples, to be used in simulation
#' @param x: data.frame with at least the following columns:
#' clone,parent,num.vars,sample1,sample2,...
#' additionally x can have *.ccf and *.vaf columns corresponding to
#' the samples, but won't be used.
#' @return the same data.frame with addtionally attached vaf columns
treeccf2vaf <- function(x){
    samples = setdiff(colnames(x), c('clone', 'parent', 'num.vars',
        grep('vaf|ccf', colnames(x), value=T)))
    clones = x$clone
    # order clones from leaves up to root
    clones = 1
    ord.clones = 1
    while (length(clones) > 0){
        cl = clones[1]; clones = clones[-1]
        children = x$clone[x$parent == cl]
        if (length(children) > 0){
            ord.clones = c(ord.clones, children)
            clones = c(clones, children)
        }
    }
    
    for (s in samples){
        # traverse from leaves to root
        # for convenience, assume that parent clone is always
        # numbered smaller than children clones 
        cluster.ccfs = x[[s]]
        for (i in length(ord.clones):1){
            cl = ord.clones[i]
            children = x$clone[x$parent == cl]
            if (length(children) > 0){#parental clones, addup all children
                cluster.ccfs[cl] = cluster.ccfs[cl] + sum(cluster.ccfs[children])
            }
        }
        x[[paste0(s, '.ccf')]] = cluster.ccfs
        x[[paste0(s, '.vaf')]] = cluster.ccfs/2
    }
    return(x)
}



#' plot boxplot of clone vaf across samples
plot.vaf <- function(x, out.pdf.file){
    # plot vafs
    pdf(out.pdf.file, width=4, height=7, useDingbats=FALSE, title='')
    pp = variant.box.plot(x,
            cluster.col.name='clone',
            cluster.axis.name='',
            show.cluster.size=F,
            cluster.size.text.color='blue',
            vaf.col.names = grep('vaf', colnames(x), value=T),
            variant.class.col.name=NULL,
            sample.title.size=20,
            vaf.limits=1,
            violin=F, violin.fill.color='gray', violin.alpha=0.2,
            box=F, jitter=T, jitter.shape=1,
            #jitter.color='#80b1d3',
            jitter.size=3,
            jitter.alpha=0.75,
            jitter.center.method='median',
            jitter.center.size=1,
            #jitter.center.color='#fb8072',
            jitter.center.color='#737373',
            jitter.center.display.value='none',
            highlight='cancer.gene',
            highlight.note.col.name='gene',
            highlight.color='red',
            highlight.note.size=2.5,
            highlight.shape=16,
            order.by.total.vaf=F,
    )
    dev.off()
}


#' Generate a random tree
#' @description Generate random tree rooted at clone 1
#' cpars = constrainted parents are the prechosen parents
#' for each node that has to be obey in the random tree
generate.random.tree <- function(x, cpars, samples){
    clones = x$clone
    pars = rep(NA, length(clones))
    pars[1] = -1
    assigned = c(1)
    num.clones = length(clones)
    unassigned = sample(setdiff(clones, assigned))
    cnt = 1
    while (cnt < num.clones){
        cl = unassigned[1]
        cat('randomize parent of clone', cl, '\n')
        if (is.na(cpars[[cl]][1])){
            par = sample(assigned, 1)
        }else{
            this.cpars = cpars[[cl]]
            this.cpars = intersect(assigned, this.cpars)
            if (length(this.cpars) > 0){
                # some constrainted parents already in tree, take them
                #par = cpars[[cl]][sample(length(cpars[[cl]]), 1)]
                par = this.cpars[sample(length(this.cpars), 1)]
                cat(cl, par, paste(cpars[[cl]],collapse=','),'\n')
            }else{
                cat('WARN: None of constraint parents are in tree (clone', cl, ')\n')
                print(this.cpars)
                # let assign parent for this clone (cl) later
                unassigned = c(unassigned[-1], cl)
                next
            }
        }
        pars[cl] = par
        assigned = c(assigned, cl)
        unassigned = unassigned[-1]
        cnt = cnt + 1
    }

    x$parent = pars[x$clone]
    #print(cbind(x[, 1:2], cons))

    return(x)
}


#' make a data frame of a tree to be used in clonevol tree plotting
make.clonevol.tree <- function(x, samples){
    x$lab = x$clone
    x$parent[x$clone == 1] = -1
    x$color = get.clonevol.colors(nrow(x))
    x$sample.with.nonzero.cell.frac.ci = ''
    x$samples.with.nonzero.cell.frac = ''
    x$sample.group = 'grp1'
    x$sample.group.color = 'black'
    x$excluded = F
    for (i in 1:nrow(x)){
        x$sample.with.nonzero.cell.frac[i] = paste(samples[x[i, samples] >0], collapse=',')
    }
    x$sample.with.nonzero.cell.frac.ci = x$sample.with.nonzero.cell.frac
    x$samples.with.nonzero.cell.frac = x$sample.with.nonzero.cell.frac
    x = convert.clone.to.branch(x, branch.lens = NULL, merged.tree.node.annotation='none')
    return(x)
}


#### BEGIN

args = commandArgs(trailingOnly=T)
if (length(args) != 10){stop(USAGE)}

subpop.file = args[1]
par.cons = ifelse(args[2] == 'T', T, F)
mean.depth = as.integer(args[3])
num.trees = as.integer(args[4])
num.replicates = as.integer(args[5])
error.rate = as.numeric(args[6])
cn.rate = as.numeric(args[7])
rand.seed = as.integer(args[8])
out.dir = args[9]
prefix = args[10]
dir.create(out.dir)

cat('subpopulations.truth.file:', subpop.file,
    '\nparental constraint:', par.cons, 
    '\nmean.depth:', mean.depth, 
    '\nnum.trees:', num.trees,
    '\nnum.replicates:', num.replicates,
    '\nerror.rate:', error.rate,
    '\ncn.rate:', cn.rate,
    '\nrand.seed:', rand.seed,
    '\nout.dir:', out.dir,
    '\nprefix:', prefix,
    '\n')

if (!is.na(rand.seed)){
    set.seed(rand.seed)
}


# read ground truth CCF
x = read.table(subpop.file, header=T, sep='\t', stringsAsFactors=F, na.strings='')
if (!par.cons){
    x$parent = NA
}
x$parent[1] = '-1,-1' # make sure 1st clone is root of tree

# build seeding constraints using parent col
cons = x
cpars = list()
for (cl in cons$clone){
    cpars = c(cpars, list(as.integer(unlist(strsplit(cons$parent[cons$clone==cl], ',')))))
}

# determine sample names from input column names
samples = setdiff(colnames(x) , c(grep('ccf|vaf|parent', colnames(x), value=T),
    'clone', 'num.vars'))

# generate num.trees random trees, each with num.replicates cases
for (tr in 1:num.trees){
    cat('**** tree', tr, '\n')
    tree.out.dir = paste0(out.dir, '/', prefix, '.', tr)
    dir.create(tree.out.dir)
    if (tr > 0){
        x = generate.random.tree(x, cpars, samples)
        x = treeccf2vaf(x)
    }

    # plot vaf
    plot.vaf(x, paste0(tree.out.dir, '/vaf.truth.pdf'))

    # plot the ground truth tree just generated
    t = make.clonevol.tree(x, samples)
    pdf(paste0(tree.out.dir, '/tree.truth.pdf'), width=6, height=6, useDingbats=F)
    plot.tree.clone.as.branch(t, node.label.size=1, node.text.size=0.75,
        branch.width=0.75, angle=20, tree.label=paste0(prefix, '.', tr))
    dev.off()
    write.table(x, file=paste0(tree.out.dir, '/clones.truth.w-vaf.tsv'),
        row.names=F, sep='\t', quote=F)
    write.table(x[, c('clone', 'parent', samples)],
        file=paste0(tree.out.dir, '/tree.truth.tsv'),
        row.names=F, col.names=T, sep='\t', quote=F)


    clones = x$clone
    num.vars = x$num.vars
    samples = setdiff(colnames(x) , c(grep('ccf|vaf|parent', colnames(x), value=T),
        'clone', 'num.vars'))


    stat.vs.truth = data.frame(tree=0, rep=0, sample='-', status='-',
        violated.clones='-', stringsAsFactors=F)
    for (r in 1:num.replicates){
        cat('rep = ', r, ':')
        v = NULL
        for (cl in clones){
            cat(' ', cl)
            vi = generate.variants(cl, x$num.vars[x$clone == cl], samples,
                true.vafs=unlist(x[x$clone == cl, paste0(samples, '.vaf')]),
                mean.depth=mean.depth, err.rate=error.rate, cn.rate=cn.rate)
            if (is.null(v)){
                v = vi
            }else{
                v = rbind(v, vi)
            }
        }

        # some faked normal and chromosomal position (some tools require these)
        v$normal.vaf = 0
        v$normal.ref.cnt = mean.depth
        v$normal.var.cnt = 0
        v = cbind(gene_name='None', v)
        v = cbind(type='SNP', v)
        v = cbind(tier='tier1', v)
        v = cbind(stop=seq(1, nrow(v))*10, v)
        v = cbind(start=seq(1, nrow(v))*10, v)
        v = cbind(chromosome_name=1, v)

        write.table(v, file=paste0(tree.out.dir, '/variants.', r, '.tsv'), sep='\t',
                    row.names=F, quote=F)

        cmp = check.sum.violation(v, x, paste0(samples, '.vaf'))
        cmp = cbind(tree=paste0(prefix, '.', tr), rep=r, cmp)
        cat('\n')
        if (is.null(stat.vs.truth)){stat.vs.truth = cmp}
            else{stat.vs.truth = rbind(stat.vs.truth, cmp)}
    }
    stat.vs.truth = stat.vs.truth[-1,]
    write.table(stat.vs.truth, file=paste0(tree.out.dir, '/sim.stat.tsv'),
        sep='\t', row.names=F, quote=F)
}


