# EXAMPLE CODE

## Variants list with VAF and cluster identity provided
run.AML31.bootstrap <- function(){
  # read AML31 variant + cluster
  aml31.variants = read.table('samples/AML31.variants.tsv',
                              header=TRUE, sep='\t',
                              stringsAsFactors=FALSE, quote='')
  aml31.variants = aml31.variants[!is.na(aml31.variants$cluster),]
  aml31.variants = aml31.variants[, c('cluster',
                                      grep('.vaf', colnames(aml31.variants),
                                           value=TRUE))]
  colnames(aml31.variants) = gsub('AML31.|.vaf', '', colnames(aml31.variants))

  variants = aml31.variants
  vaf.col.names = c('Tumor', 'Relapse')
  #clone.vafs = estimate.clone.vaf(variants, 'cluster', vaf.col.names)
  model = 'normal-truncated'
  x = infer.clonal.models(variants=variants, vaf.col.names=vaf.col.names,
                          subclonal.test='bootstrap',
                          subclonal.test.model=model,
                          num.boots=1000,
                          founding.cluster=1,
                          p.value.cutoff=0.1)
  plot.clonal.models(models=x$models,
                     out.dir=paste0('test-out/bootstrap-test-aml31-', model),
                     matched=x$matched,
                     box.plot=FALSE,
                     variants=variants,
                     out.format='png', overwrite.output=TRUE,
                     tree.node.shape='square',
                     scale.monoclonal.cell.frac=TRUE,
                     cell.frac.ci=TRUE)

  for (s in vaf.col.names){
    draw.sample.clones.all(x$models[[s]], paste0('test-out/bootstrap-test-aml31/', s))
  }

}

## Fixed VAF, not test, no variant provided
# AML31 two samples (primary vs. relapse) deepSeq
run.AML31 <- function(){
  c = read.table('samples/AML31.clusters.tsv', header=TRUE)
  x = infer.clonal.models(c)
  plot.clonal.models(x$models, out.dir='test-out/clonevol-AML31-output',
                     matched=x$matched, out.format='png',
                     overwrite.output=TRUE)
}

# a makeup example of 3 samples
run.3samples.example <- function(){
  c = read.table('clusters.3samples.tsv', header=TRUE)
  x = infer.clonal.models(c)
  plot.clonal.models(x$models, out.dir='clonevol-3samples-output',
                     matched=x$matched, out.format='pdf', overwrite.output=TRUE)
}

# a makeup example of 1 sample
run.1sample.example <- function(){
  c = read.table('clusters.1sample.tsv', header=TRUE)
  x = infer.clonal.models(c)
  plot.clonal.models(x$models, out.dir='clonevol-1sample-output',
                     matched=x$matched, out.format='png', overwrite.output=TRUE)
}
