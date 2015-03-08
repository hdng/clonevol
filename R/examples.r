# EXAMPLE CODE

## Variants list with VAF and cluster identity provided
run.AML31.bootstrap <- function(){
  # read AML31 variant + cluster
  aml31.variants = read.table('samples/AML31.variants.tsv',
                              header=T, sep='\t', stringsAsFactors=F, quote='')
  aml31.variants = aml31.variants[!is.na(aml31.variants$cluster),]
  aml31.variants = aml31.variants[, c('cluster',
                                      grep('.vaf', colnames(aml31.variants), value=T))]
  colnames(aml31.variants) = gsub('AML31.|.vaf', '', colnames(aml31.variants))

  variants = aml31.variants
  vaf.col.names = c('Tumor', 'Relapse')
  #clone.vafs = estimate.clone.vaf(variants, 'cluster', vaf.col.names)
  x = infer.clonal.models(variants=variants, vaf.col.names=vaf.col.names,
                          subclonal.test='bootstrap', num.boots=1000,
                          founding.cluster=1,
                          p.value.cutoff=0.25)
  plot.clonal.models(x$models,
                     out.dir='test-out/bootstrap-test-aml31',
                     matched=x$matched,
                     out.format='png', overwrite.output=T,
                     tree.node.shape='square',
                     scale.monoclonal.cell.frac=TRUE,
                     cell.frac.ci=T)

  for (s in vaf.col.names){
    draw.sample.clones.all(x$models[[s]], paste0('test-out/bootstrap-test-aml31/', s))
  }

}

## Fixed VAF, not test, no variant provided
# AML31 two samples (primary vs. relapse) deepSeq
run.AML31 <- function(){
  c = read.table('samples/AML31.clusters.tsv', header=T)
  x = infer.clonal.models(c)
  plot.clonal.models(x$models, out.dir='test-out/clonevol-AML31-output',
                     matched=x$matched, out.format='png',
                     overwrite.output=T)
}

# a makeup example of 3 samples
run.3samples.example <- function(){
  c = read.table('clusters.3samples.tsv', header=T)
  x = infer.clonal.models(c)
  plot.clonal.models(x$models, out.dir='clonevol-3samples-output',
                     matched=x$matched, out.format='pdf', overwrite.output=T)
}

# a makeup example of 1 sample
run.1sample.example <- function(){
  c = read.table('clusters.1sample.tsv', header=T)
  x = infer.clonal.models(c)
  plot.clonal.models(x$models, out.dir='clonevol-1sample-output',
                     matched=x$matched, out.format='png', overwrite.output=T)
}
