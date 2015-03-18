# Test model inference using list of variants + VAFs + bootstrap

library(testthat)
library(clonevol)

context('Test CRC12')

test_that('CRC12 4 samples!!!!',{
  variants = crc12.variants.new
  variants = variants[, !grepl('_XT1', colnames(variants))]
  variants = variants[variants$cluster != 'c9',]
  variants = variants[variants$cluster != 'c11',]
  #variants = variants[variants$cluster != 'c23',]
  out.dir = 'test-out/CRC12-new'
  vaf.col.names = grep('.VAF', colnames(variants), value=T, fixed=T)
  vaf.col.names = vaf.col.names[!grepl('PBMC', vaf.col.names)]
  clone.vafs = estimate.clone.vaf(variants, 'cluster', vaf.col.names)
  x = infer.clonal.models(variants=variants, vaf.col.names=vaf.col.names,
                          subclonal.test='bootstrap', num.boots=1000,
                          founding.cluster='c1', min.cluster.vaf=0.025,
                          p.value.cutoff=0.25, verbose=F)
  num.models = nrow(x$matched$index)
  expect_equal(16, num.models)
})


test_that('AML31 tumor/relapse paired samples!!!!',{
  variants = aml31.variants
  vaf.col.names = c('Tumor', 'Relapse')
  clone.vafs = estimate.clone.vaf(variants, 'cluster', vaf.col.names)
  x = infer.clonal.models(variants=variants, vaf.col.names=vaf.col.names,
                          subclonal.test='bootstrap', num.boots=1000,
                          founding.cluster=1,
                          p.value.cutoff=0.25, verbose=F)

  num.models = nrow(x$matched$index)
  expect_equal(5, num.models)
  #plot.clonal.models(x$models,
  #                   out.dir='test-out/bootstrap-test-aml31',
  #                   matched=x$matched,
  #                   out.format='png', overwrite.output=T,
  #                   tree.node.shape='square',
  #                   scale.monoclonal.cell.frac=TRUE,
  #                   cell.frac.ci=T)
  #for (s in vaf.col.names){
  #  draw.sample.clones.all(x$models[[s]], paste0('test-out/bootstrap-test-aml31/', s))
  #}
})


test_that('AML1 primary/relapse tumors',{
  variants = aml1.variants
  vaf.col.names = c('Tumor', 'Relapse')
  clone.vafs = estimate.clone.vaf(variants, 'cluster', vaf.col.names, vaf.in.percent=F)
  x = infer.clonal.models(variants=variants, vaf.col.names=vaf.col.names,
                          vaf.in.percent=F,
                          subclonal.test='bootstrap', num.boots=10000,
                          min.cluster.vaf=0.01,
                          founding.cluster=1, verbose=F)
  num.models = nrow(x$matched$index)
  expect_equal(1, num.models)
  #plot.clonal.models(x$models,
  #                   out.dir='test-out/bootstrap-test-aml1',
  #                   matched=x$matched,
  #                   out.format='png', overwrite.output=T,
  #                   tree.node.shape='square',
  #                   scale.monoclonal.cell.frac=TRUE,
  #                   width=10,
  #                   cell.frac.ci=T)
  #for (s in vaf.col.names){
  #  draw.sample.clones.all(x$models[[s]], paste0('test-out/bootstrap-test-aml1/', s))
  #}
})
