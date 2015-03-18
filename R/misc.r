
# generic function to call clonevol
run.clonevol <- function(c, out.dir){
  x = infer.clonal.models(c)
  plot.clonal.models(x$models, out.dir=out.dir,
                     matched=x$matched, out.format='png', overwrite.output=T)
}

run.AML31.2 <- function(){
  run.clonevol(aml31, 'clonevol-out-AML31')
}


# run.AML31()
# run.3samples.example()
# run.1sample.example()


# Prepare some example data
include.data.in.package <- function(){
  aml31 = read.table('samples/AML31.clusters.tsv', header=T)
  devtools::use_data(aml31, overwrite=T)
  threeSamples = read.table('samples/clusters.3samples.tsv', header=T)
  devtools::use_data(threeSamples, overwrite=T)

  crc10.variants = read.table('samples/CRC10.output.clusters.anno.processed.tsv',
                              header=T, sep='\t', stringsAsFactors=F, quote='')
  devtools::use_data(crc10.variants, overwrite=T)

  crc12.variants = read.table('samples/CRC12.output.clusters.anno.processed.prev.tsv',
                              header=T, sep='\t', stringsAsFactors=F, quote='')
  devtools::use_data(crc12.variants, overwrite=T)

  crc12.variants.new = read.table('samples/CRC12.new.tsv',
                              header=T, sep='\t', stringsAsFactors=F, quote='')
  colnames(crc12.variants.new) = gsub('CRC12_322_', '', colnames(crc12.variants.new))
  devtools::use_data(crc12.variants.new, overwrite=T)


  aml31.variants = read.table('samples/AML31.variants.tsv',
                              header=T, sep='\t', stringsAsFactors=F, quote='')
  aml31.variants = aml31.variants[!is.na(aml31.variants$cluster),]
  aml31.variants = aml31.variants[, c('cluster',
        grep('.vaf', colnames(aml31.variants), value=T))]
  colnames(aml31.variants) = gsub('AML31.|.vaf', '', colnames(aml31.variants))
  devtools::use_data(aml31.variants, overwrite=T)

  aml1.variants = read.table('samples/AML1.variants.tsv',
                              header=T, sep='\t', stringsAsFactors=F, quote='')
  aml1.variants[is.na(aml1.variants$Relapse),]$Relapse = 0
  devtools::use_data(aml1.variants, overwrite=T)

  msclc3151.ccf = read.table('samples/MSCLC3151.CCF.tsv',
                             header=T, sep='\t', stringsAsFactors=F, quote='')
  devtools::use_data(msclc3151.ccf, overwrite=T)

  msclc3151.ccf.touched = read.table('samples/MSCLC3151.CCF.touched.tsv',
                             header=T, sep='\t', stringsAsFactors=F, quote='')
  devtools::use_data(msclc3151.ccf.touched, overwrite=T)


  msclc3588.ccf = read.table('samples/MSCLC3588.CCF.tsv',
                             header=T, sep='\t', stringsAsFactors=F, quote='')
  devtools::use_data(msclc3588.ccf, overwrite=T)

  msclc984.ccf = read.table('samples/MSCLC984.CCF.tsv',
                             header=T, sep='\t', stringsAsFactors=F, quote='')
  devtools::use_data(msclc984.ccf, overwrite=T)

  crc8.variants = read.table('samples/CRC8.tsv',
                            header=T, sep='\t', stringsAsFactors=F, quote='')
  colnames(crc8.variants) = gsub('CRC8_237_', '', colnames(crc8.variants))
  colnames(crc8.variants) = gsub('_\\d+\\.', '.', colnames(crc8.variants))
  devtools::use_data(crc8.variants, overwrite=T)


}


# TEST CASES


# run a test on one clonality analysis output
run.test <- function(clones, models, num.true.models, out.dir=NULL,
                     clone.width=NULL, sample.height=NULL,
                     text.size=1){
  test.out.dir = 'test-out'
  if (!is.null(out.dir)){
    if (!file.exists(test.out.dir)){dir.create(test.out.dir)}
    this.out.dir = file.path(test.out.dir, out.dir)
    if (!file.exists(this.out.dir)){dir.create(this.out.dir)}
  }
  num.found.models = rep(NA, length(models))
  names(num.found.models) = models
  for (model in models){
    #cat('Testing ', model, ' model...\n')
    x = infer.clonal.models(c=clones, model=model, verbose=F)
    num.found.model = ifelse(is.null(x$matched), 0, nrow(x$matched$index))
    num.found.models[model] = num.found.model
    #if (!is.null(x$matched) && nrow(x$matched$index) == num.true.models[model]){
    #  cat('Test OK\n')
    #}else{
    #  message(paste0('Test FAILED. Should find ', num.true.models[model], ' models'))
    #}
    if (!is.null(out.dir)){
      plot.clonal.models(x$models, out.dir=file.path(this.out.dir, model),
                         matched=x$matched,
                         width=clone.width, height=sample.height,
                         text.size=text.size,
                         out.format='png', overwrite.output=T,
                         scale.monoclonal.cell.frac=FALSE)
      plot.clonal.models(x$models,
                         out.dir=file.path(this.out.dir, paste0(model,'-scaled')),
                         matched=x$matched,
                         width=clone.width, height=sample.height,
                         text.size=text.size,
                         out.format='png', overwrite.output=T,
                         scale.monoclonal.cell.frac=TRUE)
    }
  }
  return(num.found.models)
}

run.all.test <- function(){
  test.data = list()

  models = c('monoclonal', 'polyclonal')

  cat('*** Test 3a. 1 sample with 100% purity, two clones of same 0.5 VAF,
      expect to have 2 models in mono and 6 in polyclonal\n')
  clones = data.frame(cluster=c(1,2,3,4,5,6), primary=c(0.5, 0, 0.5, 0, 0.3, 0))
  num.true.models = c('monoclonal'=2, 'polyclonal'=2)
  run.test(clones, models, num.true.models, out.dir='1sample-2clones-sameVAF50-pu100')

  cat('*** Test 0. 1 sample with 80% purity, 1 clone, expect 1 model
      in both mono and polyclonal\n')
  clones = data.frame(cluster=c(1), primary=c(0.4))
  num.true.models = c('monoclonal'=1, 'polyclonal'=1)
  run.test(clones, models, num.true.models, out.dir='1sample-1clone-pu80')

  cat('*** Test 1. 1 sample with 100% purity, expect to have 2 models
      in both mono and polyclonal\n')
  clones = data.frame(cluster=c(2,3,4), primary=c(0.5,0.1,0.05))
  num.true.models = c('monoclonal'=2, 'polyclonal'=2)
  run.test(clones, models, num.true.models, out.dir='1sample-pu100')

  cat('*** Test 2. 1 sample with 70% purity, expect to have 2 models
      in mono and 6 in polyclonal\n')
  clones = data.frame(cluster=c(1,2,3), primary=c(0.35,0.1,0.05))
  num.true.models = c('monoclonal'=2, 'polyclonal'=6)
  run.test(clones, models, num.true.models, out.dir='1sample-pu70')

  cat('*** Test 3. 1 sample with 100% purity, two clones of same VAF,
      expect to have 3 models in mono and 6 in polyclonal\n')
  clones = data.frame(cluster=c(1,2,3), primary=c(0.5,0.2,0.2))
  num.true.models = c('monoclonal'=3, 'polyclonal'=3)
  run.test(clones, models, num.true.models, out.dir='1sample-2clones-sameVAF-pu100')


  cat('*** Test 4. 3 samples with 100% purity
      expect to have 1 model in mono and 1 in polyclonal\n')
  clones = threeSamples
  num.true.models = c('monoclonal'=1, 'polyclonal'=1)
  run.test(clones, models, num.true.models, out.dir='3sample-pu100')

  cat('*** Test 5. AML31 primary/relapse samples with VAF scaled to 0.5 (100% purity),
      expect 5 models in mono and 5 in polyclonal\n')
  #clones = data.frame(cluster=c(1,2,3), primary=c(0.5,0.2,0.2))
  clones = aml31
  num.true.models = c('monoclonal'=5, 'polyclonal'=5)
  run.test(clones, models, num.true.models, out.dir='AML31-pu100')

  cat('*** Test 6. AML31 primary/relapse samples with 90% purity in primary,
      and 40% purity in relapse. Expect 5 models in mono and 37 in polyclonal\n')
  #clones = data.frame(cluster=c(1,2,3), primary=c(0.5,0.2,0.2))
  clones = aml31
  clones$primary = clones$primary*0.9
  clones$relapse = clones$relapse*0.4
  num.true.models = c('monoclonal'=5, 'polyclonal'=37)
  run.test(clones, models, num.true.models, out.dir='AML31-pu-90-40')

}
models = c('monoclonal', 'polyclonal')

