run.clonevol.from.variants <- function(variants, models, num.true.models,
                                       out.dir, case, vaf.col.names=NULL){
  if (is.null(vaf.col.names)){
    vaf.col.names = grep('.WGS_VAF', colnames(variants), value=T, fixed=T)
  }
  cluster.col.name = 'cluster'
  var = variants[, c(cluster.col.name, vaf.col.names)]
  colnames(var) = gsub('CRC\\d+_\\d+_|_\\d+', '', colnames(var))
  var = var[, !grepl('_XT|PBMC|normal', colnames(var))]
  if (case != 'AML31'){
    vaf.col.names = grep('VAF', colnames(var), value=T)
  }

  #CRC10
  if (case == 'CRC10'){
    var = var[var$cluster != 1, ]
  }
  if (case == 'CRC12' && 'tumor.WGS_VAF' %in% vaf.col.names){# prev data
    var[[cluster.col.name]][var[[cluster.col.name]] %in% c(5,6,7)] = 5
  }


  v = estimate.clone.vaf(var, cluster.col.name, vaf.col.names)
  #v[, vaf.col.names] = v[, vaf.col.names]/100

  #CRC12
  if (case == 'CRC12'){
    colnames(v) = gsub('met7.WGS_VAF', 'Li3.VAF', colnames(v))
    colnames(v) = gsub('met8.WGS_VAF', 'Li2.VAF', colnames(v))
    colnames(v) = gsub('met9.WGS_VAF', 'Li6.VAF', colnames(v))
    colnames(v) = gsub('tumor.WGS_VAF', 'C.VAF', colnames(v))


    v[3, 'Li2.VAF'] = v[1, 'Li2.VAF']
    v[4, 'Li2.VAF'] = v[1, 'Li2.VAF']
    v[2, 'Li6.VAF'] = v[1, 'Li6.VAF']
    v[2, 'Li3.VAF'] = v[1, 'Li3.VAF']


  }

  if (case == 'AML31'){
    v[3, 'Relapse'] = v[1, 'Relapse']
  }

  clones = v
  run.test(clones, models, num.true.models, out.dir=out.dir)

}

test.CRC <- function(){
  models = c('monoclonal', 'polyclonal')
  num.true.models = c('monoclonal'=0, 'polyclonal'=247)
  run.clonevol.from.variants(crc10.variants, models, num.true.models,
                             out.dir='CRC10-4samples-impure', 'CRC10')

  num.true.models = c('monoclonal'=1, 'polyclonal'=8)
  run.clonevol.from.variants(crc12.variants, models, num.true.models,
                             out.dir='CRC12-4samples-impure', 'CRC12')



  num.true.models = c('monoclonal'=5, 'polyclonal'=38)
  run.clonevol.from.variants(aml31.variants, models, num.true.models,
                             out.dir='AML31-2samples-impure', 'AML31',
                             vaf.col.names=c('Tumor', 'Relapse'))

  variants = crc12.variants
  variants[variants$cluster %in% c(5,6,7),]$cluster = 5
  variants[variants$cluster > 5,]$cluster = variants[variants$cluster > 5,]$cluster - 2
  vaf.col.names = grep('.WGS_VAF', colnames(variants), value=T, fixed=T)
  vaf.col.names = vaf.col.names[!grepl('normal', vaf.col.names)]
  highlight = c()
  boxPlots3(variants, 'cluster',
    vaf.col.names, length(vaf.col.names), horizontal=F,
    showClusterSize=F, sumName=NULL,
    outPlotPrefix=paste('test-out/CRC12-new/boxplot', sep=''),
    yMax=70, width=0, height=0, width1=0, height1=0,
    hscale=1.25, panel.border.linetype='solid',
    panel.border.linesize=2, medianCutOff4Coloring=0,
    q3CutOff4Coloring=0, box=T, violin=F, violinLineType='dashed',
    violinLineSize=0.5, sampleTitleSize=14, highLight=c(),
    sizeName=NULL, colorName='silent')
  num.true.models = c('monoclonal'=1, 'polyclonal'=8)
  vaf.col.names = grep('.WGS_VAF', colnames(variants), value=T, fixed=T)
  vaf.col.names = vaf.col.names[!grepl('_XT|PBMC|normal', vaf.col.names)]
  clone.vafs = estimate.clone.vaf(variants, 'cluster', vaf.col.names)
  adj.clone.vafs = adjust.clone.vaf(clone.vafs, variants, 'cluster')
  run.test(adj.clone.vafs, models, num.true.models,
           out.dir='CRC12-new',
           clone.width=7, sample.height=8, text.size=2)



  variants = aml31.variants
  vaf.col.names = c('Tumor', 'Relapse')
  highlight = c()
  dir.create('test-out/AML31-new')
  boxPlots3(variants, 'cluster', sumName=NULL,
            vaf.col.names, length(vaf.col.names), horizontal=F,
            showClusterSize=F,
            outPlotPrefix=paste('test-out/AML31-new/boxplot', sep=''),
            yMax=60, width=0, height=0, width1=0, height1=0,
            hscale=1.25, panel.border.linetype='solid',
            panel.border.linesize=2, medianCutOff4Coloring=0,
            q3CutOff4Coloring=0, box=T, violin=F, violinLineType='dashed',
            violinLineSize=0.5, sampleTitleSize=12, highLight=c(),
            sizeName=NULL, colorName='silent')
  num.true.models = c('monoclonal'=5, 'polyclonal'=38)
  clone.vafs = estimate.clone.vaf(variants, 'cluster', vaf.col.names)
  adj.clone.vafs = adjust.clone.vaf(clone.vafs, variants, 'cluster')
  run.test(adj.clone.vafs, models, num.true.models,
           out.dir='AML31-new', text.size=1)


  # use vaf from CMiller
  # this did not produce any models, need to change some VAF to zero
  # to make biological sense
  v = read.table('samples/AML31.vaf.cmiller-final.tsv', header=T)
  v[2, 'Relapse'] = 0
  v[4, 'Relapse'] = 0
  v[6, 'Relapse'] = 0
  run.test(v, models, num.true.models,
           out.dir='AML31-cmiller-final', text.size=1)



}
