# ClonEvol
Inferring and visualizing clonal evolution in multi-sample cancer sequencing

##Installation instructions

###Requirements:
- R 2.15 or later

###Install ClonEvol
```
> install.packages("devtools")
> library(devtools)
> install_github(“hdng/clonevol”)
```

###Install dependencies:

```
> install.packages(“ggplot2”)
> install.packages(“igraph”)
```

###Running ClonEvol

ClonEvol infers clonal evolution models in single sample or multiple samples using the clusters of variants identified previously using other methods such as sciClone or PyClone.

####Prepare input file
An input file typically has the following columns (* indicated mandatory):

1. cluster*: the cluster identity of the variant (make sure do not name cluster as “-1”. This value is reserved for ClonEvol internal use.
2. sample1.VAF*: VAF of the variant in sample1
3. sample1.Depth: depth of the variant in sample1
4. sample2.VAF: VAF of the variant in sample2
5. sample2.Depth: depth of the variant in sample2
6. Additional variant annotation columns

Example input file:

| cluster  |  prim.vaf  |  met1.vaf  |  met2.vaf |
|----------|------------|------------|-----------|
| 1        |  51        |  44        |  52       |
| 1        |  45        |  56        |  47       |
| 1        |  55        |  50        |  49       |
| 2        |  31        |  47        |  0        |
| 2        |  28        |  38        |  0        |
| 2        |  31        |  45        |  0        |
| 2        |  30        |  47        |  0        |
| 2        |  31        |  53        |  0        |
| 2        |  38        |  48        |  0        |

| ….

####Run ClonEvol

You can read your data into a data frame (eg. using read.table). Here let's use AML1 data (Ding et al., 2012) included in ClonEvol.

**Load AML1 data**
```{r}
library(clonevol)
data(aml1)
vaf.col.names = grep('.vaf', colnames(aml1), value=TRUE)
```

**Infer clonal evolution models**
```{r}
x = infer.clonal.models(variants=aml1,
            cluster.col.name="cluster",
            vaf.col.names=vaf.col.names,
            sample.names=c("primary", "relapse"),
            subclonal.test="bootstrap",
            subclonal.test.model="non-parametric",
            cluster.center="mean",
            num.boots=1000,
            founding.cluster=1,
            p.value.cutoff=0.01,
            alpha=0.1,
            random.seed=63108)
```

**Plot clonal evolution models**
```
> plot.clonal.models(x$models,
                       out.dir=”output”,
                       matched=x$matched,
                       variants=v,
                       box.plot=T,
                       out.format="pdf",
                       overwrite.output=T,
                       scale.monoclonal.cell.frac=T,
                       cell.frac.ci=T,
                       tree.node.shape="circle",
                       tree.node.size=40,
                       tree.node.text.size=0.65,
                       width=7, height=10)
```

**Plot box/violin/jitter of VAFs**
```
> variant.box.plot(aml1, vaf.col.names = vaf.col.names,
                          variant.class.col.name=NULL,
                          cluster.axis.name="",
                          vaf.limits=70,
                          violin=F,
                          box=F,
                          order.by.total.vaf=F,
                          jitter=T,
                          jitter.center.method=jitter.center.method,
                          jitter.center.size=0.5,
                          jitter.center.color='darkgray,
                          jitter.shape=1,
                          jitter.color=get.clonevol.colors(num.clusters),
                          jitter.size=2,
                          jitter.alpha=1)

```



