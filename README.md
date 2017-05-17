# ClonEvol
Inferring and visualizing clonal evolution in multi-sample cancer sequencing.

ClonEvol is currently under some major improvement to allow simulatenous analysis and intuitive visualization of a great number of samples across many types of tissues/organs. Join <a href="https://groups.google.com/forum/#!forum/clonevol">clonEvol mailing list</a> for annoucements, feature requests, Q/A, etc. The following figure is a reanalysis of a relapse AML case published in Griffith et al, Cell Systems 2015. Panels a and g are from the original publication. The rest of the figures are direct output of ClonEvol.

Technical documentation is under construction.

<img src="images/AML31-small.png" width="800">

## Installation

### Requirements:
- R 3.0.2 or later

### Install ClonEvol
```{r}
install.packages("devtools")
library(devtools)
install_github("clonevol", "hdng") or install_github("hdng/clonevol") if the former does not work
```

### Install dependencies

```{r}
install.packages("ggplot2")
install.packages("igraph")
```

## Run ClonEvol

ClonEvol infers clonal evolution models in single sample or multiple samples using the clusters of variants identified previously using other methods such as sciClone or PyClone. Variant clusters must be biologically interpretable for ClonEvol to be able to infer some models. Most of the time you will find yourself iteratively refining the clustering of the variants and running ClonEvol, until some reasonable models are found.

### Prepare input file
An input file typically has the following columns (* indicates mandatory):

1. cluster*: the cluster identity of the variant (make sure do not name cluster as “-1”. This value is reserved for ClonEvol internal use.)
2. sample1.VAF*: VAF of the variant in sample1
3. sample1.Depth: depth of the variant in sample1
4. sample2.VAF: VAF of the variant in sample2
5. sample2.Depth: depth of the variant in sample2
6. Additinal samples' VAF and depth columns
7. Additional variant annotation columns

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
...

### Run ClonEvol

You can read your data into a data frame (eg. using read.table). Here let's use AML1 data (Ding et al., 2012) included in ClonEvol.

**Load AML1 data**
```{r}
library(clonevol)
data(aml1)
x = aml1
vaf.col.names <- grep('.vaf', colnames(x), value=T)
sample.names = gsub('.vaf', '', vaf.col.names)
x[, sample.names] = x[, vaf.col.names]
vaf.col.names = sample.names
sample.groups = c('P', 'R');
names(sample.groups) = vaf.col.names
x = x[order(x$cluster),]

```

***Set up the colors***
```{r}
colors = c('#999793', '#8d4891', '#f8e356', '#fe9536', '#d7352e')
#colors = get.clonevol.colors(length(unique(v$cluster)))

```

***Prepare output directory***
```{r}
output.dir = 'output';
dir.create(out.dir)
```

**Visualize the clustering (and clean-up as needed, before running ClonEvol)**
ClonEvol takes clustering of variants and perform clonal ordering to infer the trees. Although it can tolerate errors in clustering, it is important to have the best clustering results possible to feed to ClonEvol. Plot them and see if they are reasonably good. If not, recluster and/or clean up. The following code will plot the clustering results for you to investigate. This plot is very powerful as it can visualize lots of samples and clusters at once. ClonEvol calls this the "boxplot", as the very first version only plot the box plots, but it now can plot jitter, box, and violin plots to allow close investigation of the clustering.

```{r}
# box plot
pdf(paste0(out.dir, '/box.pdf'), width=3, height=5, useDingbats=FALSE, title='')
pp = variant.box.plot(x,
    cluster.col.name = 'cluster',
    show.cluster.size = FALSE,
    cluster.size.text.color = 'blue',
    vaf.col.names = vaf.col.names,
    vaf.limits=70,
    sample.title.size=20,
    violin=F,
    box=F, jitter=T, jitter.shape=1,
    jitter.color=colors,
    jitter.size=3,
    jitter.alpha=1,
    jitter.center.method='median',
    jitter.center.size=1,
    jitter.center.color='darkgray',
    jitter.center.display.value='none',
    highlight='is.driver',
    highlight.note.col.name='gene',
    highlight.note.size=2,
    highlight.shape=16,
    order.by.total.vaf=F
)
dev.off()

```


**Infer clonal evolution models**
At this step, we assume that you already have reasonable good clustering results. Let's tell ClonEvol to perfrom clonal ordering and construct the consensus trees.

```{r}
y = infer.clonal.models(variants = x,
        vaf.col.names = vaf.col.names,
        sample.groups = sample.groups,
        subclonal.test = 'bootstrap',
        subclonal.test.model = 'non-parametric',
        num.boots = 1000,
        founding.cluster = '1',
        cluster.center='mean',
        ignore.clusters=NULL,
        clone.colors=colors,
        min.cluster.vaf=0.01,
        p.value.cutoff=0.05,
        alpha=0.05)

```

** Mapping driver events onto the trees
If the previous step succeeds (congrats!) and predicts some models, we can map some driver events onto the tree
```{r}
y = transfer.events.to.consensus.trees(y,
    x[x$is.driver,],
    cluster.col.name='cluster',
    event.col.name='gene')

and create node-based trees to branch-based trees as follows.
```{r} 
y = convert.consensus.tree.clone.to.branch(y, branch.scale='sqrt')
```

**Plot clonal evolution models**
Now it is exciting time, visualzing the clonal evolution models.
```{r}

```

Output should look like this:

<img src="images/model.png" width="600">


**Plot box/violin/jitter of VAFs with cancer gene variants highlighted**
```{r}
num.clusters <- length(unique(aml1$cluster))
pdf("variants.jitter.pdf", width=5, height=5, useDingbats=FALSE)
pp = variant.box.plot(aml1,
                 vaf.col.names=vaf.col.names,
                 variant.class.col.name=NULL,
                 cluster.axis.name="",
                 vaf.limits=70,
                 violin=FALSE,
                 box=FALSE,
                 order.by.total.vaf=FALSE,
                 jitter=TRUE,
                 jitter.center.method="mean",
                 jitter.center.size=0.5,
                 jitter.center.color="darkgray",
                 jitter.shape=1,
                 jitter.color=get.clonevol.colors(num.clusters),
                 jitter.size=2,
                 jitter.alpha=1,
                 highlight="is.cancer.gene",
                 highlight.note.col.name="gene",
                 highlight.shape=19,
                 display.plot=TRUE)
dev.off()
```

Output figure should look like this:

<img src="images/variants.jitter.png" width="400">


**Plot pairwise VAFs across samples**
```{r}
plot.pairwise(aml1, col.names=vaf.col.names,
                  out.prefix="variants.pairwise.plot",
                  colors=get.clonevol.colors(num.clusters))
```

**Plot mean/median of clusters across samples (cluster flow)**
```{r}
plot.cluster.flow(aml1, vaf.col.names=vaf.col.names,
                      sample.names=c("Primary", "Relapse"),
                      out.file="flow.pdf",
                      colors=get.clonevol.colors(num.clusters))

```

## Known issues
Bell plots sometimes do not position nicely in plot.clonal.models function (eg. when there is a clone of extremely low cellular fraction together with complex clonal structure). Setting bell.curve.step=x where x is a small value (eg. x=0) or clone.shape="polygon" in plot.clonal.models function will fix it.

If you encounter this error: "Error: evaluation nested too deeply: infinite recursion / options(expressions=)?", increase recursive stack size by:

```{r}
options(expressions=10000)
```

## How to cite clonEvol

Ha X. Dang, Brian S. White, Steven M. Foltz, Christopher A. Miller, Jingqin Luo, Ryan C. Fields, Christopher A. Maher. ClonEvol: clonal ordering and visualization in cancer sequencing (under review)

Ha X. Dang, Julie G. Grossman, Brian S. White, Steven M. Foltz, Christopher A. Miller, Jingqin Luo, Timothy J. Ley, Richard K. Wilson, Elaine R. Mardis, Ryan C. Fields, Christopher A. Maher. Clonal evolution inference and visualization in metastatic colorectal cancer. Late Breaking Research Track. Intelligent Systems for Molecular Biology (ISMB) 2016. Orlando, Florida, USA. Jul. 2016.

## Contact
Ha X. Dang @ haxdang (at) gmail (dot) com
