# ClonEvol
Inferring and visualizing clonal evolution in multi-sample cancer sequencing.

- What's new?: Major overhaul of code and features, including sphere of cells, annotated branch-based clonal evolution tree visualizations. Make sure you are running the latest version of ClonEvol in the following example.
- Join <a href="https://groups.google.com/forum/#!forum/clonevol">clonEvol mailing list</a> for annoucements, feature requests, Q/A, etc.

The following figure demonstrates the reanalysis of a relapse acute myeloid leukemia case (AML1) published in Ding et al., Nature (2012). Top panel shows the original published fishplot, and the bottom panel shows the model inferred and visualized by ClonEvol.

![](images/fig1-AML1.jpg)
*Fig. 1. ClonEvol reanalysis of AML1 (a) Original model, represented by a fishplot (b-f) Matching model predicted and visualized by ClonEvol.*

## How to install and use ClonEvol?

A tutorial is available <a href="https://github.com/hdng/clonevol/vignettes/clonevol.pdf">here</a>

## Known issues

Bell plots sometimes do not position nicely in plot.clonal.models function (eg. when there is a clone of extremely low cellular fraction together with complex clonal structure). Setting bell.curve.step=x where x is a small value (eg. x=0) or clone.shape='polygon' in plot.clonal.models function will fix it.

If you encounter this error: "Error: evaluation nested too deeply: infinite recursion / options(expressions=)?", increase recursive stack size by:

```{r}
options(expressions=100000)
```

## How to cite ClonEvol

Ha X. Dang, Brian S. White, Steven M. Foltz, Christopher A. Miller, Jingqin Luo, Ryan C. Fields, Christopher A. Maher. ClonEvol: clonal ordering and visualization in cancer sequencing (under review)

Ha X. Dang, Julie G. Grossman, Brian S. White, Steven M. Foltz, Christopher A. Miller, Jingqin Luo, Timothy J. Ley, Richard K. Wilson, Elaine R. Mardis, Ryan C. Fields, Christopher A. Maher. Clonal evolution inference and visualization in metastatic colorectal cancer. Late Breaking Research Track. Intelligent Systems for Molecular Biology (ISMB) 2016. Orlando, Florida, USA. Jul. 2016.

## Contact
Ha X. Dang @ haxdang (at) gmail (dot) com
