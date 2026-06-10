# qPCRtools

An R package for qPCR data processing and visualization.

## Installation

```r
install.packages("qPCRtools")
```

Development version:

```r
# install.packages("devtools")
devtools::install_github("lixiang117423/qPCRtools")
```

## Calculate volume for reverse transcription

Calculate the volume of RNA needed for reverse transcription to cDNA.

```r
library(qPCRtools)

df.1.path <- system.file("examples", "crtv.data.txt", package = "qPCRtools")
df.2.path <- system.file("examples", "crtv.template.txt", package = "qPCRtools")
df.1 <- read.table(df.1.path, sep = "\t", header = TRUE)
df.2 <- read.table(df.2.path, sep = "\t", header = TRUE)
result <- CalRTable(data = df.1, template = df.2, rna_weight = 2)
head(result)
```

## Calculate standard curve

Calculate the standard curve and obtain the amplification efficiency of primer(s).

```r
df.1.path <- system.file("examples", "calsc.cq.txt", package = "qPCRtools")
df.2.path <- system.file("examples", "calsc.info.txt", package = "qPCRtools")
df.1 <- read.table(df.1.path, header = TRUE)
df.2 <- read.table(df.2.path, header = TRUE)
res <- CalCurve(
  cq_table = df.1,
  concen_table = df.2,
  lowest_concen = 4,
  highest_concen = 4096,
  dilution = 4,
  by_mean = TRUE
)

res[["table"]]
res[["figure"]]
```

## Calculate expression using standard curve

```r
df1.path <- system.file("examples", "cal.exp.curve.cq.txt", package = "qPCRtools")
df2.path <- system.file("examples", "cal.expre.curve.sdc.txt", package = "qPCRtools")
df3.path <- system.file("examples", "cal.exp.curve.design.txt", package = "qPCRtools")

cq_table <- read.table(df1.path, header = TRUE)
curve_table <- read.table(df2.path, sep = "\t", header = TRUE)
design_table <- read.table(df3.path, header = TRUE)

res <- CalExpCurve(
  cq_table,
  curve_table,
  design_table,
  correction = TRUE,
  ref_gene = "OsUBQ",
  stat_method = "t.test",
  ref_group = "CK",
  fig_type = "box",
  fig_ncol = NULL
)

res[["table"]]
res[["figure"]]
```

## Calculate expression using 2-dCt

```r
df1.path <- system.file("examples", "dct.cq.txt", package = "qPCRtools")
df2.path <- system.file("examples", "dct.design.txt", package = "qPCRtools")
cq_table <- read.table(df1.path, sep = ",", header = TRUE)
design_table <- read.table(df2.path, sep = ",", header = TRUE)

res <- CalExp2dCt(cq_table, design_table, ref_gene = "Actin")
head(res)
```

## Calculate expression using 2-ddCt

```r
df1.path <- system.file("examples", "ddct.cq.txt", package = "qPCRtools")
df2.path <- system.file("examples", "ddct.design.txt", package = "qPCRtools")
cq_table <- read.table(df1.path, header = TRUE)
design_table <- read.table(df2.path, header = TRUE)

res <- CalExp2ddCt(
  cq_table,
  design_table,
  ref_gene = "OsUBQ",
  ref_group = "CK",
  stat_method = "t.test",
  remove_outliers = TRUE,
  fig_type = "box",
  fig_ncol = NULL
)

res[["table"]]
res[["figure"]]
```

## Calculate expression using RqPCR

The method from [SATQPCR](http://satqpcr.sophia.inra.fr/cgi/home.cgi) can identify the most stable reference genes across biological and technical replicates.

```r
df1.path <- system.file("examples", "cal.expre.rqpcr.cq.txt", package = "qPCRtools")
df2.path <- system.file("examples", "cal.expre.rqpcr.design.txt", package = "qPCRtools")
cq_table <- read.table(df1.path, header = TRUE)
design_table <- read.table(df2.path, header = TRUE)

res <- CalExpRqPCR(
  cq_table,
  design_table,
  ref_gene = NULL,
  ref_group = "CK",
  stat_method = "t.test",
  fig_type = "box",
  fig_ncol = NULL
)

res[["table"]]
res[["figure"]]
```

## Reference

If this package is used in your publication, please cite:

> [Li X, Wang Y, Li J, et al. qPCRtools: An R package for qPCR data processing and visualization[J]. Frontiers in Genetics, 2022, 13: 1002704.](https://www.frontiersin.org/articles/10.3389/fgene.2022.1002704/full)
