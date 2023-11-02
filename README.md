## Install

```{r
install.packages("qPCRtools")

# or
devtools::install_github("https://github.com/lixiang117423/qPCRtools/tree/main/CRAN/qPCRtools")
```

## Calculate volume for reverse transcription

The first step of qPCR is usually the preparation of cDNA. We need to calculate the column of RNA for reverse transcription to cDNA. So, if we have the concentration of RNA, we can use the function `CalRTable` to do that. The function have three patameters:

- `data`: The table of RNA concentration. The unit of concentration is ng/μl. The demo data can be found at [GitHub](https://github.com/lixiang117423/qPCRtools/blob/main/CRAN/qPCRtools/inst/examples/crtv.data.txt).
- `template`: The table of reagent for reverse transcription. The demo data can be found at [GitHub](https://github.com/lixiang117423/qPCRtools/blob/main/CRAN/qPCRtools/inst/examples/crtv.template.txt). The column `All` is the total volume for 1 μg RNA.
- `RNA.weight`: The mass of RNA. The unit is μg. The default value is 2.

```{r
suppressMessages(library(tidyverse))
library(qPCRtools)

df.1.path <- system.file("examples", "crtv.data.txt", package = "qPCRtools")
df.2.path <- system.file("examples", "crtv.template.txt", package = "qPCRtools")
df.1 <- data.table::fread(df.1.path)
df.2 <- data.table::fread(df.2.path)
result <- CalRTable(data = df.1, template = df.2, RNA.weight = 2)

result %>% 
  dplyr::slice(1:6) %>% 
  kableExtra::kable(format = "html") %>% 
  kableExtra::kable_styling("striped")
```

## Calculate standard curve

The function can calculate the standard curve. At the same time, it can get the amplification efficiency of primer(s). Based on the amplification efficiency, we can know which method can be used to calculate the expression level. The function has 6 parameters:

- `cq.table`: The table of Cq. It must contain at least two columns:One `Position` and `Cq`. The demo data can be found at [GitHub](https://github.com/lixiang117423/qPCRtools/blob/main/CRAN/qPCRtools/inst/examples/calsc.cq.txt).
- `concen.table`: The table of gene(s) and concentration. It must contain at least three columns: `Position`, `Gene` and `Conc`. The demo data can be found at [GitHub](https://github.com/lixiang117423/qPCRtools/blob/main/CRAN/qPCRtools/inst/examples/calsc.info.txt).
- `lowest.concen`: The lowest concentration used to calculate the standard curve.
- `highest.concen`: The highest concentration used to calculate the standard curve.
- `dilu`: The dilution factor of cDNA template. The default value is 4.
- `by`: Calculate the standard curve by average data or the full data. The default value is `mean`.

```{r
library(qPCRtools)

df.1.path <- system.file("examples", "calsc.cq.txt", package = "qPCRtools")
df.2.path <- system.file("examples", "calsc.info.txt", package = "qPCRtools")
df.1 <- data.table::fread(df.1.path)
df.2 <- data.table::fread(df.2.path)
CalCurve(
  cq.table = df.1,
  concen.table = df.2,
  lowest.concen = 4,
  highest.concen = 4096,
  dilu = 4,
  by = "mean"
) -> p

p[["table"]] %>% 
  dplyr::slice(1:6) %>% 
  kableExtra::kable(format = "html") %>% 
  kableExtra::kable_styling("striped")

p[["figure"]]

```

## Calculate expression using standard curve

After we calculated the standard curve, we can use the standard curve to calculate the expression level of genes. In `qPCRtools`, function `CalExpCurve` can get the expression using standard curve. There are several parameters in this function:

- `cq.table`: The table of Cq. It must contain at least two columns:One `Position` and `Cq`. The demo data can be found at [GitHub](https://github.com/lixiang117423/qPCRtools/blob/main/CRAN/qPCRtools/inst/examples/cal.exp.curve.cq.csv).
- `curve.table`: The table of standard curve calculated by `CalCurve`.
- `design.table`: The design information including three columns: `Position`,	`Treatment` and	`Gene`. The demo table can be found at [GitHub](https://github.com/lixiang117423/qPCRtools/blob/main/CRAN/qPCRtools/inst/examples/cal.exp.curve.design.txt).
- `correction`: Expression level is corrected or not with internal reference genes. The default value is `TRUE`.
- `ref.gene`: The name of reference gene.
- `stat.method`: The method used to calculate differential expression of genes. If we want to calculate the difference between target group and reference group, one of `t.test` or `wilcox.test` can be used. `anova` is for all groups. The default value is `t.test`.
- `ref.group`: The name of reference group. If `stat.method` is `t.test` or `wilcox.test`, the function need a `ref.group`.
- `fig.type`: The type of figure, `box` or `bar`. `box` represents `boxplot`. `bar` represents `barplot`. The default value is `box`.
- `fig.ncol`: The column of figure. The default value is `NULL`.

```{r}
df1.path = system.file("examples", "cal.exp.curve.cq.txt", package = "qPCRtools")
df2.path = system.file("examples", "cal.expre.curve.sdc.txt", package = "qPCRtools")
df3.path = system.file("examples", "cal.exp.curve.design.txt", package = "qPCRtools")

cq.table = data.table::fread(df1.path)
curve.table = data.table::fread(df2.path)
design.table = data.table::fread(df3.path)

CalExpCurve(
  cq.table,
  curve.table,
  design.table,
  correction = TRUE,
  ref.gene = "OsUBQ",
  stat.method = "t.test",
  ref.group = "CK",
  fig.type = "box",
  fig.ncol = NULL) -> res

res[["table"]] %>% 
  dplyr::slice(1:6) %>% 
  kableExtra::kable(format = "html") %>% 
  kableExtra::kable_styling("striped")
res[["figure"]]
```

## Calculate expression using 2-ΔΔCt

$2^{-{Δ}{Δ}{C_t }} $is a widely used method to calculate qPCR data[@livak2001analysis]. Our function `CalExp2ddCt` can do it. Seven parameters are required for this function:

- `cq.table`: The demo file can be found at [GitHub](https://github.com/lixiang117423/qPCRtools/blob/main/CRAN/qPCRtools/inst/examples/ddct.cq.txt).
- `design.table`: The demo data can be found at [GitHub](https://github.com/lixiang117423/qPCRtools/blob/main/CRAN/qPCRtools/inst/examples/ddct.design.txt).
  Other parameters are same as the function `CalExpCurve`.
- `ref.gene`: The name of reference gene.
- `ref.group`: The name of reference group. If `stat.method` is `t.test` or `wilcox.test`, the function need a `ref.group`.
- `stat.method`: The method used to calculate differential expression of genes. If we want to calculate the difference between target group and reference group, one of `t.test` or `wilcox.test` can be used. `anova` is for all groups. The default value is `t.test`.
- `fig.type`: The type of figure, `box` or `bar`. `box` represents `boxplot`. `bar` represents `barplot`. The default value is `box`.
- `fig.ncol`: The column of figure. The default value is `NULL`.

```{r
df1.path = system.file("examples", "ddct.cq.txt", package = "qPCRtools")
df2.path = system.file("examples", "ddct.design.txt", package = "qPCRtools")

cq.table = data.table::fread(df1.path)
design.table = data.table::fread(df2.path)

CalExp2ddCt(cq.table,
            design.table,
            ref.gene = "OsUBQ",
            ref.group = "CK",
            stat.method = "t.test",
            fig.type = "bar",
            fig.ncol = NULL) -> res

res[["table"]] %>% 
  dplyr::slice(1:6) %>% 
  kableExtra::kable(format = "html") %>% 
  kableExtra::kable_styling("striped")

res[["figure"]]
```

## Calculate expression using RqPCR

The method from [SATQPCR](http://satqpcr.sophia.inra.fr/cgi/home.cgi) can identify the most stable reference genes (REF) across biological replicates and technical replicates[@rancurel2019satqpcr]. Our package provides a function, `CalExpRqPCR`, to achieve it. In the `design.table`, `BioRep`,	`TechRep` and `Eff` are required. `BioRep` is the `biological replicates`. `TechRep` is the `technical replicates`. `Eff` is the amplification efficiency of genes. The `cq.table` can be found at [GitHub](https://github.com/lixiang117423/qPCRtools/blob/main/CRAN/qPCRtools/inst/examples/cal.expre.rqpcr.cq.txt) and the `design,table` can be found at [GitHub](https://github.com/lixiang117423/qPCRtools/blob/main/CRAN/qPCRtools/inst/examples/cal.expre.rqpcr.design.txt). If user want to give reference gene, `ref.gene` can be used (The default is `NULL`).

- `ref.group`: The name of reference group. If `stat.method` is `t.test` or `wilcox.test`, the function need a `ref.group`.
- `stat.method`: The method used to calculate differential expression of genes. If we want to calculate the difference between target group and reference group, one of `t.test` or `wilcox.test` can be used. `anova` is for all groups. The default value is `t.test`.
- `fig.type`: The type of figure, `box` or `bar`. `box` represents `boxplot`. `bar` represents `barplot`. The default value is `box`.
- `fig.ncol`: The column of figure. The default value is `NULL`.

```{r
df1.path <- system.file("examples", "cal.expre.rqpcr.cq.txt", package = "qPCRtools")
df2.path <- system.file("examples", "cal.expre.rqpcr.design.txt", package = "qPCRtools")

cq.table <- data.table::fread(df1.path, header = TRUE)
design.table <- data.table::fread(df2.path, header = TRUE)

CalExpRqPCR(cq.table,
            design.table,
            ref.gene = NULL,
            ref.group = "CK",
            stat.method = "t.test",
            fig.type = "bar",
            fig.ncol = NULL
            ) -> res

res[["table"]] %>% 
  dplyr::slice(1:6) %>% 
  kableExtra::kable(format = "html") %>% 
  kableExtra::kable_styling("striped")

res[["figure"]]
```

## References

If this package is used in your publication, please cite `qPCRtools` paper:

>[Li X, Wang Y, Li J, et al. qPCRtools: An R package for qPCR data processing and visualization[J]. Frontiers in Genetics, 2022, 13: 1002704.](https://www.frontiersin.org/articles/10.3389/fgene.2022.1002704/full)
