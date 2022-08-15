## ----eval=FALSE---------------------------------------------------------------
#  install.packages("qPCRtools")

## ----echo=TRUE----------------------------------------------------------------
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

## ----echo=TRUE----------------------------------------------------------------
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


## -----------------------------------------------------------------------------
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

## ----echo=TRUE----------------------------------------------------------------
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

## ----echo=TRUE----------------------------------------------------------------
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

