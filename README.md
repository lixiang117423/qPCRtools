# qPCRtools

# CalRTable
Calculate RNA and other reagent volume required for reverse transcription.
```{r}
df.1.path <- system.file("examples", "crtv.data.txt", package = "qPCRtools")
df.2.path <- system.file("examples", "template", package = "qPCRtools")
df.1 <- data.table::fread(df.1.path)
df.2 <- data.table::fread(df.2.path)
result <- CalRTable(data = df.1, template = df.2, RNA.weight = 2)
head(result)
```
# CalCurve
Standard Curve Calculation
```{r}
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
p[["table"]]
p[["figure"]]
```
# calculate expression using curve
```{r}
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
res[["table"]]
res[["figure"]]
```
