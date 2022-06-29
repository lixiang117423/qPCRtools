pkgname <- "qPCRtools"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
base::assign(".ExTimings", "qPCRtools-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('qPCRtools')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("CalCurve")
### * CalCurve

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: CalCurve
### Title: Standard Curve Calculation.
### Aliases: CalCurve

### ** Examples

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



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("CalCurve", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("CalExp2ddCt")
### * CalExp2ddCt

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: CalExp2ddCt
### Title: Calculate expression using standard curve.
### Aliases: CalExp2ddCt

### ** Examples

df1.path = system.file("examples", "ddct.cq.txt", package = "qPCRtools")
df2.path = system.file("examples", "ddct.design.txt", package = "qPCRtools")

cq.table = data.table::fread(df1.path)
design.table = data.table::fread(df2.path)

CalExp2ddCt(cq.table,
            design.table,
            ref.gene = "OsUBQ",
            ref.group = "CK",
            stat.method = "t.test",
            fig.type = "box",
            fig.ncol = NULL) -> res

res[["table"]]
res[["figure"]]




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("CalExp2ddCt", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("CalExpCurve")
### * CalExpCurve

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: CalExpCurve
### Title: Calculate expression using standard curve.
### Aliases: CalExpCurve

### ** Examples

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

res[["table"]]
res[["figure"]]




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("CalExpCurve", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("CalRTable")
### * CalRTable

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: CalRTable
### Title: Calculate volume.
### Aliases: CalRTable

### ** Examples

df.1.path <- system.file("examples", "crtv.data.txt", package = "qPCRtools")
df.2.path <- system.file("examples", "crtv.template.txt", package = "qPCRtools")
df.1 <- data.table::fread(df.1.path)
df.2 <- data.table::fread(df.2.path)
result <- CalRTable(data = df.1, template = df.2, RNA.weight = 2)
head(result)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("CalRTable", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
