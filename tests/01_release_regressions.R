library(qPCRtools)

read.example <- function(name) {
  read.table(system.file("examples", name, package = "qPCRtools"), header = TRUE)
}

cq <- data.frame(
  Position = paste0("A", 1:12),
  Cq = c(22, 20, 22.2, 20.1, 21.8, 19.9,
         21, 20, 21.1, 20.0, 20.9, 19.9)
)
design <- data.frame(
  Position = paste0("A", 1:12),
  Group = rep(c("CK", "CK", "CK", "T", "T", "T"), each = 2),
  Gene = rep(c("G1", "REF"), 6),
  BioRep = rep(1:3, each = 2, times = 2)
)
ddct <- CalExp2ddCt(
  cq, design,
  ref_gene = "REF",
  ref_group = "CK",
  stat_method = "t.test",
  remove_outliers = FALSE,
  fig_type = "box"
)
tb <- ddct[["table"]]
stopifnot(abs(mean(tb$mean.expre[tb$Treatment == "CK"]) - 1) < 0.05)
stopifnot(all(abs(tb$se.expre - tb$sd.expre / sqrt(tb$n)) < 1e-12))
stopifnot(all(c("RawPValue", "PValue", "PAdjustMethod") %in% names(tb)))
stopifnot(is.finite(tb$PValue[tb$Treatment == "T"]))
stopifnot(unique(tb$PAdjustMethod[tb$Treatment == "T"]) == "none")

text.layers <- vapply(
  ddct[["figure"]]$layers,
  function(layer) inherits(layer$geom, "GeomText"),
  logical(1)
)
stopifnot(sum(text.layers) == 1)

bar <- CalExp2ddCt(
  cq, design,
  ref_gene = "REF",
  ref_group = "CK",
  remove_outliers = FALSE,
  fig_type = "bar"
)
expected <- unique(bar[["figure"]]$data[, c("gene", "Treatment", "mean.expre")])
expected <- expected[order(expected$gene, expected$Treatment), ]
bar.data <- ggplot2::ggplot_build(bar[["figure"]])$data[[1]]
stopifnot(max(abs(sort(bar.data$y) - sort(expected$mean.expre))) < 1e-8)

anova.design <- expand.grid(
  BioRep = 1:3,
  Group = c("CK", "T1", "T2"),
  Gene = c("G1", "REF"),
  stringsAsFactors = FALSE
)
anova.cq <- data.frame(
  Position = paste0("B", seq_len(nrow(anova.design))),
  Cq = ifelse(
    anova.design$Gene == "REF",
    20,
    22 + 2 * as.numeric(factor(anova.design$Group))
  )
)
anova.design$Position <- anova.cq$Position
anova.res <- CalExp2ddCt(
  anova.cq, anova.design,
  ref_gene = "REF",
  ref_group = "CK",
  stat_method = "anova",
  remove_outliers = FALSE,
  fig_type = "box"
)
anova.tb <- anova.res[["table"]]
nonref <- anova.tb$Treatment != "CK"
stopifnot(all(c("RawPValue", "PValue", "PAdjustMethod") %in% names(anova.tb)))
stopifnot(all(is.finite(anova.tb$RawPValue[nonref])))
stopifnot(all(is.finite(anova.tb$PValue[nonref])))
stopifnot(all(anova.tb$PValue[nonref] >= anova.tb$RawPValue[nonref] - 1e-12))
stopifnot(unique(anova.tb$PAdjustMethod[nonref]) == "Tukey")
invisible(ggplot2::ggplot_build(anova.res[["figure"]]))

rq.cq <- read.example("cal.expre.rqpcr.cq.txt")
rq.design <- read.example("cal.expre.rqpcr.design.txt")
rq <- CalExpRqPCR(
  rq.cq, rq.design,
  ref_gene = "OsUBQ",
  ref_group = "CK",
  stat_method = "t.test"
)
stopifnot(nrow(rq[["table"]]) > 0)
stopifnot(any(is.finite(rq[["table"]][["Expression"]])))
stopifnot(all(c("RawPValue", "PValue", "PAdjustMethod") %in% names(rq[["table"]])))

sc.cq <- read.example("calsc.cq.txt")
sc.info <- read.example("calsc.info.txt")
curve <- CalCurve(
  sc.cq, sc.info,
  lowest_concen = 4, highest_concen = 4096,
  dilution = 4, by_mean = TRUE
)
deprecated <- FALSE
invisible(withCallingHandlers(
  ggplot2::ggplot_build(curve[["figure"]]),
  warning = function(w) {
    if (grepl("deprecated", conditionMessage(w), ignore.case = TRUE)) {
      deprecated <<- TRUE
    }
    invokeRestart("muffleWarning")
  }
))
stopifnot(!deprecated)

cat("OK: root package release regression tests passed\n")
