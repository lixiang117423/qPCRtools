library(qPCRtools)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("usage: 01_run_version_case.R <version-label> <output-dir>")
}
label <- args[[1]]
out_dir <- args[[2]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

example <- function(name) {
  sep <- if (grepl("\\.csv$", name, ignore.case = TRUE)) "," else "\t"
  read.table(system.file("examples", name, package = "qPCRtools"),
             header = TRUE, check.names = FALSE, sep = sep)
}

write_case <- function(x, name) {
  utils::write.csv(x, file.path(out_dir, paste0(label, "_", name, ".csv")),
                   row.names = FALSE, fileEncoding = "UTF-8")
}

sc_cq <- example("calsc.cq.txt")
sc_info <- example("calsc.info.txt")
curve <- CalCurve(
  sc_cq, sc_info,
  lowest.concen = 4, highest.concen = 4096,
  dilution = 4, by.mean = TRUE
)
write_case(curve[["table"]], "calcurve")

group_cq <- example("cal.exp.curve.cq.txt")
group_design <- example("cal.exp.curve.design.txt")
bio_design <- example("ddct.design.txt")
bio_map <- unique(bio_design[, c("Position", "BioRep")])
group_design <- merge(group_design, bio_map, by = "Position", all.x = TRUE)
curve_table <- example("cal.expre.curve.sdc.txt")

curve_args <- list(
  cq.table = group_cq,
  curve.table = curve_table,
  design.table = group_design,
  correction = TRUE,
  ref.gene = "OsUBQ",
  stat.method = "t.test",
  ref.group = "CK",
  fig.type = "box"
)
if ("dilution" %in% names(formals(CalExpCurve))) {
  curve_args$dilution <- 4
}
curve_expr <- do.call(CalExpCurve, curve_args)
write_case(curve_expr[["table"]], "calexpcurve")

ddct_cq <- example("ddct.cq.txt")
ddct_design <- example("ddct.design.txt")
ddct <- CalExp2ddCt(
  ddct_cq, ddct_design,
  ref.gene = "OsUBQ",
  ref.group = "CK",
  stat.method = "t.test",
  remove.outliers = FALSE,
  fig.type = "box"
)
write_case(ddct[["table"]], "calexp2ddct")

rq_cq <- example("cal.expre.rqpcr.cq.txt")
rq_design <- example("cal.expre.rqpcr.design.txt")
rq <- CalExpRqPCR(
  rq_cq, rq_design,
  ref.gene = NULL,
  ref.group = "CK",
  stat.method = "t.test",
  fig.type = "box"
)
write_case(rq[["table"]], "calexprqpcr")

capture.output(sessionInfo(), file = file.path(out_dir, paste0(label, "_sessionInfo.txt")))
cat("wrote", label, "case to", out_dir, "\n")
