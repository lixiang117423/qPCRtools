library(qPCRtools)

design <- rbind(
  data.frame(
    Group = "CK",
    BioRep = 1:3,
    Gene = rep(c("G1", "REF"), each = 3),
    stringsAsFactors = FALSE
  ),
  data.frame(
    Group = "T",
    BioRep = 1:4,
    Gene = rep(c("G1", "REF"), each = 4),
    stringsAsFactors = FALSE
  )
)
design$Position <- paste0("A", seq_len(nrow(design)))

cq <- data.frame(Position = design$Position, Cq = 20)
cq$Cq[design$Gene == "G1" & design$Group == "CK"] <- c(22, 23, 21)
cq$Cq[design$Gene == "G1" & design$Group == "T"] <- c(25, 24, 23, 10)

res <- CalExp2ddCt(
  cq,
  design,
  ref_gene = "REF",
  ref_group = "CK",
  stat_method = "t.test",
  remove_outliers = TRUE,
  fig_type = "box"
)

table <- res[["table"]]
stopifnot(nrow(table) == 6)
stopifnot(all(table$n == 3))
stopifnot(setequal(table$Treatment, c("CK", "T")))
stopifnot(is.finite(table$PValue[table$Treatment == "T"]))
stopifnot(nrow(res[["figure"]]$data) == 6)

cat("OK: unequal BioRep counts and IQR outlier removal passed\n")
