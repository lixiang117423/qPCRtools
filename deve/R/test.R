rm(list = ls())


df1 <- data.table::fread("inst/examples/calsc.cq.txt")
df2 <- data.table::fread("inst/examples/calsc.info.txt")

df1 %>%
  dplyr::left_join(df2, by = "Position") %>%
  dplyr::filter(Conc >= 4 & Conc <= 4096) %>%
  dplyr::group_by(Gene, Conc) %>%
  dplyr::mutate(mean.cq = mean(Conc)) %>%
  dplyr::ungroup() -> df

CalCurve(cq.table = df1, concen.table = df2,
         lowest.concen = 4, hightest.concen = 4096,
         dilution = 4, by.mean = T) -> p

CalCurve(cq.table = df1, concen.table = df2,
         lowest.concen = 4, hightest.concen = 4096,
         dilu = 4, by = "all") -> p1

p[["figure"]] +
  scale_color_jama()
