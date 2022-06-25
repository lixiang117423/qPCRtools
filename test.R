rm(list = ls())

cq.table = data.table::fread("CRAN/qPCRtools/inst/examples/cal.exp.curve.cq.txt", header = TRUE)
curve.table = data.table::fread("CRAN/qPCRtools/inst/examples/cal.expre.curve.sdc.txt", header = TRUE)
design.table = data.table::fread("CRAN/qPCRtools/inst/examples/cal.exp.curve.design.txt", header = TRUE)


cq.table %>%
  dplyr::left_join(design.table, by = "Position") %>%
  dplyr::left_join(curve.table, by = "Gene") %>%
  dplyr::mutate(out = dplyr::case_when(
    Cq > max.Cq | Cq < min.Cq ~ "yes",
    TRUE ~ "no"
  )) %>%
  dplyr::mutate(expre = (Cq - Intercept) / Slope) -> df

if (isTRUE(correction)) {
  df %>%
    dplyr::filter(Gene == ref.gene) %>%
    dplyr::group_by(Treatment) %>%
    dplyr::summarise(mean.ref = mean(expre)) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(df, by = "Treatment") %>%
    dplyr::filter(Gene != ref.gene) %>%
    dplyr::mutate(expre = expre / mean.ref) %>%
    dplyr::select(Treatment, Gene, expre) -> df
}else{
  df %>%
    dplyr::select(Treatment, Gene, expre) -> df
}

