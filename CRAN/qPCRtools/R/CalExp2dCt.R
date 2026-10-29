#' @name CalExp2dCt
#' @author Xiang LI <lixiang117423@@gmail.com>
#'
#' @title Calculate expression using standard curve.
#' @description Calculate expression using standard curve.
#'
#' @param cq.table The data frame of the position and cq value.
#' @param design.table The data frame of the position and corresponding information.
#' @param ref.gene The name of reference gene.
#'
#' @export
#' @return A list contain a table and a figure.
#' @examples
#' df1.path <- system.file("examples", "dct.cq.txt", package = "qPCRtools")
#' df2.path <- system.file("examples", "dct.design.txt", package = "qPCRtools")
#' cq.table <- read.table(df1.path, sep = ",", header = TRUE)
#' design.table <- read.table(df2.path, sep = ",", header = TRUE)
#' CalExp2dCt(cq.table,
#'            design.table,
#'            ref.gene = "Actin"
#' ) -> res
#'
globalVariables(c(
  "cq.table",
  "curve.table",
  "design.table",
  "correction",
  "ref.gene",
  "stat.method",
  "ref.group",
  "fig.type",
  "fig.ncol",
  "out",
  "cq",
  "max.cq",
  "min.cq",
  "expre",
  "Intercept",
  "Slope",
  "Treatment",
  "element_text",
  "group2",
  "max.temp",
  "mean.expre",
  "mean.ref",
  "n",
  "sd",
  "sd.expre",
  "Treatment",
  "mean.ref",
  "group2",
  "temp",
  "n",
  "sd.expre",
  "n",
  "mean.expre",
  "element_text",
  "max.temp",
  "gene",
  "group",
  "biorep",
  "Target",
  "Reference",
  "ddct1",
  "mean.expression",
  "n.biorep",
  "sd.expression",
  "se.expression",
  "Reference",
  "Target",
  "biorep",
  "ddct1",
  "gene",
  "group",
  "mean.expression",
  "n.biorep",
  "sd.expression",
  "se.expression",
  "BioRep",
  "Eff",
  "Group",
  "TechRep"
))
CalExp2dCt <- function(cq.table,
                       design.table,
                       ref.gene = "Actin") {
  # merge data
  cq.table %>%
    dplyr::left_join(design.table, by = "Position") %>%
    dplyr::rename(
      position = Position,
      cq = Cq,
      group = Group,
      gene = Gene,
      biorep = BioRep
    ) -> df

  # reference gene
  df %>%
    dplyr::filter(gene == ref.gene) %>%
    dplyr::group_by(group, biorep) %>%
    dplyr::mutate(
      mean.cq = mean(cq),
      temp = paste0(group, biorep)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(temp, mean.cq) %>%
    dplyr::distinct_all() -> df.ref

  # target gene
  df %>%
    dplyr::filter(gene != ref.gene) %>%
    dplyr::mutate(temp = paste0(group, biorep)) %>%
    dplyr::left_join(df.ref, by = "temp") %>%
    dplyr::select(-temp) %>%
    dplyr::mutate(expre = 2^(mean.cq - cq)) %>%
    dplyr::group_by(group, gene) %>%
    dplyr::mutate(
      n = dplyr::n(),
      mean.expre = mean(expre),
      sd.expre = sd(expre),
      se.expre = sd.expre / sqrt(n)
    ) %>%
    dplyr::ungroup() -> res

  return(res)
}
