globalVariables(c(
  "Position", "Cq", "Group", "Gene", "BioRep",
  "cq", "group", "gene", "biorep", "mean.cq",
  "expre", "n", "mean.expre", "sd.expre", "se.expre"
))

#' Calculate Expression Using 2-dCt Method
#'
#' Calculate relative gene expression using the 2-dCt method
#' with a reference gene for normalization.
#'
#' @param cq.table A data frame containing position and Cq values.
#'   Must have columns: Position, Gene, Cq.
#' @param design.table A data frame containing position and group information.
#'   Must have columns: Position, Group, BioRep.
#' @param ref.gene Character. The name of the reference gene (default: "Actin").
#'
#' @return A data frame with expression values, including columns:
#'   position, cq, group, gene, biorep, mean.cq, expre, n, mean.expre,
#'   sd.expre, se.expre.
#'
#' @importFrom magrittr %>%
#' @importFrom stats sd
#'
#' @export
#'
#' @examples
#' \dontrun{
#' df1.path <- system.file("examples", "dct.cq.txt", package = "qPCRtools")
#' df2.path <- system.file("examples", "dct.design.txt", package = "qPCRtools")
#' cq.table <- read.table(df1.path, sep = ",", header = TRUE)
#' design.table <- read.table(df2.path, sep = ",", header = TRUE)
#' res <- CalExp2dCt(cq.table, design.table, ref.gene = "Actin")
#' head(res)
#' }
#'
#' @author Xiang LI <lixiang117423@gmail.com>
CalExp2dCt <- function(cq.table,
                       design.table,
                       ref.gene = "Actin") {
  if (!is.data.frame(cq.table)) {
    stop("'cq.table' must be a data frame")
  }
  if (!is.data.frame(design.table)) {
    stop("'design.table' must be a data frame")
  }

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
      sd.expre = stats::sd(expre),
      se.expre = sd.expre / sqrt(n)
    ) %>%
    dplyr::ungroup()
}
