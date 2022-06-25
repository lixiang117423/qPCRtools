#' Calculate expression using standard curve.
#' @author Xiang LI <lixiang117423@@gmail.com>
#' @description A shiny Module.
#'
#' @param cq.table The data frame of the position and Cq value.
#' @param design.table The data frame of the position and corresponding information.
#' @param correction Correct expression value by reference gene.
#' @param ref.gene The name of reference gene.
#' @param stat.method Statistical method.
#' @param ref.group The name of reference group.
#' @param fig.type Calculation by mean Cq value or not.
#'
#' @importFrom dplyr left_join filter group_by mutate ungroup
#' @importFrom stats lm
#' @importFrom broom glance
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth labs theme_bw
#' @importFrom ggpmisc stat_poly_eq
#'
#' @export
#'
#' @examples
#' df.1.path <- system.file("examples", "calsc.cq.txt", package = "qPCRtools")
#' df.2.path <- system.file("examples", "calsc.info.txt", package = "qPCRtools")
#' df.1 <- data.table::fread(df.1.path)
#' df.2 <- data.table::fread(df.2.path)
#' CalCurve(
#'   cq.table = df.1,
#'   concen.table = df.2,
#'   lowest.concen = 4,
#'   highest.concen = 4096,
#'   dilu = 4,
#'   by = "mean"
#' ) -> p
#'
#' p[["table"]]
#' p[["figure"]]
#'
globalVariables(c(
  "cq.table",
  "concen.table",
  "highest.concen",
  "lowest.concen",
  "dilution",
  "by.mean",
  "Conc",
  "Gene",
  "Cq",
  "Formula",
  "Slope",
  "Intercept",
  "R2",
  "P.value",
  "max.Cq",
  "min.Cq",
  "E",
  "Date",
  "..rr.label..",
  "..p.value.label..",
  "mean.cq"
))

CalExpCurve = function(cq.table,
                       design.table,
                       correction,
                       ref.gene,
                       stat.method,
                       ref.group,
                       fig.type){
  # merge data

}

