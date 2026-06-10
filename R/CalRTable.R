globalVariables(c("sample", "concentration", "all", "volume.RNA", "mean"))

#' Calculate RNA Volume for Reverse Transcription
#'
#' The first step of qPCR is usually the preparation of cDNA.
#' This function calculates the volume of RNA needed for reverse
#' transcription based on RNA concentration.
#'
#' @param data A data frame containing sample names and concentration
#'   values (default unit: ng/uL). Must have columns: sample, concentration.
#' @param template A data frame containing reverse transcription information.
#'   Must have a column called `all`.
#' @param rna_weight Numeric. RNA weight required for reverse transcription
#'   in micrograms (default: 1).
#'
#' @return A data frame with calculated RNA and water volumes for each sample.
#'
#' @importFrom magrittr %>%
#' @importFrom dplyr select rename group_by mutate ungroup summarise
#'
#' @export
#'
#' @examples
#' \dontrun{
#' df.1.path <- system.file("examples", "crtv.data.txt", package = "qPCRtools")
#' df.2.path <- system.file("examples", "crtv.template.txt", package = "qPCRtools")
#' df.1 <- read.table(df.1.path, sep = "\t", header = TRUE)
#' df.2 <- read.table(df.2.path, sep = "\t", header = TRUE)
#' result <- CalRTable(data = df.1, template = df.2, rna_weight = 2)
#' head(result)
#' }
#'
#' @author Xiang LI <lixiang117423@gmail.com>
CalRTable <- function(data, template, rna_weight = 1) {
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame")
  }
  if (!is.data.frame(template)) {
    stop("'template' must be a data frame")
  }
  if (!is.numeric(rna_weight) || rna_weight <= 0) {
    stop("'rna_weight' must be a positive number")
  }

  df.1 <- template * rna_weight

  sum.temp <- rowSums(df.1[1, ]) - df.1$all

  df.2 <- data %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise(mean = mean(concentration)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(volume.RNA = rna_weight / mean * 1000) %>%
    cbind(df.1) %>%
    dplyr::mutate(volume.h2o = all - sum.temp - volume.RNA)

  return(df.2)
}
