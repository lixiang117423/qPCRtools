#' @name CalRTable
#' @author Xiang LI <lixiang117423@@gmail.com>
#' @title Calculate RNA volume for reverse transcription.
#' @description The first step of qPCR is usually the preparation of cDNA.
#' We need to calculate the column of RNA for reverse transcription to cDNA.
#' So, if we have the concentration of RNA, we can use the function `CalRTable` to do that.
#'
#' @param data A data.frame contained the sample names and the concentration value. The default unit of concentration is ng/uL.
#' @param template A data.frame contained the information of reverse transcription. In this data.frame there must be a column called `all`.
#' @param RNA.weight RNA weight required for reverse transcription. Default is 1 ug.
#'
#' @importFrom dplyr select rename group_by mutate ungroup
#' @importFrom magrittr %>%
#'
#' @export
#'
#' @examples
#' df.1.path <- system.file("examples", "crtv.data.txt", package = "qPCRtools")
#' df.2.path <- system.file("examples", "crtv.template.txt", package = "qPCRtools")
#' df.1 <- read.table(df.1.path, sep = "\t", header = TRUE)
#' df.2 <- read.table(df.2.path, sep = "\t", header = TRUE)
#' result <- CalRTable(data = df.1, template = df.2, RNA.weight = 2)
#' head(result)
#' @return A list contain a table and a figure.
globalVariables(c("data", "template", "RNA.weight",
                  "df.1", "sum.temp", "sample",
                  "concentration", "volume.RNA",
                  "mean","volume.h2o","all",
                  "sum.temp","volume.RNA",
                  "df.2"))
CalRTable <- function(data, template, RNA.weight = 1) {
  df.1 <- template * RNA.weight

  sum.temp <- rowSums(df.1[1, ]) - df.1$all

  df.2 <- data %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise(mean = mean(concentration)) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(volume.RNA = RNA.weight / mean * 1000) %>%
    cbind(df.1) %>%
    dplyr::mutate(volume.h2o = all - sum.temp - volume.RNA)
  return(df.2)
}


