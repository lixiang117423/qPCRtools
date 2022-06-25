#' @name CalRTable
#' @author Xiang LI <lixiang117423@@gmail.com>
#' @title Calculate volume.
#' @description Calculate RNA and other reagent volume required for reverse transcription.
#'
#' @param data A data.frame contained the sample names and the concentration value. The default unit of concentration is ng/uL.
#' @param template A data.frame contained the information of reverse transcription.
#' @param RNA.weight RNA weight required for reverse transcription. Default is 1 ug.
#'
#' @importFrom dplyr select rename group_by mutate ungroup
#' @importFrom magrittr %>%
#'
#' @export
#'
#' @examples
#' df.1.path <- system.file("examples", "crtv.data.txt", package = "qPCRtools")
#' df.2.path <- system.file("examples", "template", package = "qPCRtools")
#' df.1 <- data.table::fread(df.1.path)
#' df.2 <- data.table::fread(df.2.path)
#' result <- CalRTable(data = df.1, template = df.2, RNA.weight = 2)
#' head(result)
#' @return A data frame.
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


