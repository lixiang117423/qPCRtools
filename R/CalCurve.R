globalVariables(c(
  "Conc", "Gene", "Cq", "mean.cq", "max.Cq", "min.Cq",
  "Position", "rr.label", "p.value.label"
))

#' Standard Curve Calculation
#'
#' Calculate the standard curve and obtain the amplification efficiency
#' of primer(s). Based on the amplification efficiency, we can determine
#' which method to use for expression level calculation.
#'
#' @param cq.table A data frame containing position and Cq values.
#'   Must have columns: Position, Gene, Cq.
#' @param concen.table A data frame containing position and concentration.
#'   Must have columns: Position, Conc.
#' @param highest.concen Numeric. The highest concentration for calculation.
#' @param lowest.concen Numeric. The lowest concentration for calculation.
#' @param dilution Numeric. Dilution factor of cDNA template (default: 4).
#' @param by.mean Logical. Calculate by mean Cq value or not (default: TRUE).
#'
#' @return A list containing:
#'   \item{table}{Data frame with standard curve parameters per gene
#'     (Formula, Slope, Intercept, R2, P.value, max.Cq, min.Cq, E, Date)}
#'   \item{figure}{ggplot object of the standard curve}
#'
#' @importFrom magrittr %>%
#' @importFrom stats lm sd
#' @importFrom broom glance
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth labs theme_bw
#' @importFrom ggpmisc stat_poly_eq
#' @importFrom kableExtra kable_styling kable
#'
#' @export
#'
#' @examples
#' \dontrun{
#' df.1.path <- system.file("examples", "calsc.cq.txt", package = "qPCRtools")
#' df.2.path <- system.file("examples", "calsc.info.txt", package = "qPCRtools")
#' df.1 <- read.table(df.1.path, header = TRUE)
#' df.2 <- read.table(df.2.path, header = TRUE)
#' res <- CalCurve(
#'   cq.table = df.1,
#'   concen.table = df.2,
#'   lowest.concen = 4,
#'   highest.concen = 4096,
#'   dilution = 4,
#'   by.mean = TRUE
#' )
#' res[["table"]]
#' res[["figure"]]
#' }
#'
#' @author Xiang LI <lixiang117423@gmail.com>
CalCurve <- function(cq.table,
                     concen.table,
                     highest.concen,
                     lowest.concen,
                     dilution = 4,
                     by.mean = TRUE) {
  # input validation
  if (!is.data.frame(cq.table)) {
    stop("'cq.table' must be a data frame")
  }
  if (!is.data.frame(concen.table)) {
    stop("'concen.table' must be a data frame")
  }
  if (!is.numeric(highest.concen) || !is.numeric(lowest.concen)) {
    stop("'highest.concen' and 'lowest.concen' must be numeric")
  }

  cq.table %>%
    dplyr::left_join(concen.table, by = "Position") %>%
    dplyr::filter(Conc >= lowest.concen & Conc <= highest.concen) %>%
    dplyr::group_by(Gene, Conc) %>%
    dplyr::mutate(
      mean.cq = mean(Cq),
      Conc = log(Conc, base = dilution)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(Gene) %>%
    dplyr::mutate(
      max.Cq = max(Cq),
      min.Cq = min(Cq)
    ) %>%
    dplyr::ungroup() -> df

  if (isTRUE(by.mean)) {
    res <- cal_curve_by_mean(df, dilution)
  } else {
    res <- cal_curve_by_raw(df, dilution)
  }

  return(res)
}

#' Build standard curve model using mean Cq values
#' @keywords internal
cal_curve_by_mean <- function(df, dilution) {
  fit.res <- NULL

  for (i in unique(df$Gene)) {
    df.sub <- df %>%
      dplyr::filter(Gene == i)

    fit <- stats::lm(mean.cq ~ Conc, data = df.sub)
    intercept <- fit[["coefficients"]][["(Intercept)"]] %>%
      round(2)
    slope <- fit[["coefficients"]][["Conc"]] %>%
      round(2)

    formula <- paste0("y = ", slope, "*Conc", " + ", intercept)

    r.2 <- broom::glance(fit)[1, 1] %>%
      round(4) %>%
      as.numeric()

    p.value <- broom::glance(fit)[1, 5] %>%
      round(5) %>%
      as.numeric()

    df.temp <- data.frame(
      Gene = i,
      Formula = formula,
      Slope = slope,
      Intercept = intercept,
      R2 = r.2,
      P.value = p.value,
      max.Cq = unique(df.sub$max.Cq),
      min.Cq = unique(df.sub$min.Cq),
      E = round(dilution^(-1 / slope) - 1, 3),
      Date = as.character(Sys.Date())
    )
    fit.res <- rbind(fit.res, df.temp) -> res.table
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(Conc, mean.cq, color = Gene)) +
    ggplot2::geom_point() +
    ggplot2::geom_smooth(method = "lm", show.legend = FALSE) +
    ggpmisc::stat_poly_eq(
      ggplot2::aes(
        label = paste(
          ggplot2::after_stat(rr.label),
          ggplot2::after_stat(p.value.label),
          sep = "~~~~"
        )
      ),
      show.legend = FALSE,
      formula = y ~ x,
      parse = TRUE,
      rr.digits = 5,
      coef.digits = 3,
      label.x = c(0.05),
      label.y = seq(
        0.05,
        0.06 * (length(unique(df$Gene)) + 1),
        0.06
      )
    ) +
    ggplot2::labs(
      x = "log Relative Concentration",
      y = "Mean Cq Value"
    ) +
    ggplot2::theme_bw()

  list(table = res.table, figure = p)
}

#' Build standard curve model using raw Cq values
#' @keywords internal
cal_curve_by_raw <- function(df, dilution) {
  fit.res <- NULL

  for (i in unique(df$Gene)) {
    df.sub <- df %>%
      dplyr::filter(Gene == i)

    fit <- stats::lm(Cq ~ Conc, data = df.sub)
    intercept <- fit[["coefficients"]][["(Intercept)"]] %>%
      round(2)
    slope <- fit[["coefficients"]][["Conc"]] %>%
      round(2)

    formula <- paste0("y = ", slope, "*Conc", " + ", intercept)

    r.2 <- broom::glance(fit)[1, 1] %>%
      round(4) %>%
      as.numeric()

    p.value <- broom::glance(fit)[1, 5] %>%
      round(5) %>%
      as.numeric()

    df.temp <- data.frame(
      Gene = i,
      Formula = formula,
      Slope = slope,
      Intercept = intercept,
      R2 = r.2,
      P.value = p.value,
      max.Cq = unique(df.sub$max.Cq),
      min.Cq = unique(df.sub$min.Cq),
      E = round(dilution^(-1 / slope) - 1, 3),
      Date = as.character(Sys.Date())
    )
    fit.res <- rbind(fit.res, df.temp) -> res.table
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(Conc, Cq, color = Gene)) +
    ggplot2::geom_point() +
    ggplot2::geom_smooth(method = "lm", show.legend = FALSE) +
    ggpmisc::stat_poly_eq(
      ggplot2::aes(
        label = paste(
          ggplot2::after_stat(rr.label),
          ggplot2::after_stat(p.value.label),
          sep = "~~~~"
        )
      ),
      show.legend = FALSE,
      formula = y ~ x,
      parse = TRUE,
      rr.digits = 5,
      coef.digits = 3,
      label.x = c(0.05),
      label.y = seq(
        0.05,
        0.06 * (length(unique(df$Gene)) + 1),
        0.06
      )
    ) +
    ggplot2::labs(
      x = "log Relative Concentration",
      y = "Cq Value"
    ) +
    ggplot2::theme_bw()

  list(table = res.table, figure = p)
}
