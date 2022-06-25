#' Standard Curve Calculation
#' @description A shiny Module.
#'
#' @param cq.table The data frame of the position and Cq value.
#' @param concen.table The data frame of the position and concentration.
#' @param highest.concen The highest concentration.
#' @param lowest.concen The lowest concentration.
#' @param dilution Dilution factor of cDNA template.
#' @param by.mean Calculation by mean Cq value or not.
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

CalCurve <- function(cq.table,
                     concen.table,
                     highest.concen,
                     lowest.concen,
                     dilution = 4,
                     by.mean = TRUE) {
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
    # build model
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

    ggplot2::ggplot(df, ggplot2::aes(Conc, mean.cq, color = Gene)) +
      ggplot2::geom_point() +
      ggplot2::geom_smooth(method = "lm", show.legend = FALSE) +
      ggpmisc::stat_poly_eq(
        ggplot2::aes(
          label = paste( # ..eq.label..,
            ..rr.label..,
            ..p.value.label..,
            sep = "~~~~"
          )
        ),
        show.legend = FALSE,
        formula = y ~ x,
        parse = T,
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
      ggplot2::theme_bw() -> p

    res <- list(table = res.table, figure = p)
  } else {
    # build model
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

    ggplot2::ggplot(df, ggplot2::aes(Conc, Cq, color = Gene)) +
      ggplot2::geom_point() +
      ggplot2::geom_smooth(method = "lm", show.legend = FALSE) +
      ggpmisc::stat_poly_eq(
        ggplot2::aes(
          label = paste( # ..eq.label..,
            ..rr.label..,
            ..p.value.label..,
            sep = "~~~~"
          )
        ),
        show.legend = FALSE,
        formula = y ~ x,
        parse = T,
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
      ggplot2::theme_bw() -> p

    res <- list(table = res.table, figure = p)
  }
  return(res)
}
