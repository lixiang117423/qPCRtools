globalVariables(c(
  "Position", "Cq", "Gene",
  "max.Cq", "min.Cq", "Intercept", "Slope",
  "Treatment", "expre", "mean.ref",
  "mean.expre", "sd.expre", "se", "n",
  "signif", "temp", "group2", "out", "max.temp", "p"
))

#' Calculate Expression Using Standard Curve
#'
#' Calculate relative gene expression using a standard curve method
#' with optional reference gene correction and statistical testing.
#'
#' @param cq.table A data frame containing position and Cq values.
#'   Must have columns: Position, Gene, Cq.
#' @param curve.table A data frame with standard curve parameters per gene.
#'   Must have columns: Gene, Slope, Intercept, max.Cq, min.Cq.
#' @param design.table A data frame containing position and group information.
#'   Must have columns: Position, Treatment, Gene.
#' @param correction Logical. Correct expression by reference gene
#'   (default: TRUE).
#' @param ref.gene Character. The name of the reference gene (default: "OsUBQ").
#' @param stat.method Character. Statistical method for group comparison.
#'   One of "t.test", "wilcox.test", or "anova" (default: "t.test").
#' @param ref.group Character. The name of the reference/control group
#'   (default: "CK").
#' @param fig.type Character. Plot type: "box" for boxplot, "bar" for barplot
#'   (default: "box").
#' @param fig.ncol Integer. Number of columns in facet plot (default: NULL).
#'
#' @return A list containing:
#'   \item{table}{Data frame with expression values and statistics}
#'   \item{figure}{ggplot object}
#'
#' @importFrom magrittr %>%
#' @importFrom stats sd aov
#' @importFrom broom glance
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth labs theme_bw
#' @importFrom ggpmisc stat_poly_eq
#'
#' @export
#'
#' @examples
#' \dontrun{
#' df1.path <- system.file("examples", "cal.exp.curve.cq.txt", package = "qPCRtools")
#' df2.path <- system.file("examples", "cal.expre.curve.sdc.txt", package = "qPCRtools")
#' df3.path <- system.file("examples", "cal.exp.curve.design.txt", package = "qPCRtools")
#' cq.table <- read.table(df1.path, header = TRUE)
#' curve.table <- read.table(df2.path, sep = "\t", header = TRUE)
#' design.table <- read.table(df3.path, header = TRUE)
#' res <- CalExpCurve(
#'   cq.table, curve.table, design.table,
#'   correction = TRUE,
#'   ref.gene = "OsUBQ",
#'   stat.method = "t.test",
#'   ref.group = "CK",
#'   fig.type = "box",
#'   fig.ncol = NULL
#' )
#' res[["table"]]
#' res[["figure"]]
#' }
#'
#' @author Xiang LI <lixiang117423@gmail.com>
CalExpCurve <- function(cq.table,
                        curve.table,
                        design.table,
                        correction = TRUE,
                        ref.gene = "OsUBQ",
                        stat.method = "t.test",
                        ref.group = "CK",
                        fig.type = "box",
                        fig.ncol = NULL) {
  if (!is.data.frame(cq.table)) {
    stop("'cq.table' must be a data frame")
  }
  if (!is.data.frame(curve.table)) {
    stop("'curve.table' must be a data frame")
  }
  if (!is.data.frame(design.table)) {
    stop("'design.table' must be a data frame")
  }
  if (!stat.method %in% c("t.test", "wilcox.test", "anova")) {
    stop("'stat.method' must be one of 't.test', 'wilcox.test', or 'anova'")
  }
  if (!fig.type %in% c("box", "bar")) {
    stop("'fig.type' must be 'box' or 'bar'")
  }

  # merge data
  cq.table %>%
    dplyr::left_join(design.table, by = "Position") %>%
    dplyr::left_join(curve.table, by = "Gene") %>%
    dplyr::mutate(out = dplyr::case_when(
      Cq > max.Cq | Cq < min.Cq ~ "yes",
      TRUE ~ "no"
    )) %>%
    dplyr::mutate(expre = (Cq - Intercept) / Slope) -> df

  # warn about out-of-range Cq values
  df.out <- df %>%
    dplyr::filter(out == "yes")
  if (nrow(df.out) != 0) {
    warning(paste0("Cq of ", as.character(df.out$Position), " out of curve range!"))
  }

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
  } else {
    df %>%
      dplyr::select(Treatment, Gene, expre) -> df
  }

  # statistical tests
  df <- cal_curve_stat_test(df, stat.method, ref.group)

  # summary for plot
  df %>%
    dplyr::group_by(Gene, Treatment) %>%
    dplyr::mutate(
      mean.expre = mean(expre),
      sd.expre = stats::sd(expre),
      n = dplyr::n(),
      se = sd.expre / sqrt(n)
    ) %>%
    dplyr::ungroup() -> df.plot

  p <- build_curve_plot(df.plot, fig.type, fig.ncol)

  res <- list(table = df.plot, figure = p)
  return(res)
}

#' Run statistical tests for curve-based expression
#' @keywords internal
cal_curve_stat_test <- function(df, stat.method, ref.group) {
  if (stat.method == "t.test") {
    df %>%
      dplyr::group_by(Gene) %>%
      rstatix::t_test(expre ~ Treatment, ref.group = ref.group) %>%
      dplyr::select(Gene, group2, p) %>%
      dplyr::mutate(temp = paste0(Gene, group2)) %>%
      dplyr::select(temp, p) %>%
      dplyr::mutate(signif = dplyr::case_when(
        p < 0.001 ~ "***",
        p > 0.001 & p < 0.01 ~ "**",
        p > 0.01 & p < 0.05 ~ "*",
        TRUE ~ "NS"
      )) %>%
      dplyr::select(temp, signif) -> df.stat

    df %>%
      dplyr::mutate(temp = paste0(Gene, Treatment)) %>%
      dplyr::left_join(df.stat, by = "temp")
  } else if (stat.method == "wilcox.test") {
    df %>%
      dplyr::group_by(Gene) %>%
      rstatix::wilcox_test(expre ~ Treatment, ref.group = ref.group) %>%
      dplyr::select(Gene, group2, p) %>%
      dplyr::mutate(temp = paste0(Gene, group2)) %>%
      dplyr::select(temp, p) %>%
      dplyr::mutate(signif = dplyr::case_when(
        p < 0.001 ~ "***",
        p > 0.001 & p < 0.01 ~ "**",
        p > 0.01 & p < 0.05 ~ "*",
        TRUE ~ "NS"
      )) %>%
      dplyr::select(temp, signif) -> df.stat

    df %>%
      dplyr::mutate(temp = paste0(Gene, Treatment)) %>%
      dplyr::left_join(df.stat, by = "temp")
  } else if (stat.method == "anova") {
    df.stat <- NULL
    for (i in unique(df$Gene)) {
      df.sub <- df %>%
        dplyr::filter(Gene == i) %>%
        dplyr::mutate(Treatment = factor(Treatment))
      fit <- stats::aov(expre ~ Treatment, data = df.sub)
      tuk <- multcomp::glht(fit, linfct = multcomp::mcp(Treatment = "Tukey"))
      multcomp::cld(tuk, level = 0.95, ddecreasing = TRUE)[["mcletters"]][["Letters"]] %>%
        as.data.frame() %>%
        dplyr::mutate(Gene = i) %>%
        tibble::rownames_to_column(var = "Treatment") %>%
        magrittr::set_colnames(c("Treatment", "signif", "Gene")) %>%
        dplyr::select(Gene, Treatment, signif) %>%
        dplyr::mutate(temp = paste0(Gene, Treatment)) %>%
        dplyr::select(temp, signif) %>%
        rbind(df.stat) -> df.stat

      df %>%
        dplyr::mutate(temp = paste0(Gene, Treatment)) %>%
        dplyr::left_join(df.stat, by = "temp") -> df
    }
    df
  }
}

#' Build curve expression plot
#' @keywords internal
build_curve_plot <- function(df.plot, fig.type, fig.ncol) {
  if (fig.type == "box") {
    df.plot %>%
      ggplot2::ggplot(ggplot2::aes(Treatment, expre, fill = Treatment)) +
      ggplot2::geom_boxplot(width = 0.6) +
      ggplot2::facet_wrap(. ~ Gene, scales = "free_y", ncol = fig.ncol) +
      ggplot2::geom_text(
        ggplot2::aes(Treatment, mean.expre, label = "."),
        check_overlap = TRUE, size = 15, color = "red"
      ) +
      ggplot2::geom_text(
        ggplot2::aes(Treatment, min(expre), label = signif),
        check_overlap = TRUE, size = 3, color = "black"
      ) +
      ggthemes::theme_pander() +
      ggplot2::labs(y = "Relative expression") +
      ggplot2::theme(
        legend.position = "none",
        strip.text.x = ggplot2::element_text(face = "italic")
      )
  } else if (fig.type == "bar") {
    df.plot %>%
      dplyr::group_by(Gene) %>%
      dplyr::mutate(max.temp = max(expre)) %>%
      ggplot2::ggplot(ggplot2::aes(Treatment, mean.expre / n, fill = Treatment)) +
      ggplot2::geom_bar(stat = "identity", width = 0.6) +
      ggplot2::geom_errorbar(ggplot2::aes(
        Treatment,
        ymin = mean.expre - sd.expre,
        ymax = mean.expre + sd.expre
      ), width = 0.2) +
      ggplot2::geom_jitter(ggplot2::aes(Treatment, expre), width = 0.1, alpha = 0.4) +
      ggplot2::geom_hline(ggplot2::aes(yintercept = max.temp * 1.1), color = NA) +
      ggplot2::facet_wrap(. ~ Gene, scales = "free_y", ncol = fig.ncol) +
      ggplot2::geom_text(
        ggplot2::aes(Treatment, max.temp * 1.08, label = signif),
        check_overlap = TRUE, size = 4, color = "black"
      ) +
      ggthemes::theme_pander() +
      ggplot2::labs(y = "Relative expression") +
      ggplot2::theme(
        legend.position = "none",
        strip.text.x = ggplot2::element_text(face = "italic")
      )
  }
}
