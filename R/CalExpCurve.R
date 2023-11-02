#' @name CalExpCurve
#' @author Xiang LI <lixiang117423@@gmail.com>
#'
#' @title Calculate expression using standard curve.
#' @description Calculate expression using standard curve.
#'
#' @param cq.table The data frame of the position and Cq value.
#' @param design.table The data frame of the position and corresponding information.
#' @param correction Correct expression value by reference gene.
#' @param ref.gene The name of reference gene.
#' @param stat.method Statistical method.
#' @param ref.group The name of reference group.
#' @param fig.type Output image type, `box` represents `boxplot`, `bar` represents `barplot`.
#' @param fig.ncol Number of columes of figure.
#'
#' @importFrom dplyr left_join filter group_by mutate ungroup
#' @importFrom stats lm
#' @importFrom broom glance
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth labs theme_bw
#' @importFrom ggpmisc stat_poly_eq
#'
#' @export
#'
#' @return A list contain a table and a figure.
#'
#' @examples
#' df1.path = system.file("examples", "cal.exp.curve.cq.txt", package = "qPCRtools")
#' df2.path = system.file("examples", "cal.expre.curve.sdc.txt", package = "qPCRtools")
#' df3.path = system.file("examples", "cal.exp.curve.design.txt", package = "qPCRtools")
#'
#' cq.table = data.table::fread(df1.path)
#' curve.table = data.table::fread(df2.path)
#' design.table = data.table::fread(df3.path)
#'
#' CalExpCurve(
#'   cq.table,
#'   curve.table,
#'   design.table,
#'   correction = TRUE,
#'   ref.gene = "OsUBQ",
#'   stat.method = "t.test",
#'   ref.group = "CK",
#'   fig.type = "box",
#'   fig.ncol = NULL) -> res
#'
#' res[["table"]]
#' res[["figure"]]
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
  "Cq",
  "max.Cq",
  "min.Cq",
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
  "max.temp"
))
CalExpCurve <- function(cq.table,
                        curve.table,
                        design.table,
                        correction = TRUE,
                        ref.gene = "OsUBQ",
                        stat.method = "t.test",
                        ref.group = "CK",
                        fig.type = "box",
                        fig.ncol = NULL) {
  # merge data
  cq.table %>%
    dplyr::left_join(design.table, by = "Position") %>%
    dplyr::left_join(curve.table, by = "Gene") %>%
    dplyr::mutate(out = dplyr::case_when(
      Cq > max.Cq | Cq < min.Cq ~ "yes",
      TRUE ~ "no"
    )) %>%
    dplyr::mutate(expre = (Cq - Intercept) / Slope) -> df

  # print warning message
  df.out <- df %>%
    dplyr::filter(out == "yes")
  if (dim(df.out)[1] != 0) {
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

  # statistics
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
      dplyr::left_join(df.stat, by = "temp") -> df
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
      dplyr::left_join(df.stat, by = "temp") -> df
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
  }

  # plot
  df %>%
    dplyr::group_by(Gene, Treatment) %>%
    dplyr::mutate(
      mean.expre = mean(expre),
      sd.expre = stats::sd(expre),
      n = dplyr::n(),
      se = sd.expre / sqrt(n)
    ) %>%
    dplyr::ungroup() -> df.plot
  if (fig.type == "box") {
    df.plot %>%
      ggplot2::ggplot(ggplot2::aes(Treatment, expre, fill = Treatment)) +
      ggplot2::geom_boxplot(width = 0.6) +
      ggplot2::facet_wrap(. ~ Gene, scales = "free_y", ncol = fig.ncol) +
      ggplot2::geom_text(ggplot2::aes(Treatment, mean.expre, label = "."),
        check_overlap = TRUE, size = 15, color = "red"
      ) +
      ggplot2::geom_text(ggplot2::aes(Treatment, min(expre), label = signif),
        check_overlap = TRUE, size = 3, color = "black"
      ) +
      ggthemes::theme_pander() +
      ggplot2::labs(y = "Relative expression") +
      ggplot2::theme(
        legend.position = "none",
        strip.text.x = ggplot2::element_text(face = "italic")
      ) -> p
  } else if (fig.type == "bar") {
    df.plot %>%
      dplyr::group_by(Gene) %>%
      dplyr::mutate(max.temp = max(expre)) %>%
      ggplot2::ggplot(ggplot2::aes(Treatment, mean.expre / n, fill = Treatment)) +
      ggplot2::geom_bar(stat = "identity", width = 0.6) +
      ggplot2::geom_errorbar(ggplot2::aes(Treatment,
        ymin = mean.expre - sd.expre,
        ymax = mean.expre + sd.expre
      ),
      width = 0.2
      ) +
      ggplot2::geom_jitter(ggplot2::aes(Treatment, expre), width = 0.1, alpha = 0.4) +
      ggplot2::geom_hline(ggplot2::aes(yintercept = max.temp * 1.1), color = NA) +
      ggplot2::facet_wrap(. ~ Gene, scales = "free_y", ncol = fig.ncol) +
      ggplot2::geom_text(ggplot2::aes(Treatment, max.temp * 1.08, label = signif),
        check_overlap = TRUE, size = 4, color = "black"
      ) +
      ggthemes::theme_pander() +
      ggplot2::labs(y = "Relative expression") +
      ggplot2::theme(
        legend.position = "none",
        strip.text.x = ggplot2::element_text(face = "italic")
      ) -> p
  }

  res <- list(table = df.plot, figure = p)
  return(res)
}
