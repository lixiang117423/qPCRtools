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
#' @param dilution Dilution factor used when the standard curve was fitted; it is
#'   the base of the log transform in `CalCurve()` (default 4). Expression is
#'   reported as concentration, i.e. `dilution^((Cq - Intercept) / Slope)`.
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
#' cq.table = read.table(df1.path, header = TRUE)
#' curve.table = read.table(df2.path, sep = "\t", header = TRUE)
#' design.table = read.table(df3.path, header = TRUE)
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
  "max.temp",
  "RawPValue",
  "PValue",
  "PAdjustMethod",
  "contrast",
  "letters"
))
CalExpCurve <- function(cq.table,
                        curve.table,
                        design.table,
                        correction = TRUE,
                        ref.gene = "OsUBQ",
                        stat.method = "t.test",
                        ref.group = "CK",
                        fig.type = "box",
                        fig.ncol = NULL,
                        dilution = 4) {
  # merge data
  cq.table %>%
    dplyr::left_join(design.table, by = "Position") %>%
    dplyr::left_join(curve.table, by = "Gene") %>%
    dplyr::mutate(out = dplyr::case_when(
      Cq > max.Cq | Cq < min.Cq ~ "yes",
      TRUE ~ "no"
    )) %>%
    # 反对数（2026-10-15 改正）：(Cq - Intercept) / Slope 是**对数浓度**
    # （CalCurve 里 Conc 先取 log(Conc, base = dilution)），必须还原成浓度再用于
    # 内参校正与统计；原实现直接拿对数值相除，两个对数相除没有意义。
    dplyr::mutate(expre = dilution^((Cq - Intercept) / Slope)) -> df

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
      # 保留 BioRep：技术重复层级（见下方）需要它把技术孔归并到生物学重复
      dplyr::select(dplyr::any_of(c("Treatment", "Gene", "BioRep")), expre) -> df
  } else {
    # 同上：保留 BioRep，否则组间检验会退回按孔计算（伪重复）
    df %>%
      dplyr::select(dplyr::any_of(c("Treatment", "Gene", "BioRep")), expre) -> df
  }

  # 技术重复层级（2026-10-15，导师审稿意见第 5 条）：
  #   技术重复只用来估计同一个生物学样本的测量值；组间检验必须以**生物学重复**为单位，
  #   否则同一 BioRep 的多个技术孔会被当成独立样本（伪重复），p 值被算得过小。
  #   设计表含 BioRep 时先在 (Gene, Treatment, BioRep) 内平均技术孔，再以生物学重复为单位
  #   做统计与作图；不含 BioRep 时保持按孔计算并明确告警。
  if ("BioRep" %in% colnames(df)) {
    df.rep <- df %>%
      dplyr::group_by(Gene, Treatment, BioRep) %>%
      dplyr::summarise(expre = mean(expre), .groups = "drop")
  } else {
    warning(paste0(
      "CalExpCurve: the design table has no 'BioRep' column, so technical ",
      "replicates are treated as independent observations ",
      "(pseudoreplication); add BioRep to test at the biological-replicate level."
    ))
    df.rep <- df %>% dplyr::select(Treatment, Gene, expre)
  }

  # statistics（以生物学重复为单位）
  if (stat.method == "t.test") {
    df.rep %>%
      dplyr::group_by(Gene) %>%
      rstatix::t_test(expre ~ Treatment, ref.group = ref.group) %>%
      dplyr::select(Gene, group2, p) %>%
      dplyr::mutate(temp = paste0(Gene, group2)) %>%
      dplyr::select(temp, p) %>%
      dplyr::mutate(signif = dplyr::case_when(
        p < 0.001 ~ "***",
        p < 0.01 ~ "**",
        p < 0.05 ~ "*",
        TRUE ~ "NS"
      ), RawPValue = p, PValue = p, PAdjustMethod = "none") %>%
      dplyr::select(temp, signif, RawPValue, PValue, PAdjustMethod) -> df.stat

    df.rep %>%
      dplyr::mutate(temp = paste0(Gene, Treatment)) %>%
      dplyr::left_join(df.stat, by = "temp") -> df.rep
  } else if (stat.method == "wilcox.test") {
    df.rep %>%
      dplyr::group_by(Gene) %>%
      rstatix::wilcox_test(expre ~ Treatment, ref.group = ref.group) %>%
      dplyr::select(Gene, group2, p) %>%
      dplyr::mutate(temp = paste0(Gene, group2)) %>%
      dplyr::select(temp, p) %>%
      dplyr::mutate(signif = dplyr::case_when(
        p < 0.001 ~ "***",
        p < 0.01 ~ "**",
        p < 0.05 ~ "*",
        TRUE ~ "NS"
      ), RawPValue = p, PValue = p, PAdjustMethod = "none") %>%
      dplyr::select(temp, signif, RawPValue, PValue, PAdjustMethod) -> df.stat

    df.rep %>%
      dplyr::mutate(temp = paste0(Gene, Treatment)) %>%
      dplyr::left_join(df.stat, by = "temp") -> df.rep
  } else if (stat.method == "anova") {
    df.stat <- NULL
    for (i in unique(df.rep$Gene)) {
      df.sub <- df.rep %>%
        dplyr::filter(Gene == i) %>%
        dplyr::mutate(Treatment = factor(Treatment))
      fit <- stats::aov(expre ~ Treatment, data = df.sub)
      tuk <- multcomp::glht(fit, linfct = multcomp::mcp(Treatment = "Tukey"))
      # cld.glht 的 level 是显著性水平（默认 0.05）；原写成 0.95 = α=0.95，几乎全标成不同字母
      # （2026-10-15 改正；decreasing 原拼写 ddecreasing 被忽略）
      letters <- multcomp::cld(
        tuk, level = 0.05, decreasing = TRUE
      )[["mcletters"]][["Letters"]] %>%
        as.data.frame() %>%
        dplyr::mutate(Gene = i) %>%
        tibble::rownames_to_column(var = "Treatment") %>%
        magrittr::set_colnames(c("Treatment", "signif", "Gene"))

      raw.test <- summary(
        tuk, test = multcomp::adjusted("none")
      )[["test"]]
      raw.p <- if (is.list(raw.test) && !is.null(raw.test[["pvalues"]])) {
        as.numeric(raw.test[["pvalues"]])
      } else if (is.null(dim(raw.test))) {
        as.numeric(raw.test)
      } else {
        raw.test[, "pvalue"]
      }
      adjusted.test <- summary(tuk)[["test"]]
      adjusted.p <- if (
        is.list(adjusted.test) && !is.null(adjusted.test[["pvalues"]])
      ) {
        as.numeric(adjusted.test[["pvalues"]])
      } else if (is.null(dim(adjusted.test))) {
        as.numeric(adjusted.test)
      } else {
        adjusted.test[, "pvalue"]
      }
      contrasts <- if (
        is.list(adjusted.test) && !is.null(adjusted.test[["coefficients"]])
      ) {
        names(adjusted.test[["coefficients"]])
      } else {
        rownames(adjusted.test)
      }
      p.values <- data.frame(
        contrast = contrasts,
        RawPValue = raw.p,
        PValue = adjusted.p,
        stringsAsFactors = FALSE
      ) %>%
        dplyr::mutate(Treatment = dplyr::case_when(
          endsWith(contrast, paste0(" - ", ref.group)) ~
            sub(paste0(" - ", ref.group, "$"), "", contrast),
          startsWith(contrast, paste0(ref.group, " - ")) ~
            sub(paste0("^", ref.group, " - "), "", contrast),
          TRUE ~ NA_character_
        )) %>%
        dplyr::filter(!is.na(Treatment)) %>%
        dplyr::mutate(PAdjustMethod = "Tukey") %>%
        dplyr::select(Treatment, RawPValue, PValue, PAdjustMethod)

      letters %>%
        dplyr::left_join(p.values, by = "Treatment") %>%
        dplyr::mutate(temp = paste0(Gene, Treatment)) %>%
        dplyr::select(temp, signif, RawPValue, PValue, PAdjustMethod) %>%
        rbind(df.stat) -> df.stat
    }
    df.rep %>%
      dplyr::mutate(temp = paste0(Gene, Treatment)) %>%
      dplyr::left_join(df.stat, by = "temp") %>%
      dplyr::select(-temp) -> df.rep
  }

  # plot（同样以生物学重复为单位）
  df.rep %>%
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
      dplyr::distinct(Treatment, Gene, .keep_all = TRUE) %>%
      dplyr::group_by(Gene) %>%
      dplyr::mutate(max.temp = max(expre)) %>%
      # 柱高应为组均值；原实现写成 mean.expre / n（除以重复数没有意义），2026-10-15 改正
      ggplot2::ggplot(ggplot2::aes(Treatment, mean.expre, fill = Treatment)) +
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
