globalVariables(c(
  "Position", "Cq", "Group", "Gene", "BioRep",
  "cq", "group", "gene", "biorep",
  "expression", "mean.expression", "sd.expression", "se.expression",
  "n.biorep", "group2", "signif", "temp",
  "Treatment", "expre", "mean.expre", "sd.expre", "se.expre", "n",
  "max.temp", "is.out", "Target", "Reference", "ddct1", "p"
))

#' Calculate Expression Using 2-ddCt Method
#'
#' Calculate relative gene expression using the 2-ddCt method
#' with a reference gene and reference group for normalization.
#' Supports statistical testing and outlier removal.
#'
#' @param cq_table A data frame containing position and Cq values.
#'   Must have columns: Position, Gene, Cq.
#' @param design_table A data frame containing position and group information.
#'   Must have columns: Position, Group, BioRep.
#' @param ref_gene Character. The name of the reference gene (default: "OsUBQ").
#' @param ref_group Character. The name of the reference/control group
#'   (default: "CK").
#' @param stat_method Character. Statistical method for group comparison.
#'   One of "t.test", "wilcox.test", or "anova" (default: "t.test").
#' @param remove_outliers Logical. Remove outliers using IQR method
#'   (default: TRUE).
#' @param fig_type Character. Plot type: "box" for boxplot, "bar" for barplot
#'   (default: "box").
#' @param fig_ncol Integer. Number of columns in facet plot (default: NULL).
#'
#' @return A list containing:
#'   \item{table}{Data frame with expression values and statistics}
#'   \item{figure}{ggplot object}
#'
#' @importFrom magrittr %>%
#' @importFrom stats sd
#'
#' @export
#'
#' @examples
#' \dontrun{
#' df1.path <- system.file("examples", "ddct.cq.txt", package = "qPCRtools")
#' df2.path <- system.file("examples", "ddct.design.txt", package = "qPCRtools")
#' cq_table <- read.table(df1.path, header = TRUE)
#' design_table <- read.table(df2.path, header = TRUE)
#' res <- CalExp2ddCt(
#'   cq_table,
#'   design_table,
#'   ref_gene = "OsUBQ",
#'   ref_group = "CK",
#'   stat_method = "t.test",
#'   remove_outliers = TRUE,
#'   fig_type = "box",
#'   fig_ncol = NULL
#' )
#' res[["table"]]
#' res[["figure"]]
#' }
#'
#' @author Xiang LI <lixiang117423@gmail.com>
CalExp2ddCt <- function(cq_table,
                        design_table,
                        ref_gene = "OsUBQ",
                        ref_group = "CK",
                        stat_method = "t.test",
                        remove_outliers = TRUE,
                        fig_type = "box",
                        fig_ncol = NULL) {
  if (!is.data.frame(cq_table)) {
    stop("'cq_table' must be a data frame")
  }
  if (!is.data.frame(design_table)) {
    stop("'design_table' must be a data frame")
  }
  if (!stat_method %in% c("t.test", "wilcox.test", "anova")) {
    stop("'stat_method' must be one of 't.test', 'wilcox.test', or 'anova'")
  }
  if (!fig_type %in% c("box", "bar")) {
    stop("'fig_type' must be 'box' or 'bar'")
  }

  res.all <- NULL

  # merge data
  cq_table %>%
    dplyr::left_join(design_table, by = "Position") %>%
    dplyr::rename(
      position = Position,
      cq = Cq,
      group = Group,
      gene = Gene,
      biorep = BioRep
    ) -> df

  # for each target gene
  target.genes <- setdiff(unique(df$gene), ref_gene)

  for (genes in target.genes) {
    df.sub <- df %>%
      dplyr::filter(gene %in% c(genes, ref_gene))

    # reference gene in CK
    df.sub.ck.ref_gene <- df.sub %>%
      dplyr::filter(gene == ref_gene)
    mean.ck.ref_gene <- mean(df.sub.ck.ref_gene$cq)

    df.sub.ck.target.gene <- df.sub %>%
      dplyr::filter(gene != ref_gene)
    mean.ck.target.gene <- mean(df.sub.ck.target.gene$cq)

    dct1 <- mean.ck.target.gene - mean.ck.ref_gene

    # for each treatment
    for (groups in unique(df.sub$group)) {
      df.sub.group <- df.sub %>%
        dplyr::filter(group == groups) %>%
        dplyr::select(biorep, gene, cq) %>%
        dplyr::group_by(gene, biorep) %>%
        dplyr::mutate(cq = mean(cq)) %>%
        dplyr::ungroup() %>%
        dplyr::distinct_all() %>%
        tidyr::pivot_wider(id_cols = "biorep", names_from = "gene", values_from = "cq") %>%
        dplyr::mutate(dct1 = dct1)

      if (ncol(df.sub.group) == 3 & !genes %in% colnames(df.sub.group)) {
        stop(paste0("Data of target gene ", genes, " has some problem, please check it and try again!"))
      }

      # handle column order depending on ref gene position
      if (colnames(df.sub.group)[3] == ref_gene) {
        df.sub.group %>%
          magrittr::set_names(c("biorep", "Target", "Reference", "ddct1")) %>%
          dplyr::mutate(expression = 2^(-(Target - Reference - ddct1))) %>%
          dplyr::mutate(
            group = groups,
            gene = genes
          ) %>%
          dplyr::select(group, gene, biorep, expression) %>%
          rbind(res.all) -> res.all
      } else {
        df.sub.group %>%
          dplyr::select(1, 3, 2, 4) %>%
          magrittr::set_names(c("biorep", "Target", "Reference", "ddct1")) %>%
          dplyr::mutate(expression = 2^(-(Target - Reference - ddct1))) %>%
          dplyr::mutate(
            group = groups,
            gene = genes
          ) %>%
          dplyr::select(group, gene, biorep, expression) %>%
          rbind(res.all) -> res.all
      }
    }
  }

  # outlier removal
  if (remove_outliers) {
    res.all %>%
      dplyr::group_by(group, gene) %>%
      dplyr::mutate(is.out = find_outlier(expression)) %>%
      dplyr::ungroup() %>%
      dplyr::filter(is.out == "no") -> res.all
  }

  # summary statistics
  res.all %>%
    dplyr::group_by(group, gene) %>%
    dplyr::mutate(
      mean.expression = mean(expression),
      sd.expression = stats::sd(expression),
      n.biorep = dplyr::n(),
      se.expression = mean.expression / sqrt(n.biorep)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(temp = paste0(gene, group)) -> res.all

  # statistical tests
  res.all <- cal_stat_test(res.all, stat_method, ref_group)

  # plot
  df.plot <- res.all %>%
    dplyr::rename(
      Treatment = group,
      expre = expression,
      mean.expre = mean.expression,
      sd.expre = sd.expression,
      se.expre = se.expression,
      n = n.biorep
    )

  p <- build_exp_plot(df.plot, fig_type, fig_ncol)

  res <- list(table = df.plot, figure = p)
  return(res)
}

#' Detect outliers using IQR method
#' @keywords internal
find_outlier <- function(x) {
  ifelse(
    x < stats::quantile(x, 0.25) - 1.5 * stats::IQR(x) |
      x > stats::quantile(x, 0.75) + 1.5 * stats::IQR(x),
    "yes",
    "no"
  )
}

#' Run statistical tests on expression data
#' @keywords internal
cal_stat_test <- function(res.all, stat_method, ref_group) {
  if (stat_method == "t.test") {
    res.all %>%
      dplyr::group_by(gene) %>%
      rstatix::t_test(expression ~ group, ref.group = ref_group) %>%
      dplyr::ungroup() %>%
      dplyr::select(gene, group2, p) %>%
      dplyr::mutate(signif = dplyr::case_when(
        p < 0.001 ~ "***",
        p > 0.001 & p < 0.01 ~ "**",
        p > 0.01 & p < 0.05 ~ "*",
        TRUE ~ "NS"
      )) %>%
      dplyr::add_row(group2 = ref_group, p = NA, signif = NA) %>%
      dplyr::rename(group = group2) %>%
      dplyr::mutate(temp = paste0(gene, group)) %>%
      dplyr::select(temp, signif) -> df.stat

    res.all %>%
      dplyr::left_join(df.stat, by = "temp")
  } else if (stat_method == "wilcox.test") {
    res.all %>%
      dplyr::group_by(gene) %>%
      rstatix::wilcox_test(expression ~ group, ref.group = ref_group) %>%
      dplyr::ungroup() %>%
      dplyr::select(gene, group2, p) %>%
      dplyr::mutate(signif = dplyr::case_when(
        p < 0.001 ~ "***",
        p > 0.001 & p < 0.01 ~ "**",
        p > 0.01 & p < 0.05 ~ "*",
        TRUE ~ "NS"
      )) %>%
      dplyr::add_row(group2 = ref_group, p = NA, signif = NA) %>%
      dplyr::rename(group = group2) %>%
      dplyr::mutate(temp = paste0(gene, group)) %>%
      dplyr::select(temp, signif) -> df.stat

    res.all %>%
      dplyr::left_join(df.stat, by = "temp")
  } else {
    df.stat <- NULL
    for (i in unique(res.all$gene)) {
      df.sub <- res.all %>%
        dplyr::filter(gene == i) %>%
        dplyr::mutate(group = factor(group))
      fit <- stats::aov(expression ~ group, data = df.sub)
      tuk <- multcomp::glht(fit, linfct = multcomp::mcp(group = "Tukey"))
      multcomp::cld(tuk, level = 0.95, ddecreasing = TRUE)[["mcletters"]][["Letters"]] %>%
        as.data.frame() %>%
        dplyr::mutate(gene = i) %>%
        tibble::rownames_to_column(var = "group") %>%
        magrittr::set_colnames(c("group", "signif", "gene")) %>%
        dplyr::select(group, gene, signif) %>%
        dplyr::mutate(temp = paste0(group, gene)) %>%
        dplyr::select(temp, signif) %>%
        rbind(df.stat) -> df.stat
    }
    res.all %>%
      dplyr::mutate(temp = paste0(group, gene)) %>%
      dplyr::left_join(df.stat, by = "temp")
  }
}

#' Build expression plot (box or bar)
#' @keywords internal
build_exp_plot <- function(df.plot, fig_type, fig_ncol) {
  if (fig_type == "box") {
    df.plot %>%
      ggplot2::ggplot(ggplot2::aes(Treatment, expre, fill = Treatment)) +
      ggplot2::geom_boxplot(width = 0.6) +
      ggplot2::facet_wrap(. ~ gene, scales = "free_y", ncol = fig_ncol) +
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
  } else if (fig_type == "bar") {
    df.plot %>%
      dplyr::group_by(Treatment, gene) %>%
      dplyr::mutate(max.temp = max(expre)) %>%
      dplyr::ungroup() %>%
      ggplot2::ggplot(ggplot2::aes(Treatment, mean.expre / n, fill = Treatment)) +
      ggplot2::geom_bar(stat = "identity", width = 0.6) +
      ggplot2::geom_errorbar(ggplot2::aes(
        Treatment,
        ymin = mean.expre - sd.expre,
        ymax = mean.expre + sd.expre
      ), width = 0.2) +
      ggplot2::geom_jitter(ggplot2::aes(Treatment, expre), width = 0.1, alpha = 0.4) +
      ggplot2::geom_hline(ggplot2::aes(yintercept = max.temp * 1.1), color = NA) +
      ggplot2::facet_wrap(. ~ gene, scales = "free_y", ncol = fig_ncol) +
      ggplot2::geom_text(
        ggplot2::aes(Treatment, max.temp * 1.08, label = signif),
        check_overlap = TRUE, size = 4, color = "red"
      ) +
      ggthemes::theme_pander() +
      ggplot2::labs(y = "Relative expression") +
      ggplot2::theme(
        legend.position = "none",
        strip.text.x = ggplot2::element_text(face = "italic")
      )
  }
}
