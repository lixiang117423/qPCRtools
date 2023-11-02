#' @name CalExp2ddCt
#' @author Xiang LI <lixiang117423@@gmail.com>
#'
#' @title Calculate expression using standard curve.
#' @description Calculate expression using standard curve.
#'
#' @param cq.table The data frame of the position and cq value.
#' @param design.table The data frame of the position and corresponding information.
#' @param correction Correct expression value by reference gene.
#' @param ref.gene The name of reference gene.
#' @param ref.group The name of reference group.
#' @param stat.method Statistical method.
#' @param remove.outliers Remove the outliers of each group and gene, or not.
#' @param fig.type Output image type, `box` represents `boxplot`, `bar` represents `barplot`.
#' @param fig.ncol Number of columes of figure.
#'
#' @export
#' @return A list contain a table and a figure.
#' @examples
#' df1.path = system.file("examples", "ddct.cq.txt", package = "qPCRtools")
#' df2.path = system.file("examples", "ddct.design.txt", package = "qPCRtools")
#'
#' cq.table = read.table(df1.path, header = TRUE)
#' design.table = read.table(df2.path, header = TRUE)
#'
#' CalExp2ddCt(cq.table,
#'             design.table,
#'             ref.gene = "OsUBQ",
#'             ref.group = "CK",
#'             stat.method = "t.test",
#'             remove.outliers = TRUE,
#'             fig.type = "box",
#'             fig.ncol = NULL) -> res
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
  "remove.outliers",
  "ref.group",
  "fig.type",
  "fig.ncol",
  "out",
  "cq",
  "max.cq",
  "min.cq",
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
  'gene',
  'group',
  'biorep',
  'Target',
  'Reference',
  'ddct1',
  'mean.expression',
  'n.biorep',
  'sd.expression',
  'se.expression',
  'Reference',
  'Target',
  'biorep',
  'ddct1',
  'gene',
  'group',
  'mean.expression',
  'n.biorep',
  'sd.expression',
  'se.expression',
  'BioRep',
  'Eff',
  'Group',
  'TechRep',
  'rr.label',
  'p.value.label',
  'IQR',
  'quantile',
  'is.out'

))
CalExp2ddCt <- function(cq.table,
                        design.table,
                        ref.gene = "OsUBQ",
                        ref.group = "CK",
                        stat.method = "t.test",
                        remove.outliers = TRUE,
                        fig.type = "box",
                        fig.ncol = NULL) {

  # res
  res.all <- NULL

  # merge data
  cq.table %>%
    dplyr::left_join(design.table, by = "Position") %>%
    dplyr::rename(position = Position,
                  cq = Cq,
                  group = Group,
                  gene = Gene,
                  biorep = BioRep) -> df

  # for each gene
  target.genes <- setdiff(unique(df$gene), ref.gene)

  for (genes in target.genes) {
    df.sub <- df %>%
      dplyr::filter(gene %in% c(genes, ref.gene))

    df.sub.ck <- df.sub %>%
      dplyr::filter(group == ref.group)

    # reference gene in CK
    df.sub.ck.ref.gene <- df.sub %>%
      dplyr::filter(gene == ref.gene)
    mean.ck.ref.gene <- mean(df.sub.ck.ref.gene$cq)

    df.sub.ck.target.gene <- df.sub %>%
      dplyr::filter(gene != ref.gene)
    mean.ck.target.gene <- mean(df.sub.ck.target.gene$cq)

    dct1 <- mean.ck.target.gene - mean.ck.ref.gene

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
      # including ref gene or not
      if (ncol(df.sub.group) == 3 & !genes %in% colnames(df.sub.group)) {
        stop(paste0("Data of target gene ", genes, " has some problem, please check it and try again!"))
      }
      if (colnames(df.sub.group)[3] == ref.gene) {
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

  # find outliner function
  findoutliner <- function(x) {
    return(
      ifelse(
        x < quantile(x, .25) - 1.5 * IQR(x) | x > quantile(x, .75) + 1.5 * IQR(x),
        "yes",
        "no"
      )
    )
  }

  if (remove.outliers) {
    res.all %>%
      dplyr::group_by(group, gene) %>%
      dplyr::mutate(is.out = findoutliner(expression)) %>%
      dplyr::ungroup() %>%
      dplyr::filter(is.out == "no") -> res.all
  }else{
    res.all -> res.all
  }

  # group and mean and sd
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

  # statistics
  if (stat.method == "t.test") {
    res.all %>%
      dplyr::group_by(gene) %>%
      rstatix::t_test(expression ~ group, ref.group = ref.group) %>%
      dplyr::ungroup() %>%
      dplyr::select(gene, group2, p) %>%
      dplyr::mutate(signif = dplyr::case_when(
        p < 0.001 ~ "***",
        p > 0.001 & p < 0.01 ~ "**",
        p > 0.01 & p < 0.05 ~ "*",
        TRUE ~ "NS"
      )) %>%
      dplyr::add_row(group2 = ref.group, p = NA, signif = NA) %>%
      dplyr::rename(group = group2) %>%
      dplyr::mutate(temp = paste0(gene, group)) %>%
      dplyr::select(temp, signif) -> df.stat

    res.all <- res.all %>%
      dplyr::left_join(df.stat, by = "temp")
  } else if (stat.method == "wilcox.test") {
    res.all %>%
      dplyr::group_by(gene) %>%
      rstatix::wilcox_test(expression ~ group, ref.group = ref.group) %>%
      dplyr::ungroup() %>%
      dplyr::select(gene, group2, p) %>%
      dplyr::mutate(signif = dplyr::case_when(
        p < 0.001 ~ "***",
        p > 0.001 & p < 0.01 ~ "**",
        p > 0.01 & p < 0.05 ~ "*",
        TRUE ~ "NS"
      )) %>%
      dplyr::add_row(group2 = ref.group, p = NA, signif = NA) %>%
      dplyr::rename(group = group2) %>%
      dplyr::mutate(temp = paste0(gene, group)) %>%
      dplyr::select(temp, signif) -> df.stat

    res.all <- res.all %>%
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
      dplyr::left_join(df.stat, by = "temp") -> res.all
  }

  # plot
  df.plot <- res.all %>%
    dplyr::rename(
      Treatment = group,
      gene = gene,
      expre = expression,
      mean.expre = mean.expression,
      sd.expre = sd.expression,
      se.expre = se.expression,
      n = n.biorep
    )


  if (fig.type == "box") {
    df.plot %>%
      ggplot2::ggplot(ggplot2::aes(Treatment, expre, fill = Treatment)) +
      ggplot2::geom_boxplot(width = 0.6) +
      ggplot2::facet_wrap(. ~ gene, scales = "free_y", ncol = fig.ncol) +
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
      dplyr::group_by(Treatment, gene) %>%
      dplyr::mutate(max.temp = max(expre)) %>%
      dplyr::ungroup() %>%
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
      ggplot2::facet_wrap(. ~ gene, scales = "free_y", ncol = fig.ncol) +
      ggplot2::geom_text(ggplot2::aes(Treatment, max.temp * 1.08, label = signif),
                         check_overlap = TRUE, size = 4, color = "red"
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
