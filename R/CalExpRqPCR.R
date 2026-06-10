globalVariables(c(
  "Position", "Cq", "Group", "Gene", "BioRep", "TechRep", "Eff",
  "cq", "group", "gene", "biorep", "techrep", "eff",
  "mean.cq", "sd.cq", "min.mean.cq", "QCq", "SD_QCq",
  "expression", "mean.expression", "sd.expression", "se.expression",
  "min.expression", "Expre4Stat", "Expression", "SD", "SE",
  "signif", "temp", "group2", "max.temp", "n", "p",
  "Treatment", "expre", "mean.expre", "sd.expre",
  "temp_2", "factor", "SD.factor", "SD_1", "SE_1"
))

#' Calculate Expression Using RqPCR Method
#'
#' Calculate relative gene expression using the RqPCR method with
#' amplification efficiency correction. Can auto-select reference genes
#' using the GeNorm algorithm when ref_gene is NULL.
#'
#' @param cq_table A data frame containing position and Cq values.
#'   Must have columns: Position, Gene, Cq, BioRep, TechRep, Eff.
#' @param design_table A data frame containing position and group information.
#'   Must have columns: Position, Group, BioRep, TechRep, Eff.
#' @param ref_gene Character. The name(s) of reference gene(s).
#'   If NULL, reference genes are auto-selected via GeNorm (default: NULL).
#' @param ref_group Character. The name of the reference/control group
#'   (default: "CK").
#' @param stat_method Character. Statistical method for group comparison.
#'   One of "t.test", "wilcox.test", or "anova" (default: "t.test").
#' @param fig_type Character. Plot type: "box" for boxplot, "bar" for barplot
#'   (default: "box").
#' @param fig_ncol Integer. Number of columns in facet plot (default: NULL).
#'
#' @return A list containing:
#'   \item{table}{Data frame with expression values and statistics}
#'   \item{figure}{ggplot object}
#'
#' @importFrom magrittr %>%
#' @importFrom stats sd aov
#'
#' @export
#'
#' @examples
#' \dontrun{
#' df1.path <- system.file("examples", "cal.expre.rqpcr.cq.txt", package = "qPCRtools")
#' df2.path <- system.file("examples", "cal.expre.rqpcr.design.txt", package = "qPCRtools")
#' cq_table <- read.table(df1.path, header = TRUE)
#' design_table <- read.table(df2.path, header = TRUE)
#' res <- CalExpRqPCR(
#'   cq_table, design_table,
#'   ref_gene = NULL,
#'   ref_group = "CK",
#'   stat_method = "t.test",
#'   fig_type = "box",
#'   fig_ncol = NULL
#' )
#' res[["table"]]
#' res[["figure"]]
#' }
#'
#' @author Xiang LI <lixiang117423@gmail.com>
CalExpRqPCR <- function(cq_table,
                        design_table,
                        ref_gene = NULL,
                        ref_group = "CK",
                        stat_method = "t.test",
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

  # merge data
  cq_table %>%
    dplyr::left_join(design_table, by = "Position") %>%
    dplyr::rename(
      position = Position,
      cq = Cq,
      group = Group,
      gene = Gene,
      biorep = BioRep,
      techrep = TechRep,
      eff = Eff
    ) -> df

  # calculate QCq
  df.expre <- df %>%
    dplyr::group_by(biorep, group, gene) %>%
    dplyr::mutate(
      mean.cq = mean(cq, na.rm = TRUE),
      sd.cq = sd(cq, na.rm = TRUE),
      sd.cq = ifelse(is.na(sd.cq), 0, sd.cq)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(biorep, gene) %>%
    dplyr::mutate(
      min.mean.cq = min(mean.cq),
      QCq = eff^(min.mean.cq - mean.cq),
      SD_QCq = sd.cq * QCq * log(eff)
    ) %>%
    dplyr::ungroup()

  # auto-select reference genes if not provided
  if (is.null(ref_gene)) {
    ref_gene <- find_ref_gene(df.expre)
  }

  # calculate normalization factor
  df.factor <- df.expre %>%
    dplyr::filter(gene %in% ref_gene) %>%
    dplyr::mutate(temp_2 = paste0(group, biorep, gene))
  df.factor <- df.factor[!duplicated(df.factor$temp_2), ] %>% as.data.frame()

  norm_factor <- data.frame()
  for (i in unique(df.factor$biorep)) {
    df.temp <- df.factor %>% dplyr::filter(biorep == i)
    for (j in unique(df.temp$group)) {
      df.temp.2 <- df.temp %>%
        dplyr::filter(group == j) %>%
        dplyr::select(group, QCq)
      fac <- data.frame(group = j, biorep = i, factor = geometric_mean(df.temp.2$QCq))
      norm_factor <- rbind(norm_factor, fac)
    }
  }

  norm_factor <- norm_factor %>%
    dplyr::mutate(temp_2 = paste0(group, biorep)) %>%
    dplyr::select(temp_2, factor)
  df.factor <- df.factor %>%
    dplyr::mutate(temp_2 = paste0(group, biorep)) %>%
    merge(norm_factor, by = "temp_2") %>%
    dplyr::mutate(SD.factor = (SD_QCq / (length(ref_gene) * QCq))^2) %>%
    dplyr::group_by(biorep, group) %>%
    dplyr::mutate(SD.factor = sqrt(sum(SD.factor)) * factor)

  # corrected expression
  df.goi <- df.expre %>%
    dplyr::filter(!gene %in% ref_gene) %>%
    dplyr::mutate(temp_2 = paste0(group, biorep)) %>%
    merge(df.factor[, c("temp_2", "factor", "SD.factor")], by = "temp_2") %>%
    dplyr::mutate(
      expression = QCq / factor,
      SD_1 = expression * sqrt((SD_QCq / QCq)^2 + (SD.factor / factor)^2),
      SE_1 = SD_1 / sqrt(2)
    )

  # mean expression
  res.all <- df.goi %>%
    dplyr::ungroup() %>%
    dplyr::group_by(gene, group) %>%
    dplyr::mutate(
      mean.expression = mean(unique(expression), na.rm = TRUE),
      sd.expression = sd(unique(expression), na.rm = TRUE),
      se.expression = sd.expression / sqrt(length(unique(biorep)))
    ) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(gene) %>%
    dplyr::mutate(
      min.expression = min(mean.expression),
      mean.expression = mean.expression / min.expression,
      sd.expression = sd.expression / min.expression,
      se.expression = se.expression / min.expression
    ) %>%
    dplyr::select(
      group, gene, eff, expression, biorep,
      mean.expression, sd.expression, se.expression
    ) %>%
    dplyr::rename(
      Expre4Stat = expression,
      Expression = mean.expression,
      SD = sd.expression,
      SE = se.expression
    ) %>%
    dplyr::mutate(temp = paste0(group, gene, biorep)) %>%
    dplyr::filter(!duplicated(temp)) %>%
    dplyr::select(-temp) %>%
    dplyr::mutate(temp = paste0(gene, group))

  # statistical tests
  res.all <- cal_rqpcr_stat_test(res.all, stat_method, ref_group)

  # plot
  df.plot <- res.all %>%
    dplyr::rename(
      Treatment = group,
      expre = Expre4Stat,
      mean.expre = Expression,
      sd.expre = SD,
      se.expre = SE
    ) %>%
    dplyr::group_by(gene, Treatment) %>%
    dplyr::mutate(n = dplyr::n()) %>%
    dplyr::ungroup()

  p <- build_rqpcr_plot(df.plot, fig_type, fig_ncol)

  res.all <- res.all %>%
    dplyr::select(-temp, -eff) %>%
    dplyr::select(1, 4, 2, 3, 5:7)

  res <- list(table = res.all, figure = p)
  return(res)
}

#' Calculate geometric mean
#' @keywords internal
geometric_mean <- function(x) {
  x <- x[!is.na(x)]
  if (any(x < 0)) {
    stop("'x' contains negative value(s)")
  }
  prod(x)^(1 / length(x))
}

#' Find reference genes using GeNorm algorithm
#' @keywords internal
find_ref_gene <- function(df.expre) {
  df.ref <- df.expre %>%
    dplyr::select(group, gene, cq, biorep, techrep) %>%
    dplyr::mutate(Treatment = paste0(group, biorep, techrep)) %>%
    tidyr::spread(key = gene, value = cq)

  df.temp <- df.ref[, 5:ncol(df.ref)] %>% as.data.frame()
  n <- length(unique(df.expre$gene))

  gene.stable <- function(data, na.rm = TRUE) {
    if (!is.data.frame(data) & !is.matrix(data)) {
      stop("'data' has to be of class matrix or data.frame")
    }
    n.cols <- ncol(data)
    if (n.cols == 1) stop("you need at least two genes for this computation")
    M <- numeric(n.cols)
    for (j in 1:n.cols) {
      A <- log2(data[, j] / data[, -j])
      if (n.cols > 2) {
        M[j] <- mean(apply(A, 2, stats::sd, na.rm = na.rm))
      } else {
        M[j] <- stats::sd(A, na.rm = na.rm)
      }
    }
    if (is.data.frame(data)) {
      names(M) <- names(data)
    } else {
      names(M) <- colnames(data)
    }
    M
  }

  geneSymbol <- colnames(df.temp)
  n.genes <- ncol(df.temp)
  num.ref <- 2
  R <- character(n.genes)
  names(R) <- as.character(c(rep(1, num.ref), (num.ref + 1):length(R)))

  for (i in n.genes:num.ref) {
    M <- gene.stable(df.temp, na.rm = TRUE)
    ind <- which.max(M)
    if (i == num.ref) {
      R[1:num.ref] <- geneSymbol
    } else {
      R[i] <- geneSymbol[ind]
    }
    df.temp <- df.temp[, -ind]
    geneSymbol <- geneSymbol[-ind]
  }

  as.character(R[1:num.ref])
}

#' Run statistical tests for RqPCR expression
#' @keywords internal
cal_rqpcr_stat_test <- function(res.all, stat_method, ref_group) {
  if (stat_method == "t.test") {
    res.all %>%
      dplyr::group_by(gene) %>%
      rstatix::t_test(Expre4Stat ~ group, ref.group = ref_group) %>%
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
      rstatix::wilcox_test(Expre4Stat ~ group, ref.group = ref_group) %>%
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
      fit <- stats::aov(Expre4Stat ~ group, data = df.sub)
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

#' Build RqPCR expression plot
#' @keywords internal
build_rqpcr_plot <- function(df.plot, fig_type, fig_ncol) {
  if (fig_type == "box") {
    df.plot %>%
      ggplot2::ggplot(ggplot2::aes(Treatment, expre, fill = Treatment)) +
      ggplot2::geom_boxplot(width = 0.6) +
      ggplot2::facet_wrap(. ~ gene, scales = "free_y", ncol = fig_ncol) +
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
      dplyr::group_by(gene) %>%
      dplyr::mutate(max.temp = max(mean.expre)) %>%
      dplyr::ungroup() %>%
      ggplot2::ggplot(ggplot2::aes(Treatment, mean.expre / n, fill = Treatment)) +
      ggplot2::geom_bar(stat = "identity", width = 0.6) +
      ggplot2::geom_errorbar(ggplot2::aes(
        Treatment,
        ymin = mean.expre - sd.expre,
        ymax = mean.expre + sd.expre
      ), width = 0.2) +
      ggplot2::geom_hline(ggplot2::aes(yintercept = max.temp * 1.15), color = NA) +
      ggplot2::facet_wrap(. ~ gene, scales = "free_y", ncol = fig_ncol) +
      ggplot2::geom_text(
        ggplot2::aes(Treatment, (mean.expre + sd.expre) * 1.08, label = signif),
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
