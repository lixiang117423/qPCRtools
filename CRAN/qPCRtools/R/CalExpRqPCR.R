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
  "gene",
  "group",
  "biorep",
  "Target",
  "Reference",
  "ddct1",
  "mean.expression",
  "n.biorep",
  "sd.expression",
  "se.expression",
  "Reference",
  "Target",
  "biorep",
  "ddct1",
  "gene",
  "group",
  "mean.expression",
  "n.biorep",
  "sd.expression",
  "se.expression",
  "RawPValue",
  "PValue",
  "PAdjustMethod",
  "contrast",
  "letters",
  'sd.cq',
  'eff',
  'min.mean.cq',
  'QCq',
  'techrep',
  'temp_2',
  'SD_QCq',
  'SD.factor',
  'SD_1',
  'min.expression',
  'Expre4Stat',
  'Expression',
  'SD',
  'SE',
  'BioRep',
  'Eff',
  'Group',
  'TechRep'
))

#' @name CalExpRqPCR
#' @author Xiang LI <lixiang117423@@gmail.com>
#'
#' @title Calculate expression using standard curve.
#' @description Calculate expression using standard curve.
#'
#' @param cq.table The data frame of the position and cq value.
#' @param design.table The data frame of the position and corresponding information.
#' @param ref.gene The name of reference gene.
#' @param ref.group The name of reference group.
#' @param stat.method Statistical method.
#' @param fig.type Output image type, `box` represents `boxplot`, `bar` represents `barplot`.
#' @param fig.ncol Number of columes of figure.
#'
#' @export
#'
#' @return A list contain a table and a figure.
#'
#' @examples
#' df1.path <- system.file("examples", "cal.expre.rqpcr.cq.txt", package = "qPCRtools")
#' df2.path <- system.file("examples", "cal.expre.rqpcr.design.txt", package = "qPCRtools")
#'
#' cq.table <- read.table(df1.path, header = TRUE)
#' design.table <- read.table(df2.path, header = TRUE)
#'
#' CalExpRqPCR(cq.table,
#'            design.table,
#'            ref.gene = NULL,
#'            ref.group = "CK",
#'            stat.method = "t.test",
#'            fig.type = "box",
#'            fig.ncol = NULL
#'            ) -> res
#'
#' res[["table"]]
#' res[["figure"]]
#'
CalExpRqPCR <- function(cq.table,
                        design.table,
                        ref.gene = NULL,
                        ref.group = "CK",
                        stat.method = "t.test",
                        fig.type = "box",
                        fig.ncol = NULL) {

  # merge data
  cq.table %>%
    dplyr::left_join(design.table, by = "Position") %>%
    dplyr::rename(position = Position,
                  cq = Cq,
                  group = Group,
                  gene = Gene,
                  biorep = BioRep,
                  techrep = TechRep,
                  eff = Eff) -> df

  # start cal
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

  # 几何均值（参考基因 QCq 的归一化因子用）。2026-10-15：原来只在"自动选内参"分支里定义，
  # 导致显式给 ref.gene 时在本函数后半段报错 "没有 geometric.mean 这个函数"；现改为两种路径共用。
  geometric.mean <- function(x) {
    x <- x[!is.na(x)]
    if (any(x < 0)) {
      stop("'x' contains negative value(s)")
    } else {
      return(prod(x)^(1 / length(x)))
    }
  }

  # find ref gene
  if (!is.null(ref.gene)) {
    ref.gene <- ref.gene
  } else {
    df.ref <- df.expre %>%
      dplyr::select(group, gene, cq, biorep, techrep) %>%
      # dplyr::group_by(Treatment, Gene, bio_rep, rep) %>%
      dplyr::mutate(
        Treatment = paste0(group, biorep, techrep)
      ) %>%
      tidyr::spread(key = gene, value = cq)

    df.temp <- df.ref[, 5:ncol(df.ref)] %>% as.data.frame()

    n <- length(unique(df.expre$gene))

    # 注：这里原来还有一段算 M 的代码，但下面的逐步剔除循环每轮都会重算 M（gene.stable），
    # 那段结果从未被使用；为避免"同一件事两套实现"，2026-10-15 删除（稳定性度量只在 gene.stable 里一处）。
    geneSymbol <- colnames(df.temp)
    n <- ncol(df.temp)

    num.ref <- 2

    V <- numeric(n - num.ref)
    names(V) <- paste(((n - 1):num.ref), "/", (n:(num.ref + 1)), sep = "")
    meanM <- numeric(n - num.ref + 1)
    names(meanM) <- as.character(n:num.ref)
    R <- character(n)
    names(R) <- as.character(c(rep(1, num.ref), (num.ref + 1):length(R)))

    gene.stable <- function(data, na.rm = TRUE) {
      if (!is.data.frame(data) & !is.matrix(data)) {
        stop("'data' has to of class matrix or data.frame")
      }
      n <- ncol(data)
      if (n == 1) stop("you need at least two genes for this computation")
      # geNorm 的稳定性 M 必须建立在**相对量**上：Q = 2^(minCq − Cq)（每个基因相对自身表达
      # 最高的样品）。原实现直接用原始 Cq 相除 log2(Cq_j / Cq_k)，与 geNorm 原文
      # （Vandesompele et al. 2002）不符 —— 2026-10-15 改正（导师审稿意见第 7 条）。
      # 注意：必须用 sweep() 按列减去各基因的最小 Cq；写 data - apply(data, 2, min) 会因
      # R 的列主序向量回收（周期 = 列数 ≠ 行数）而错位。
      column.min <- apply(data, 2, min, na.rm = TRUE)
      data <- 2^(-sweep(data, 2, column.min, "-"))
      M <- numeric(n)
      for (j in 1:n) {
        A <- log2(data[, j] / data[, -j])
        if (n > 2) {
          M[j] <- mean(apply(A, 2, sd, na.rm = na.rm))
        } else {
          M[j] <- sd(A, na.rm = na.rm)
        }
      }
      if (is.data.frame(data)) {
        names(M) <- names(data)
      } else {
        names(M) <- colnames(data)
      }
      return(M)
    }


    for (i in n:num.ref) {
      M <- gene.stable(df.temp, na.rm = TRUE)
      ind <- which.max(M)
      meanM[n - i + 1] <- mean(M)
      if (i == num.ref) {
        R[1:num.ref] <- geneSymbol
      } else {
        R[i] <- geneSymbol[ind]
      }
      if (i > 2) {
        NF.old <- apply(df.temp, 1, geometric.mean)
        NF.new <- apply(df.temp[, -ind], 1, geometric.mean)
        V[n - i + 1] <- sd(log2(NF.new / NF.old), na.rm = TRUE)
      }
      df.temp <- df.temp[, -ind]
      geneSymbol <- geneSymbol[-ind]
    }

    ref.gene <- as.character(R[1:num.ref])
  }

  #### cal factor
  df.factor <- df.expre %>%
    dplyr::filter(gene %in% ref.gene) %>%
    dplyr::mutate(temp_2 = paste0(group, biorep, gene))
  df.factor <- df.factor[!duplicated(df.factor$temp_2), ] %>% as.data.frame()

  factor <- data.frame()

  for (i in unique(df.factor$biorep)) {
    df.temp <- df.factor %>% dplyr::filter(biorep == i)
    for (j in unique(df.temp$group)) {
      df.temp.2 <- df.temp %>%
        dplyr::filter(group == j) %>%
        dplyr::select(group, QCq)
      fac <- data.frame(group = j, biorep = i, factor = geometric.mean(df.temp.2$QCq))
      factor <- rbind(factor, fac)
    }
  }

  factor <- factor %>%
    dplyr::mutate(temp_2 = paste0(group, biorep)) %>%
    dplyr::select(temp_2, factor)
  df.factor <- df.factor %>%
    dplyr::mutate(temp_2 = paste0(group, biorep)) %>%
    merge(factor, by = "temp_2") %>%
    dplyr::mutate(SD.factor = (SD_QCq / (length(ref.gene) * (QCq)))^2) %>%
    dplyr::group_by(biorep, group) %>%
    dplyr::mutate(SD.factor = sqrt(sum(SD.factor)) * factor)

  ###################################
  #### cal corred expresssion              ####
  ###################################
  df.goi <- df.expre %>%
    dplyr::filter(!gene %in% ref.gene) %>%
    dplyr::mutate(temp_2 = paste0(group, biorep)) %>%
    merge(df.factor[, c("temp_2", "factor", "SD.factor")], by = "temp_2") %>%
    dplyr::mutate(
      expression = QCq / factor,
      SD_1 = expression * sqrt((SD_QCq / QCq)^2 + (SD.factor / factor)^2),
      SE_1 = SD_1 / sqrt(2)
    )

  ###################################
  #### cal mean expression
  res.all <- df.goi %>%
    dplyr::ungroup() %>%
    # dplyr::mutate(Treatment = stringr::str_sub(Sample, 1, nchar(Sample) - 2)) %>%
    dplyr::group_by(gene, group) %>%
    dplyr::mutate(
      # 组内汇总按**样本身份**计：数值相同的生物学重复仍是独立观测，不能按数值去重
      # （原写 mean(unique(expression)) / sd(unique(expression))，会把 1、1、4 算成 2.5；
      #  2026-10-15 按方法原文 Hilliou & Tran 2013 §2.2.4 改正）
      mean.expression = mean(expression, na.rm = TRUE),
      sd.expression = sd(expression, na.rm = TRUE),
      se.expression = sd.expression / sqrt(dplyr::n())
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

  # statistics
  if (stat.method == "t.test") {
    res.all %>%
      dplyr::group_by(gene) %>%
      rstatix::t_test(Expre4Stat ~ group, ref.group = ref.group) %>%
      dplyr::ungroup() %>%
      dplyr::select(gene, group2, p) %>%
      dplyr::mutate(signif = dplyr::case_when(
        p < 0.001 ~ "***",
        p < 0.01 ~ "**",
        p < 0.05 ~ "*",
        TRUE ~ "NS"
      ), RawPValue = p, PValue = p, PAdjustMethod = "none") %>%
      dplyr::add_row(
        group2 = ref.group, p = NA, signif = NA,
        RawPValue = NA, PValue = NA, PAdjustMethod = NA
      ) %>%
      dplyr::rename(group = group2) %>%
      dplyr::mutate(temp = paste0(gene, group)) %>%
      dplyr::select(temp, signif, RawPValue, PValue, PAdjustMethod) -> df.stat

    res.all <- res.all %>%
      dplyr::left_join(df.stat, by = "temp")
  } else if (stat.method == "wilcox.test") {
    res.all %>%
      dplyr::group_by(gene) %>%
      rstatix::wilcox_test(Expre4Stat ~ group, ref.group = ref.group) %>%
      dplyr::ungroup() %>%
      dplyr::select(gene, group2, p) %>%
      dplyr::mutate(signif = dplyr::case_when(
        p < 0.001 ~ "***",
        p < 0.01 ~ "**",
        p < 0.05 ~ "*",
        TRUE ~ "NS"
      ), RawPValue = p, PValue = p, PAdjustMethod = "none") %>%
      dplyr::add_row(
        group2 = ref.group, p = NA, signif = NA,
        RawPValue = NA, PValue = NA, PAdjustMethod = NA
      ) %>%
      dplyr::rename(group = group2) %>%
      dplyr::mutate(temp = paste0(gene, group)) %>%
      dplyr::select(temp, signif, RawPValue, PValue, PAdjustMethod) -> df.stat

    res.all <- res.all %>%
      dplyr::left_join(df.stat, by = "temp")
  } else {
    df.stat <- NULL
    for (i in unique(res.all$gene)) {
      df.sub <- res.all %>%
        dplyr::filter(gene == i) %>%
        dplyr::mutate(group = factor(group))
      fit <- stats::aov(Expre4Stat ~ group, data = df.sub)
      tuk <- multcomp::glht(fit, linfct = multcomp::mcp(group = "Tukey"))
      # cld.glht 的 level 是显著性水平（默认 0.05）；原写成 0.95 = α=0.95，几乎全标成不同字母
      # （2026-10-15 改正；decreasing 原拼写 ddecreasing 被忽略）
      letters <- multcomp::cld(
        tuk, level = 0.05, decreasing = TRUE
      )[["mcletters"]][["Letters"]] %>%
        as.data.frame() %>%
        dplyr::mutate(gene = i) %>%
        tibble::rownames_to_column(var = "group") %>%
        magrittr::set_colnames(c("group", "signif", "gene"))

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
        dplyr::mutate(group = dplyr::case_when(
          endsWith(contrast, paste0(" - ", ref.group)) ~
            sub(paste0(" - ", ref.group, "$"), "", contrast),
          startsWith(contrast, paste0(ref.group, " - ")) ~
            sub(paste0("^", ref.group, " - "), "", contrast),
          TRUE ~ NA_character_
        )) %>%
        dplyr::filter(!is.na(group)) %>%
        dplyr::mutate(PAdjustMethod = "Tukey") %>%
        dplyr::select(group, RawPValue, PValue, PAdjustMethod)

      letters %>%
        dplyr::left_join(p.values, by = "group") %>%
        dplyr::mutate(temp = paste0(group, gene)) %>%
        dplyr::select(temp, signif, RawPValue, PValue, PAdjustMethod) %>%
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
      expre = Expre4Stat,
      mean.expre = Expression,
      sd.expre = SD,
      se.expre = SE
    ) %>%
    dplyr::group_by(gene, Treatment) %>%
    dplyr::mutate(n = dplyr::n()) %>%
    dplyr::ungroup()

  if (fig.type == "box") {
    df.plot %>%
      ggplot2::ggplot(ggplot2::aes(Treatment, expre, fill = Treatment)) +
      ggplot2::geom_boxplot(width = 0.6) +
      ggplot2::facet_wrap(. ~ gene, scales = "free_y", ncol = fig.ncol) +
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
      dplyr::distinct(Treatment, gene, .keep_all = TRUE) %>%
      dplyr::group_by(gene) %>%
      dplyr::mutate(max.temp = max(mean.expre)) %>%
      dplyr::ungroup() %>%
      ggplot2::ggplot(ggplot2::aes(Treatment, mean.expre, fill = Treatment)) +
      ggplot2::geom_bar(stat = "identity", width = 0.6) +
      ggplot2::geom_errorbar(ggplot2::aes(Treatment,
        ymin = mean.expre - sd.expre,
        ymax = mean.expre + sd.expre
      ),
      width = 0.2
      ) +
      # ggplot2::geom_jitter(ggplot2::aes(Treatment, expre), width = 0.1, alpha = 0.4) +
      ggplot2::geom_hline(ggplot2::aes(yintercept = max.temp * 1.15), color = NA) +
      ggplot2::facet_wrap(. ~ gene, scales = "free_y", ncol = fig.ncol) +
      ggplot2::geom_text(ggplot2::aes(Treatment, (mean.expre + sd.expre) * 1.08, label = signif),
        check_overlap = TRUE, size = 4, color = "black"
      ) +
      ggthemes::theme_pander() +
      ggplot2::labs(y = "Relative expression") +
      ggplot2::theme(
        legend.position = "none",
        strip.text.x = ggplot2::element_text(face = "italic")
      ) -> p
  }
  res.all <- res.all %>%
    dplyr::select(-temp, -eff) %>%
    dplyr::select(
      group, biorep, gene, Expre4Stat, Expression, SD, SE,
      signif, RawPValue, PValue, PAdjustMethod
    )
  res <- list(table = res.all, figure = p)
  return(res)
}
