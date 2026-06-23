#' 物种丰度比较分析函数（最终版 - 添加显著性标记）
#'
#' @param abundance_matrix 物种丰度矩阵，行为样本，列为物种
#' @param grouping_factor 样本分组因子，长度等于abundance_matrix的行数
#' @param species_name 要分析的物种名称（必须是abundance_matrix的列名之一）
#' @param plot_title 图表标题，默认为"物种丰度比较"
#' @param ylab y轴标签，默认为"相对丰度"
#' @param xlab x轴标签，默认为"分组"
#' @param plot_type 可视化类型，可选"boxplot"、"barplot"、"violin"、"dotplot"，默认为"boxplot"
#' @param add_points 是否添加散点，默认为TRUE
#' @param point_type 散点类型，可选"jitter"、"beeswarm"、"center"，默认为"jitter"
#' @param point_shape 散点形状，默认为21
#' @param point_size 散点大小，默认为2.5
#' @param point_alpha 散点透明度，默认为0.8
#' @param point_color 散点边框颜色，默认为"white"
#' @param point_stroke 散点边框宽度，默认为0.5
#' @param filled 主图形是否填充颜色，默认为TRUE
#' @param border_width 当filled=FALSE时，边框宽度，默认为1.5
#' @param main_alpha 主图形透明度，默认为0.7
#' @param test_method 差异性检验方法，可选"wilcox"或"t.test"，默认为"wilcox"
#' @param paired 是否为配对样本检验，默认为FALSE
#' @param conf.level 检验的置信水平，默认为0.95
#' @param p_adjust_method 多重检验校正方法，默认为"none"
#' @param show_pvalue 是否在caption中显示p值和效应量，默认为TRUE
#' @param add_sig_label 是否在图上添加显著性标记，可选"none"（不添加）、"pvalue"（显示p值）、
#'                      "asterisk"（显示星号），默认为"asterisk"
#' @param sig_label_size 显著性标记的字体大小，默认为10
#' @param sig_hide_ns 是否隐藏非显著性标记，默认为TRUE
#' @param colors 自定义配色方案，默认为NULL
#' @param effect_size 是否计算效应量，默认为TRUE
#'
#' @return 返回一个列表，包含plot、test_result、summary_stats等
#'
#' @examples
#' # 示例1：添加星号显著性标记
#' result <- compare_species_abundance_v2(
#'   abundance, groups, "species1",
#'   add_sig_label = "asterisk"
#' )
#'
#' # 示例2：添加p值标记
#' result <- compare_species_abundance_v2(
#'   abundance, groups, "species1",
#'   add_sig_label = "pvalue"
#' )
#'
#' @importFrom stats wilcox.test t.test kruskal.test oneway.test p.adjust
#' @importFrom dplyr group_by summarise n
#' @importFrom tidyplots tidyplot
#' @export
compare_species_abundance_v2 <- function(
  abundance_matrix,
  grouping_factor,
  species_name,
  plot_title = "物种丰度比较",
  ylab = "相对丰度",
  xlab = "分组",
  plot_type = c("boxplot", "barplot", "violin", "dotplot"),
  add_points = TRUE,
  point_type = c("jitter", "beeswarm", "center"),
  point_shape = 21,
  point_size = 2.5,
  point_alpha = 0.8,
  point_color = "white",
  point_stroke = 0.5,
  filled = TRUE,
  border_width = 1.5,
  main_alpha = 0.7,
  test_method = c("wilcox", "t.test"),
  paired = FALSE,
  conf.level = 0.95,
  p_adjust_method = c(
    "none",
    "bonferroni",
    "holm",
    "hochberg",
    "hommel",
    "BH",
    "BY",
    "fdr"
  ),
  show_pvalue = TRUE,
  add_sig_label = c("asterisk", "pvalue", "none"),
  sig_label_size = 10,
  sig_hide_ns = TRUE,
  colors = NULL,
  effect_size = TRUE
) {
  # ============================================================================
  # 1. 加载必要的包
  # ============================================================================
  required_packages <- c("tidyplots", "dplyr", "showtext")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste0("请安装'", pkg, "'包: install.packages('", pkg, "')"))
    }
  }

  library(showtext)
  showtext_auto()
  library(tidyplots)
  library(dplyr)

  # ============================================================================
  # 2. 参数验证和匹配
  # ============================================================================
  plot_type <- match.arg(plot_type)
  test_method <- match.arg(test_method)
  p_adjust_method <- match.arg(p_adjust_method)
  point_type <- match.arg(point_type)
  add_sig_label <- match.arg(add_sig_label)

  # 检查输入
  if (!species_name %in% colnames(abundance_matrix)) {
    stop(paste0("物种'", species_name, "'不在丰度矩阵的列名中"))
  }

  if (length(grouping_factor) != nrow(abundance_matrix)) {
    stop("分组因子的长度必须等于丰度矩阵的行数")
  }

  # ============================================================================
  # 3. 数据准备
  # ============================================================================
  species_abundance <- abundance_matrix[[species_name]]

  data <- data.frame(
    Abundance = species_abundance,
    Group = factor(grouping_factor)
  )

  # ============================================================================
  # 4. 描述性统计
  # ============================================================================
  summary_stats <- data %>%
    group_by(Group) %>%
    summarise(
      n = n(),
      mean = mean(Abundance, na.rm = TRUE),
      median = median(Abundance, na.rm = TRUE),
      sd = sd(Abundance, na.rm = TRUE),
      se = sd / sqrt(n),
      min = min(Abundance, na.rm = TRUE),
      max = max(Abundance, na.rm = TRUE),
      q25 = quantile(Abundance, 0.25, na.rm = TRUE),
      q75 = quantile(Abundance, 0.75, na.rm = TRUE),
      .groups = "drop"
    )

  # ============================================================================
  # 5. 统计检验
  # ============================================================================
  n_groups <- length(unique(data$Group))

  if (n_groups < 2) {
    stop("需要至少两组数据进行比较")
  }

  # 执行统计检验
  if (n_groups == 2) {
    group_levels <- levels(data$Group)
    group1_data <- data$Abundance[data$Group == group_levels[1]]
    group2_data <- data$Abundance[data$Group == group_levels[2]]

    if (test_method == "wilcox") {
      test_result <- wilcox.test(
        group1_data,
        group2_data,
        paired = paired,
        conf.level = conf.level,
        exact = FALSE
      )
    } else {
      test_result <- t.test(
        group1_data,
        group2_data,
        paired = paired,
        conf.level = conf.level
      )
    }

    raw_pvalue <- test_result$p.value

    # 计算效应量
    effect_size_result <- NULL
    if (effect_size) {
      effect_size_result <- calculate_effect_size(
        group1_data,
        group2_data,
        test_method,
        paired
      )
    }
  } else {
    # 多组比较
    if (test_method == "wilcox") {
      test_result <- kruskal.test(Abundance ~ Group, data = data)
    } else {
      test_result <- oneway.test(
        Abundance ~ Group,
        data = data,
        var.equal = TRUE
      )
    }

    raw_pvalue <- test_result$p.value

    # 多组比较的效应量
    effect_size_result <- NULL
    if (effect_size) {
      if (test_method == "wilcox") {
        effect_size_result <- calculate_kruskal_effect_size(
          test_result,
          nrow(data)
        )
      } else {
        effect_size_result <- calculate_anova_effect_size(data)
      }
    }
  }

  # ============================================================================
  # 6. 多重检验校正
  # ============================================================================
  adjusted_pvalue <- NULL
  if (p_adjust_method != "none") {
    adjusted_pvalue <- p.adjust(raw_pvalue, method = p_adjust_method)
    display_pvalue <- adjusted_pvalue
  } else {
    display_pvalue <- raw_pvalue
  }

  # ============================================================================
  # 7. 创建可视化
  # ============================================================================

  # 创建基础tidyplot
  p <- data %>%
    tidyplot(x = Group, y = Abundance, color = Group)

  # 根据plot_type添加主图层
  if (plot_type == "boxplot") {
    if (filled) {
      p <- p %>% add_boxplot(alpha = main_alpha)
    } else {
      p <- p %>% add_boxplot(alpha = 0, linewidth = border_width)
    }
  } else if (plot_type == "violin") {
    if (filled) {
      p <- p %>% add_violin(alpha = main_alpha)
    } else {
      p <- p %>% add_violin(alpha = 0, linewidth = border_width)
    }
  } else if (plot_type == "barplot") {
    if (filled) {
      p <- p %>%
        add_mean_bar(alpha = main_alpha) %>%
        add_sem_errorbar()
    } else {
      p <- p %>%
        add_mean_bar(alpha = 0.1) %>%
        add_sem_errorbar()
    }
  } else if (plot_type == "dotplot") {
    p <- p %>%
      add_mean_dot(size = 4) %>%
      add_sem_errorbar()
  }

  # ============================================================================
  # 8. 添加散点图层
  # ============================================================================
  if (add_points) {
    if (point_type == "jitter") {
      p <- p %>%
        add_data_points_jitter(
          size = point_size,
          alpha = point_alpha,
          shape = point_shape,
          color = point_color,
          stroke = point_stroke
        )
    } else if (point_type == "beeswarm") {
      p <- p %>%
        add_data_points_beeswarm(
          size = point_size,
          alpha = point_alpha,
          shape = point_shape,
          color = point_color,
          stroke = point_stroke
        )
    } else if (point_type == "center") {
      p <- p %>%
        add_data_points(
          size = point_size,
          alpha = point_alpha,
          shape = point_shape,
          color = point_color,
          stroke = point_stroke
        )
    }
  }

  # ============================================================================
  # 9. 添加显著性标记
  # ============================================================================
  if (add_sig_label != "none") {
    # 确定检验方法的字符串
    test_method_str <- ifelse(test_method == "t.test", "t_test", "wilcox_test")

    if (add_sig_label == "asterisk") {
      # 添加星号标记
      p <- p %>%
        add_test_asterisks(
          method = test_method_str,
          p.adjust.method = p_adjust_method,
          label.size = sig_label_size / ggplot2::.pt,
          hide.ns = sig_hide_ns,
          hide_info = TRUE
        )
    } else if (add_sig_label == "pvalue") {
      # 添加p值标记
      p <- p %>%
        add_test_pvalue(
          method = test_method_str,
          p.adjust.method = p_adjust_method,
          label.size = sig_label_size / ggplot2::.pt,
          hide.ns = sig_hide_ns,
          hide_info = TRUE
        )
    }
  }

  # ============================================================================
  # 10. 格式化caption中的p值和效应量文本
  # ============================================================================
  if (show_pvalue) {
    if (p_adjust_method != "none") {
      p_prefix <- paste0("调整后p (", p_adjust_method, ")")
    } else {
      p_prefix <- "p"
    }
    p_text <- format_pvalue(display_pvalue, p_prefix)

    # 添加效应量信息
    if (!is.null(effect_size_result)) {
      effect_text <- paste0(
        effect_size_result$type,
        " = ",
        round(effect_size_result$value, 3),
        " (",
        effect_size_result$interpretation,
        ")"
      )
      subtitle_text <- paste(
        "物种:",
        species_name,
        "|",
        p_text,
        "|",
        effect_text
      )
    } else {
      subtitle_text <- paste("物种:", species_name, "|", p_text)
    }
  } else {
    subtitle_text <- paste("物种:", species_name)
  }

  # ============================================================================
  # 11. 调整图形属性
  # ============================================================================
  p <- p %>%
    adjust_title(title = plot_title) %>%
    adjust_caption(caption = subtitle_text) %>%
    adjust_x_axis_title(title = xlab) %>%
    adjust_y_axis_title(title = ylab)

  # 应用自定义配色
  if (!is.null(colors)) {
    p <- p %>% adjust_colors(new_colors = colors)
  }

  # 移除图例
  p <- p %>% remove_legend()

  # ============================================================================
  # 12. 返回结果
  # ============================================================================
  result_list <- list(
    plot = p,
    test_result = test_result,
    raw_pvalue = raw_pvalue,
    summary_stats = summary_stats
  )

  # 添加调整后的p值
  if (!is.null(adjusted_pvalue)) {
    result_list$adjusted_pvalue <- adjusted_pvalue
    result_list$adjustment_method <- p_adjust_method
  }

  # 添加效应量
  if (!is.null(effect_size_result)) {
    result_list$effect_size <- effect_size_result
  }

  return(result_list)
}

# ============================================================================
# 辅助函数（保持不变）
# ============================================================================

#' 计算两组比较的效应量
#' @keywords internal
calculate_effect_size <- function(group1, group2, method, paired) {
  n1 <- length(group1)
  n2 <- length(group2)

  if (method == "t.test") {
    if (paired) {
      diff <- group1 - group2
      d <- mean(diff) / sd(diff)
      effect_type <- "Cohen's d (paired)"
    } else {
      pooled_sd <- sqrt(
        ((n1 - 1) * var(group1) + (n2 - 1) * var(group2)) / (n1 + n2 - 2)
      )
      d <- (mean(group1) - mean(group2)) / pooled_sd
      effect_type <- "Cohen's d"
    }

    interpretation <- ifelse(
      abs(d) < 0.2,
      "negligible",
      ifelse(abs(d) < 0.5, "small", ifelse(abs(d) < 0.8, "medium", "large"))
    )

    return(list(
      value = d,
      type = effect_type,
      interpretation = interpretation
    ))
  } else {
    test_obj <- wilcox.test(group1, group2, paired = paired, exact = FALSE)
    N <- n1 + n2
    U <- test_obj$statistic
    r <- abs(1 - 2 * U / (n1 * n2))

    interpretation <- ifelse(
      r < 0.1,
      "negligible",
      ifelse(r < 0.3, "small", ifelse(r < 0.5, "medium", "large"))
    )

    return(list(
      value = r,
      type = "rank-biserial r",
      interpretation = interpretation
    ))
  }
}

#' 计算Kruskal-Wallis检验的效应量
#' @keywords internal
calculate_kruskal_effect_size <- function(kruskal_result, n) {
  H <- kruskal_result$statistic
  k <- kruskal_result$parameter + 1
  epsilon_sq <- H / ((n^2 - 1) / (n + 1))

  interpretation <- ifelse(
    epsilon_sq < 0.01,
    "negligible",
    ifelse(
      epsilon_sq < 0.06,
      "small",
      ifelse(epsilon_sq < 0.14, "medium", "large")
    )
  )

  return(list(
    value = epsilon_sq,
    type = "epsilon-squared",
    interpretation = interpretation
  ))
}

#' 计算ANOVA的效应量
#' @keywords internal
calculate_anova_effect_size <- function(data) {
  model <- aov(Abundance ~ Group, data = data)
  ss <- summary(model)[[1]]$`Sum Sq`
  eta_sq <- ss[1] / sum(ss)

  interpretation <- ifelse(
    eta_sq < 0.01,
    "negligible",
    ifelse(eta_sq < 0.06, "small", ifelse(eta_sq < 0.14, "medium", "large"))
  )

  return(list(
    value = eta_sq,
    type = "eta-squared",
    interpretation = interpretation
  ))
}

#' 格式化p值文本
#' @keywords internal
format_pvalue <- function(pvalue, prefix = "p") {
  if (is.na(pvalue)) {
    return(paste0(prefix, " = NA"))
  } else if (pvalue < 0.001) {
    return(paste0(prefix, " < 0.001***"))
  } else if (pvalue < 0.01) {
    return(paste0(prefix, " < 0.01**"))
  } else if (pvalue < 0.05) {
    return(paste0(prefix, " < 0.05*"))
  } else {
    return(paste0(prefix, " = ", round(pvalue, 3)))
  }
}

# ============================================================================
# 批量比较函数
# ============================================================================

#' 批量比较多个物种的丰度
#'
#' @param abundance_matrix 物种丰度矩阵
#' @param grouping_factor 分组因子
#' @param species_list 要分析的物种名称向量，默认为NULL（分析所有物种）
#' @param p_adjust_method 多重检验校正方法，默认为"BH"
#' @param return_plots 是否返回所有图形，默认为FALSE
#' @param ... 其他传递给compare_species_abundance_v2的参数
#'
#' @return 返回结果数据框，或包含results和plots的列表
#' @export
batch_compare_species <- function(
  abundance_matrix,
  grouping_factor,
  species_list = NULL,
  p_adjust_method = "BH",
  return_plots = FALSE,
  ...
) {
  if (is.null(species_list)) {
    species_list <- colnames(abundance_matrix)
  }

  results <- list()
  plots <- list()

  for (species in species_list) {
    tryCatch(
      {
        result <- compare_species_abundance_v2(
          abundance_matrix = abundance_matrix,
          grouping_factor = grouping_factor,
          species_name = species,
          p_adjust_method = "none",
          show_pvalue = FALSE,
          ...
        )

        results[[species]] <- list(
          species = species,
          pvalue = result$raw_pvalue,
          effect_size = if (!is.null(result$effect_size)) {
            result$effect_size$value
          } else {
            NA
          },
          effect_type = if (!is.null(result$effect_size)) {
            result$effect_size$type
          } else {
            NA
          },
          effect_interpretation = if (!is.null(result$effect_size)) {
            result$effect_size$interpretation
          } else {
            NA
          }
        )

        if (return_plots) {
          plots[[species]] <- result$plot
        }
      },
      error = function(e) {
        warning(paste0("物种 ", species, " 分析失败: ", e$message))
      }
    )
  }

  results_df <- do.call(
    rbind,
    lapply(results, function(x) {
      data.frame(
        species = x$species,
        raw_pvalue = x$pvalue,
        effect_size = x$effect_size,
        effect_type = x$effect_type,
        effect_interpretation = x$effect_interpretation,
        stringsAsFactors = FALSE
      )
    })
  )

  if (p_adjust_method != "none") {
    results_df$adjusted_pvalue <- p.adjust(
      results_df$raw_pvalue,
      method = p_adjust_method
    )
  }

  pval_col <- if (p_adjust_method != "none") "adjusted_pvalue" else "raw_pvalue"
  results_df$significance <- ifelse(
    results_df[[pval_col]] < 0.001,
    "***",
    ifelse(
      results_df[[pval_col]] < 0.01,
      "**",
      ifelse(results_df[[pval_col]] < 0.05, "*", "ns")
    )
  )

  results_df <- results_df[order(results_df[[pval_col]]), ]

  if (return_plots) {
    return(list(results = results_df, plots = plots))
  } else {
    return(results_df)
  }
}



library(ClusterGVis)

# load data
data(exps)

# check
head(exps,3)
#           zygote  t2.cell  t4.cell  t8.cell   tmorula blastocyst
# Oog4   1.3132282 1.237078 1.325978 1.262073 0.6549312  0.2067114
# Psmd9  1.0917337 1.315989 1.174417 1.064756 0.8685598  0.4845448
# Sephs2 0.9859232 1.201026 1.123076 1.084673 0.8878931  0.7174088

# check optimal cluster numbers
getClusters(obj = exps)







# 设置工作目录(改成您的实际路径)
setwd("Scripts\\r_functions\\dev\\enrich_plot")

# 加载函数
source("function.R")


library(clusterProfiler)
library(org.Hs.eg.db)

# 准备测试基因
test_genes <- c("NFKB1", "TNF", "IL1B", "IL6", "CXCL8", 
                "PTGS2", "MMP9", "TLR4", "MAPK14")

# 转换ID
gene_ids <- bitr(test_genes, fromType = "SYMBOL", 
                 toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# KEGG富集分析
kegg_test <- enrichKEGG(gene = gene_ids$ENTREZID, 
                        organism = "hsa",
                        pvalueCutoff = 1)

# 检查结果
print(paste("富集到的通路数:", nrow(kegg_test)))


# 准备模拟的差异表达数据
set.seed(123)
diff_test <- data.frame(
  gene = test_genes,
  log2FC = runif(length(test_genes), -2, 3),
  padj = runif(length(test_genes), 0.001, 0.05)
)

# 测试函数(返回ggplot对象)
p <- plot_enrich_trinity(
  enrich_result = kegg_test,
  top_n = 4,  # 只展示4个通路,测试更快
  diff_data = diff_test
)

# 查看对象类型
class(p)



plot_enrich_trinity(
  enrich_result = kegg_test,
  top_n = 4,
  diff_data = diff_test,
  output_file = "test_output.png",
  width = 12,
  height = 8
)

source('example.R')




# 1. 安装推荐包(获得最佳效果)
install.packages(c("cols4all", "colorspace"))

BiocManager::install('ggsankey')
pak::pkg_install('ggsankey')
devtools::install_github("davidsjoberg/ggsankey")
# 2. 重新加载函数
source("function.R")

# 3. 测试
plots <- plot_enrich_trinity(
  enrich_result = kegg_test,
  top_n = 6,
  diff_data = diff_test
)

# 4. 查看桑基图
plots$sankey_plot

# 5. 保存测试
ggsave("test_v1.3_improved_colors.png", plots$sankey_plot, 
       width = 10, height = 8)