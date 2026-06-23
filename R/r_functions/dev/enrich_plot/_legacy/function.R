# ==============================================================================
# plot_enrich_trinity
# 创建: 2025-04-13
# 功能: 富集分析三位一体可视化(桑基图+气泡图+条形图组合)
# ==============================================================================

#' 富集分析三位一体可视化
#'
#' 将富集分析结果以三位一体的形式展示:左侧条形图显示基因差异表达,
#' 中间桑基图展示基因-通路关系,右侧气泡图显示通路富集程度。
#'
#' @param enrich_result 富集分析结果,可以是:
#'   - enrichResult对象(来自clusterProfiler::enrichKEGG/enrichGO)
#'   - data.frame (必须包含特定列,见Details)
#' @param pathways 要展示的通路名称向量。如果为NULL,则自动选择Top N个通路
#' @param top_n 当pathways=NULL时,选择前几个通路,默认6个
#' @param diff_data 差异表达数据框,包含列: gene, log2FC, padj。如果为NULL,
#'   则不显示左侧条形图,或使用随机数据作为示例
#' @param gene_id_type 基因ID类型,"ENTREZID"或"SYMBOL",默认自动检测
#' @param convert_id 是否进行ID转换(ENTREZID <-> SYMBOL),默认TRUE
#' @param organism 物种注释数据库,如"org.Hs.eg.db"(人类),默认自动检测
#' @param color_pathway 通路颜色方案,默认使用dittoSeq::dittoColors()
#' @param color_pvalue 气泡图p值颜色渐变,向量c(low, high)
#' @param color_fc 条形图log2FC颜色渐变,向量c(low, high)
#' @param bubble_size_range 气泡大小范围,向量c(min, max)
#' @param width 输出图片宽度(英寸),默认12
#' @param height 输出图片高度(英寸),默认8
#' @param output_file 输出文件路径,如"result.png"。如果为NULL,则返回ggplot对象
#' @param combine 是否组合三个图形,默认FALSE返回独立对象列表
#' @param layout_widths 三个面板的宽度比例,向量c(left, middle, right),仅在combine=TRUE时使用
#' @param show_gene_label 是否在条形图上显示log2FC标签,默认TRUE
#' @param ... 其他参数,传递给主题设置
#'
#' @details
#'
#' **输入数据格式要求**:
#'
#' 1. enrich_result (富集分析结果) 必须包含以下列:
#'    - Description: 通路/GO term名称
#'    - geneID: 富集的基因ID,用"/"分隔的字符串
#'    - pvalue: 原始p值
#'    - p.adjust: 校正后的p值
#'    - Count: 富集的基因数量
#'    - BgRatio: 背景基因比例,格式"M/N"
#'
#' 2. diff_data (可选,差异表达数据) 需包含:
#'    - gene: 基因ID(需与enrich_result中的ID类型匹配)
#'    - log2FC: log2倍数变化
#'    - padj: 校正后的p值
#'
#' **ID类型说明**:
#' - clusterProfiler的enrichKEGG默认使用ENTREZID
#' - enrichGO可能使用SYMBOL或ENTREZID
#' - 差异表达数据通常使用SYMBOL
#' - 函数会自动尝试转换,也可以通过convert_id=FALSE关闭
#'
#' @return
#' - 如果combine=FALSE(默认): 返回包含三个ggplot对象的列表:
#'   - $bar_plot: 左侧基因条形图
#'   - $sankey_plot: 中间基因-通路桑基图
#'   - $bubble_plot: 右侧通路气泡图
#' - 如果combine=TRUE且output_file=NULL: 返回组合的patchwork对象
#' - 如果指定output_file: 保存图片并返回文件路径(此时自动combine=TRUE)
#'
#' @examples
#' # 示例1: 返回独立对象(推荐)
#' \dontrun{
#' library(clusterProfiler)
#' library(org.Hs.eg.db)
#'
#' # KEGG富集分析
#' kegg_result <- enrichKEGG(gene = my_genes, organism = "hsa")
#'
#' # 返回三个独立的ggplot对象
#' plots <- plot_enrich_trinity(
#'   enrich_result = kegg_result,
#'   top_n = 6,
#'   diff_data = my_diff_genes
#' )
#'
#' # 查看桑基图
#' print(plots$sankey_plot)
#'
#' # 查看气泡图
#' print(plots$bubble_plot)
#'
#' # 自定义组合(桑基图+气泡图,最常用)
#' library(patchwork)
#' plots$sankey_plot + plots$bubble_plot +
#'   plot_layout(widths = c(3, 1))
#' }
#'
#' # 示例2: 自动组合并保存
#' \dontrun{
#' # 指定output_file会自动组合并保存
#' plot_enrich_trinity(
#'   enrich_result = kegg_result,
#'   top_n = 6,
#'   diff_data = my_diff_genes,
#'   output_file = "trinity_plot.png"  # 自动combine=TRUE
#' )
#' }
#'
#' # 示例3: 手动设置combine=TRUE
#' \dontrun{
#' combined <- plot_enrich_trinity(
#'   enrich_result = kegg_result,
#'   top_n = 6,
#'   diff_data = my_diff_genes,
#'   combine = TRUE  # 返回组合后的对象
#' )
#' print(combined)
#' }
#'
#' @export
#' @importFrom dplyr filter select mutate as_tibble left_join
#' @importFrom stringr str_split
#' @importFrom ggplot2 ggplot aes geom_point geom_rect geom_text theme_void theme_classic scale_color_gradient scale_fill_gradient coord_flip labs theme element_text element_blank unit ggsave
#' @importFrom ggalluvial geom_alluvium geom_stratum
#' @importFrom patchwork plot_layout plot_annotation
#' @note
#' 推荐安装以下包以获得最佳效果:
#' - cols4all: 丰富的调色板
#' - colorspace: 颜色调整功能
#' - dittoSeq: 备选调色板
plot_enrich_trinity <- function(
  enrich_result,
  pathways = NULL,
  top_n = 6,
  diff_data = NULL,
  gene_id_type = c("auto", "ENTREZID", "SYMBOL"),
  convert_id = TRUE,
  organism = NULL,
  color_pathway = NULL,
  color_pvalue = c('#4393C3', '#D6604D'),
  color_fc = c("skyblue", "red"),
  bubble_size_range = c(5, 10),
  width = 12,
  height = 8,
  output_file = NULL,
  combine = FALSE,
  layout_widths = c(1, 0.8, 0.6),
  show_gene_label = TRUE,
  ...
) {
  # 1. 参数验证和数据准备 ----

  # 加载必需包
  required_packages <- c(
    "dplyr",
    "stringr",
    "ggplot2",
    "ggalluvial",
    "patchwork",
    "ggnewscale",
    "tidyr"
  )
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste0("需要安装包: ", pkg, "\n运行: install.packages('", pkg, "')"))
    }
  }

  # 转换enrichResult对象为data.frame
  if (inherits(enrich_result, "enrichResult")) {
    enrich_df <- as.data.frame(enrich_result)
  } else if (is.data.frame(enrich_result)) {
    enrich_df <- enrich_result
  } else {
    stop("enrich_result必须是enrichResult对象或data.frame")
  }

  # 检查必需列
  required_cols <- c(
    "Description",
    "geneID",
    "pvalue",
    "p.adjust",
    "Count",
    "BgRatio"
  )
  missing_cols <- setdiff(required_cols, colnames(enrich_df))
  if (length(missing_cols) > 0) {
    stop(paste0("缺少必需列: ", paste(missing_cols, collapse = ", ")))
  }

  # 计算RichFactor
  if (!"RichFactor" %in% colnames(enrich_df)) {
    enrich_df$RichFactor <-
      enrich_df$Count / as.numeric(sub("/.*", "", enrich_df$BgRatio))
  }

  # 2. 通路选择 ----

  # 筛选显著通路(可选,根据p.adjust)
  # 这里不做筛选,由用户决定

  if (is.null(pathways)) {
    # 自动选择Top N (按RichFactor或p.adjust排序)
    enrich_df <- enrich_df[order(enrich_df$p.adjust), ]
    pathways <- head(enrich_df$Description, top_n)
    message(paste0("自动选择Top ", top_n, "个通路"))
  }

  # 提取选中的通路
  enrich_subset <- enrich_df[
    match(pathways, enrich_df$Description),
    ,
    drop = FALSE
  ]
  enrich_subset <- na.omit(enrich_subset) # 移除未找到的通路

  if (nrow(enrich_subset) == 0) {
    stop("未找到指定的通路,请检查通路名称")
  }

  # 按RichFactor排序(用于图形展示顺序)
  enrich_subset <- enrich_subset[
    order(enrich_subset$RichFactor, decreasing = FALSE),
    ,
    drop = FALSE
  ]
  enrich_subset$Pathway <- factor(
    enrich_subset$Description,
    levels = enrich_subset$Description
  )

  # 3. 构建桑基图数据 ----

  library(dplyr)
  library(stringr)
  library(tidyr)

  sankey_gene <- enrich_subset %>%
    as_tibble() %>%
    mutate(gene = str_split(geneID, "/")) %>%
    unnest(cols = gene) %>%
    mutate(pathway = as.character(Pathway)) %>%
    dplyr::select(gene, pathway, pvalue, Count, RichFactor)

  # 4. ID转换(如果需要) ----

  gene_id_type <- match.arg(gene_id_type)

  if (convert_id && gene_id_type != "SYMBOL") {
    # 尝试转换为SYMBOL
    if (is.null(organism)) {
      # 尝试自动检测物种
      organism <- "org.Hs.eg.db" # 默认人类
      message("未指定organism,默认使用人类注释(org.Hs.eg.db)")
    }

    if (!requireNamespace(organism, quietly = TRUE)) {
      warning(paste0("包", organism, "未安装,无法进行ID转换"))
      warning("运行: BiocManager::install('", organism, "')")
    } else {
      # 执行ID转换
      org_db <- get(organism, asNamespace(organism))

      tryCatch(
        {
          library(clusterProfiler)
          gene_convert <- bitr(
            sankey_gene$gene,
            fromType = "ENTREZID",
            toType = "SYMBOL",
            OrgDb = org_db
          )

          sankey_gene$gene <- gene_convert$SYMBOL[match(
            sankey_gene$gene,
            gene_convert$ENTREZID
          )]
          sankey_gene <- na.omit(sankey_gene) # 移除无法转换的
        },
        error = function(e) {
          warning("ID转换失败,使用原始ID: ", e$message)
        }
      )
    }
  }

  # 5. 绘制桑基图(使用ggalluvial,改进配色) ----

  library(ggalluvial)
  library(ggplot2)

  # 改进的配色方案(参考原文)
  n_genes <- length(unique(sankey_gene$gene))
  n_pathways <- length(unique(sankey_gene$pathway))
  n_total_nodes <- n_genes + n_pathways

  if (is.null(color_pathway)) {
    # 尝试使用cols4all获取丰富的调色板
    if (
      requireNamespace("cols4all", quietly = TRUE) &&
        requireNamespace("colorspace", quietly = TRUE)
    ) {
      # 获取所有节点的颜色
      all_colors <- cols4all::c4a(
        "kovesi.rainbow_bgyr_35_85_c73",
        n = n_total_nodes
      )

      # 调整颜色:降低饱和度并调亮(参考原文)
      all_colors <- colorspace::lighten(
        colorspace::desaturate(all_colors, amount = 0.8),
        amount = 0.6
      )

      # 分配给基因和通路
      gene_colors <- all_colors[1:n_genes]
      pathway_colors <- all_colors[(n_genes + 1):n_total_nodes]

      # 创建颜色映射
      names(gene_colors) <- unique(sankey_gene$gene)
      names(pathway_colors) <- unique(sankey_gene$pathway)
    } else {
      # 备选方案:使用dittoSeq或rainbow
      if (requireNamespace("dittoSeq", quietly = TRUE)) {
        all_colors <- dittoSeq::dittoColors()[1:n_total_nodes]
      } else {
        all_colors <- rainbow(n_total_nodes)
      }

      # 分配颜色
      gene_colors <- all_colors[1:n_genes]
      pathway_colors <- all_colors[(n_genes + 1):n_total_nodes]
      names(gene_colors) <- unique(sankey_gene$gene)
      names(pathway_colors) <- unique(sankey_gene$pathway)
    }
  } else {
    # 用户提供了颜色,只用于通路
    pathway_colors <- color_pathway
    names(pathway_colors) <- unique(sankey_gene$pathway)

    # 基因使用彩虹色(低饱和度)
    gene_colors <- rainbow(n_genes, s = 0.3, v = 0.9)
    names(gene_colors) <- unique(sankey_gene$gene)
  }

  # 合并所有节点的颜色
  all_node_colors <- c(gene_colors, pathway_colors)

  p_sankey <- ggplot(
    data = sankey_gene,
    aes(
      axis1 = gene,
      axis2 = pathway,
      y = 1
    )
  ) +
    geom_alluvium(
      aes(fill = gene), # 流按基因着色,颜色更丰富
      curve_type = "quintic",
      show.legend = FALSE,
      alpha = 0.7 # 增加透明度
    ) +
    geom_stratum(
      aes(fill = after_stat(stratum)), # 节点颜色按节点名称
      show.legend = FALSE
    ) +
    geom_text(
      stat = "stratum",
      aes(label = after_stat(stratum)),
      size = 4,
      color = "black"
    ) +
    scale_x_discrete(limits = c("Gene", "Pathway")) +
    scale_fill_manual(values = all_node_colors) + # 统一的颜色映射
    theme_void() +
    theme(
      axis.title = element_text(size = 18, color = "black"),
      axis.text = element_text(size = 16, color = "black")
    )

  # 6. 提取节点位置并绘制气泡图 ----

  pbuilt <- ggplot2::ggplot_build(p_sankey)
  node_data <- pbuilt$data[[2]]

  # 提取右侧通路节点位置
  right_nodes <- node_data %>%
    filter(x == 2) %>%
    mutate(y = (ymax + ymin) / 2) %>%
    dplyr::select(stratum, y, ymin, ymax)

  enrich_subset$position_y <- right_nodes$y[match(
    as.character(enrich_subset$Pathway),
    right_nodes$stratum
  )]

  library(ggnewscale)

  p_bubble <- ggplot(enrich_subset, aes(x = Count, y = position_y)) +
    geom_point(aes(color = pvalue, size = Count)) +
    scale_size(range = bubble_size_range) +
    theme_void() +
    scale_color_gradient(
      low = color_pvalue[1],
      high = color_pvalue[2],
      name = "p-value"
    ) +
    labs(x = "", y = "") +
    theme(
      legend.key.size = unit(0.6, "cm"),
      legend.position = "right"
    )

  # 7. 绘制基因条形图(如果提供差异表达数据) ----

  if (!is.null(diff_data)) {
    # 验证差异表达数据格式
    required_diff_cols <- c("gene", "log2FC", "padj")
    if (!all(required_diff_cols %in% colnames(diff_data))) {
      warning("diff_data缺少必需列,将跳过条形图")
      p_bar <- ggplot() + theme_void() # 空白图
    } else {
      # 提取左侧基因节点位置
      left_nodes <- node_data %>%
        filter(x == 1) %>%
        mutate(y = (ymax + ymin) / 2) %>%
        dplyr::select(stratum, y, ymin, ymax)

      # 匹配差异表达数据
      left_nodes <- left_nodes %>%
        left_join(diff_data, by = c("stratum" = "gene"))

      # 处理缺失值
      left_nodes$log2FC[is.na(left_nodes$log2FC)] <- 0
      left_nodes$padj[is.na(left_nodes$padj)] <- 1

      p_bar <- ggplot(left_nodes, aes(fill = -log10(padj))) +
        geom_rect(aes(
          xmin = ymin,
          xmax = ymax,
          ymin = 0,
          ymax = log2FC
        )) +
        coord_flip() +
        scale_fill_gradient(low = color_fc[1], high = color_fc[2]) +
        scale_y_reverse() +
        theme_classic() +
        theme(
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          legend.position = "none"
        )

      if (show_gene_label) {
        p_bar <- p_bar +
          geom_text(
            aes(
              x = y,
              y = 0,
              label = paste0("FC=", round(log2FC, 2))
            ),
            hjust = 1.2,
            size = 3
          )
      }
    }
  } else {
    # 如果没有差异表达数据,创建空白图或示例数据
    message("未提供diff_data,将使用随机数据作为示例")

    # 提取左侧基因节点
    left_nodes <- node_data %>%
      filter(x == 1) %>%
      mutate(y = (ymax + ymin) / 2) %>%
      dplyr::select(stratum, y, ymin, ymax)

    # 随机生成log2FC
    set.seed(123)
    left_nodes$log2FC <- runif(nrow(left_nodes), -2, 3)
    left_nodes$padj <- runif(nrow(left_nodes), 0.001, 0.05)

    p_bar <- ggplot(left_nodes, aes(fill = -log10(padj))) +
      geom_rect(aes(
        xmin = ymin,
        xmax = ymax,
        ymin = 0,
        ymax = log2FC
      )) +
      coord_flip() +
      scale_fill_gradient(low = color_fc[1], high = color_fc[2]) +
      scale_y_reverse() +
      theme_classic() +
      theme(
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        legend.position = "none"
      )
  }

  # 8. 返回结果 ----

  # 如果指定了output_file,强制combine=TRUE并保存
  if (!is.null(output_file)) {
    combine <- TRUE
  }

  if (combine) {
    # 组合三个图形
    library(patchwork)

    # 定义布局 (左-中-右)
    layout <- "
      AACCBB
      AACCBB
      AACCBB
    "

    combined_plot <- p_bar +
      p_sankey +
      p_bubble +
      plot_layout(
        design = layout,
        widths = layout_widths
      )

    # 如果指定了输出文件,保存
    if (!is.null(output_file)) {
      ggsave(
        output_file,
        plot = combined_plot,
        width = width,
        height = height,
        dpi = 300
      )
      message(paste0("组合图已保存至: ", output_file))
      return(invisible(output_file))
    } else {
      # 返回组合的ggplot对象
      return(combined_plot)
    }
  } else {
    # 返回三个独立对象的列表
    result <- list(
      bar_plot = p_bar,
      sankey_plot = p_sankey,
      bubble_plot = p_bubble
    )

    message("返回三个独立的ggplot对象: $bar_plot, $sankey_plot, $bubble_plot")
    message("提示: 使用 plots$sankey_plot 查看桑基图")
    message("提示: 使用 plots$sankey_plot + plots$bubble_plot 组合展示")

    return(result)
  }
}



# 使用火山图进行标注 --------------------------------------------------------------



