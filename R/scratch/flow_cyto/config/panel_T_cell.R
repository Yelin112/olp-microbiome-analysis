# panel_T_cell.R
# T 细胞面板配置：通道映射、文件分类、门控层级
# 根据 01_read_explore.R 的输出建立，FL10 为活死染料，ECD/PB450 等待确认

PANEL <- list(

  # ── 荧光通道 → 标记物映射 ──────────────────────────────────────────────────
  # channel_area: 分析用通道（-A，面积）
  # marker:       抗体/染料靶点
  # fluorochrome: 荧光素名称
  # confirmed:    是否已确认标记（FALSE = 待核实）
  channels = list(
    FL1  = list(channel_area = "FL1-A",  marker = "CD3",   fluorochrome = "FITC",      confirmed = TRUE),
    FL2  = list(channel_area = "FL2-A",  marker = "IFN-γ", fluorochrome = "PE",        confirmed = TRUE),
    FL3  = list(channel_area = "FL3-A",  marker = NA,      fluorochrome = "ECD",       confirmed = FALSE, empty = TRUE), # 仪器自带探测器，本实验未用抗体
    FL4  = list(channel_area = "FL4-A",  marker = "CD4",   fluorochrome = "PC5.5",     confirmed = TRUE,  empty = FALSE),
    FL5  = list(channel_area = "FL5-A",  marker = "CD8",   fluorochrome = "PC7",       confirmed = TRUE,  empty = FALSE),
    FL6  = list(channel_area = "FL6-A",  marker = "IL-17A",fluorochrome = "APC",       confirmed = TRUE,  empty = FALSE),
    FL7  = list(channel_area = "FL7-A",  marker = NA,      fluorochrome = "APC-A700",  confirmed = FALSE, empty = TRUE), # 仪器自带，未使用
    FL8  = list(channel_area = "FL8-A",  marker = "CD45",  fluorochrome = "APC-A750",  confirmed = TRUE,  empty = FALSE),
    FL9  = list(channel_area = "FL9-A",  marker = NA,      fluorochrome = "PB450",     confirmed = FALSE, empty = TRUE), # 仪器自带，未使用
    FL10 = list(channel_area = "FL10-A", marker = "死活",  fluorochrome = "KO525",     confirmed = TRUE,  empty = FALSE), # 活死染料
    FL11 = list(channel_area = "FL11-A", marker = NA,      fluorochrome = "Violet610", confirmed = FALSE, empty = TRUE), # 仪器自带，未使用
    FL12 = list(channel_area = "FL12-A", marker = NA,      fluorochrome = "Violet660", confirmed = FALSE, empty = TRUE), # 仪器自带，未使用
    FL13 = list(channel_area = "FL13-A", marker = NA,      fluorochrome = "Violet780", confirmed = FALSE, empty = TRUE)  # 仪器自带，未使用
  ),

  # ── 文件分类 ───────────────────────────────────────────────────────────────
  file_roles = list(
    full_sample  = c("OXA-MN.fcs", "OXA-MN-RM.fcs", "Control.fcs"),  # 全量采集（3万+事件）
    gated_pop    = c("CD3.fcs", "CD4.fcs", "CD8.fcs", "CD45.fcs"),    # 仪器预门控群体
    cytokine_pop = c("IFN-γ.fcs", "IL-17A.fcs"),                       # 细胞因子群体
    blank        = c("blank.fcs"),                                      # 空白对照
    viability    = c("死活.fcs"),                                       # 活死对照
    duplicates   = c(                                                   # 待删除的重复文件
      "CD8_20260701_195644.fcs",
      "OXA-MN_20260701_195655.fcs",
      "blank_20260701_195619.fcs",
      "死活_20260701_195701.fcs"
    )
  ),

  # ── 标准门控层级（T 细胞免疫表型）──────────────────────────────────────────
  # 每一步 list：输入文件类型、X轴通道、Y轴通道、期望输出群体
  gate_hierarchy = list(
    step1_lymph = list(
      input    = "full_sample",
      x        = "FSC-A", y = "SSC-A",
      label    = "淋巴细胞 gate"
    ),
    step2_live = list(
      input    = "step1_lymph",
      x        = "FL10-A", y = "SSC-A",  # 死活 KO525
      label    = "活细胞 gate"
    ),
    step3_cd45 = list(
      input    = "step2_live",
      x        = "FL8-A", y = "SSC-A",   # CD45 APC-A750
      label    = "CD45+ gate"
    ),
    step4_cd3 = list(
      input    = "step3_cd45",
      x        = "FL1-A", y = "SSC-A",   # CD3 FITC
      label    = "CD3+ gate"
    ),
    step5_cd4_cd8 = list(
      input    = "step4_cd3",
      x        = "FL4-A", y = "FL5-A",   # CD4 PC5.5 vs CD8 PC7
      label    = "CD4+/CD8+ gate"
    ),
    step6_cytokine = list(
      input    = "step5_cd4_cd8",         # 在 CD4+ 群上
      x        = "FL6-A", y = "FL2-A",   # IL-17A APC vs IFN-γ PE
      label    = "IL-17A+/IFN-γ+ gate"
    )
  ),

  # ── 分析通道简称（用于图标签）──────────────────────────────────────────────
  axis_labels = c(
    "FL1-A"  = "CD3 (FITC)",
    "FL2-A"  = "IFN-γ (PE)",
    "FL4-A"  = "CD4 (PC5.5)",
    "FL5-A"  = "CD8 (PC7)",
    "FL6-A"  = "IL-17A (APC)",
    "FL8-A"  = "CD45 (APC-A750)",
    "FL10-A" = "活死染料 (KO525)",
    "FSC-A"  = "FSC-A",
    "SSC-A"  = "SSC-A"
  )
)
