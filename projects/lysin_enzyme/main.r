library(tidyplots)
library(dplyr)

# 读取清洗好的数据
data <- read.csv("裂解酶//enzyme_antibacterial_data_tidy.csv")

# 设置浓度因子的顺序
data <- data |>
  mutate(
    concentration_factor = factor(
      concentration_label,
      levels = c(
        "50 μg/ml",
        "40 μg/ml",
        "5 μg/ml",
        "500 ng/ml",
        "50 ng/ml",
        "5 ng/ml",
        "0.5 ng/ml",
        "0 (对照)"
      )
    )
  )


# ==================================================
# 图1: 生长曲线 - 按裂解酶分面
# ==================================================

data |>
  tidyplot(x = time_hr, y = OD, color = concentration_factor) |>
  add_mean_line(size = 2.5) |>
  add_sem_errorbar(width = 1.2, linewidth = 1.5) |>
  adjust_x_axis_title("Time", fontsize=10) |>
  adjust_y_axis_title("OD600", fontsize=10) |>
  adjust_legend_title("Concentration", fontsize=10) |>
  adjust_colors(colors_discrete_seaside) |>
  adjust_size(width = 120, height = 80) |>
  split_plot(by = enzyme) |>
  save_plot("裂解酶//图1_生长曲线.pdf")

# ==================================================
# 图2: 终点OD值比较 (14.5小时)
# ==================================================
data_endpoint <- data |>
  filter(time_hr == 14.5)

data_endpoint |>
  tidyplot(x = concentration_factor, y = OD, color = enzyme) |>
  add_mean_bar(alpha = 0.4) |>
  add_sem_errorbar() |>
  add_data_points_jitter(size = 1.5, alpha = 0.6) |>
  adjust_x_axis_title("浓度") |>
  adjust_y_axis_title("OD值") |>
  adjust_legend_title("裂解酶") |>
  adjust_colors(colors_discrete_friendly) |>
  adjust_size(width = 100, height = 70) |>
  theme_minimal_y()

# ==================================================
# 图3: 抑制率热图（仅V9943, V261, V1468有对照组）
# ==================================================
inhibition_data <- data |>
  filter(time_hr == 14.5, enzyme != "V118") |> # V118没有0对照
  group_by(enzyme, concentration_factor) |>
  summarise(mean_OD = mean(OD, na.rm = TRUE), .groups = "drop") |>
  group_by(enzyme) |>
  mutate(
    control_OD = mean_OD[concentration_factor == "0 (对照)"],
    inhibition_rate = (1 - mean_OD / control_OD) * 100
  ) |>
  filter(concentration_factor != "0 (对照)")

inhibition_data |>
  tidyplot(x = concentration_factor, y = enzyme, color = inhibition_rate) |>
  add_heatmap() |>
  adjust_x_axis_title("浓度") |>
  adjust_y_axis_title("裂解酶") |>
  adjust_legend_title("抑制率 (%)") |>
  adjust_colors(colors_continuous_viridis) |>
  adjust_size(width = 100, height = 50) |>
  save_plot("裂解酶//图3_抑制率热图.pdf")

# ==================================================
# 图3b: V118单独的相对效果图（相对于最低浓度0.5 ng/ml）
# ==================================================
inhibition_v118 <- data |>
  filter(time_hr == 14.5, enzyme == "V118") |>
  group_by(concentration_factor) |>
  summarise(mean_OD = mean(OD, na.rm = TRUE), .groups = "drop") |>
  mutate(
    baseline_OD = mean_OD[concentration_factor == "0.5 ng/ml"],
    relative_change = ((mean_OD - baseline_OD) / baseline_OD) * 100
  ) |>
  filter(concentration_factor != "0.5 ng/ml")

inhibition_v118 |>
  tidyplot(
    x = concentration_factor,
    y = relative_change,
    color = relative_change
  ) |>
  add_mean_bar() |>
  add_data_labels(format = "%.1f%%") |>
  adjust_x_axis_title("浓度") |>
  adjust_y_axis_title("相对变化 (相对于0.5 ng/ml, %)") |>
  adjust_title("V118浓度效应") |>
  adjust_colors(colors_continuous_viridis) |>
  adjust_size(width = 80, height = 50)
# ==================================================
# 图4: V9943详细分析（点+误差线）
# ==================================================
data |>
  filter(enzyme == "V9943") |>
  tidyplot(x = time_hr, y = OD, color = concentration_factor) |>
  add_mean_dot(size = 2) |>
  add_mean_line(size = 1) |>
  add_sem_errorbar() |>
  adjust_x_axis_title("时间 (小时)") |>
  adjust_y_axis_title("OD值") |>
  adjust_legend_title("浓度") |>
  adjust_colors(colors_discrete_rainbow) |>
  adjust_size(width = 90, height = 60)

# ==================================================
# 图5: 四种裂解酶对比（选择代表性浓度）
# ==================================================
data |>
  filter(concentration_ng_ml %in% c(50000, 500, 5, 0)) |>
  mutate(
    conc_simple = case_when(
      concentration_ng_ml == 50000 ~ "高 (50 μg/ml)",
      concentration_ng_ml == 500 ~ "中 (500 ng/ml)",
      concentration_ng_ml == 5 ~ "低 (5 ng/ml)",
      concentration_ng_ml == 0 ~ "对照 (0)"
    ),
    conc_simple = factor(
      conc_simple,
      levels = c("高 (50 μg/ml)", "中 (500 ng/ml)", "低 (5 ng/ml)", "对照 (0)")
    )
  ) |>
  tidyplot(x = time_hr, y = OD, color = enzyme) |>
  add_mean_line(size = 1.2) |>
  add_sem_errorbar() |>
  adjust_x_axis_title("时间 (小时)") |>
  adjust_y_axis_title("OD值") |>
  adjust_legend_title("裂解酶") |>
  adjust_colors(colors_discrete_metro) |>
  adjust_size(width = 120, height = 60) |>
  split_plot(by = conc_simple)

# ==================================================
# 图6: V118特殊浓度梯度分析（点+误差线）
# ==================================================
data |>
  filter(enzyme == "V118") |>
  tidyplot(x = time_hr, y = OD, color = concentration_factor) |>
  add_mean_dot(size = 2) |>
  add_mean_line(size = 1) |>
  add_sem_errorbar() |>
  adjust_x_axis_title("时间 (小时)") |>
  adjust_y_axis_title("OD值") |>
  adjust_legend_title("浓度") |>
  adjust_colors(colors_discrete_rainbow) |>
  adjust_size(width = 90, height = 60)

# ==================================================
# 图7: 浓度-响应曲线（终点，对数刻度）
# ==================================================
data_endpoint_nonzero <- data |>
  filter(time_hr == 14.5, concentration_ng_ml > 0) |>
  mutate(log_conc = log10(concentration_ng_ml))

data_endpoint_nonzero |>
  tidyplot(x = log_conc, y = OD, color = enzyme) |>
  add_data_points(size = 2.5, alpha = 0.7) |>
  add_mean_line(size = 1) |>
  adjust_x_axis_title("浓度 (log10 ng/ml)") |>
  adjust_y_axis_title("OD值") |>
  adjust_legend_title("裂解酶") |>
  adjust_colors(colors_discrete_friendly) |>
  adjust_size(width = 90, height = 60)

# ==================================================
# 图8: 最高浓度下的相对生长曲线（点+误差线）
# ==================================================
# 对于V9943, V261, V1468使用0对照；V118使用0.5 ng/ml作为基线
inhibition_over_time <- data |>
  group_by(enzyme, time_hr, concentration_factor, replicate) |>
  summarise(OD = mean(OD), .groups = "drop") |>
  group_by(enzyme, time_hr) |>
  mutate(
    baseline_OD = if_else(
      enzyme == "V118",
      mean(OD[concentration_factor == "0.5 ng/ml"]),
      mean(OD[concentration_factor == "0 (对照)"])
    ),
    relative_OD = OD / baseline_OD
  )

# 只显示最高浓度
inhibition_over_time |>
  filter(
    (enzyme %in%
      c("V9943", "V261", "V1468") &
      concentration_factor == "50 μg/ml") |
      (enzyme == "V118" & concentration_factor == "40 μg/ml")
  ) |>
  tidyplot(x = time_hr, y = relative_OD, color = enzyme) |>
  add_mean_dot(size = 2) |>
  add_mean_line(size = 1.2) |>
  add_sem_errorbar() |>
  adjust_x_axis_title("时间 (小时)") |>
  adjust_y_axis_title("相对OD") |>
  adjust_legend_title("裂解酶") |>
  adjust_colors(colors_discrete_friendly) |>
  adjust_size(width = 90, height = 60)

# ==================================================
# 图9: 各裂解酶在所有浓度下的对比（点+误差线）
# ==================================================
data |>
  tidyplot(x = time_hr, y = OD, color = enzyme) |>
  add_mean_dot(size = 1.5) |>
  add_mean_line(size = 0.8) |>
  add_sem_errorbar() |>
  adjust_x_axis_title("时间 (小时)") |>
  adjust_y_axis_title("OD值") |>
  adjust_legend_title("裂解酶") |>
  adjust_colors(colors_discrete_friendly) |>
  adjust_size(width = 120, height = 80) |>
  split_plot(by = concentration_factor) |>
  save_plot("裂解酶//all_enzymes_all_concentrations.pdf")
