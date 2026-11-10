#!/usr/bin/env Rscript

# 1. 读取数据 - 现在有8列
xpehh_data <- read.table("aixiao-Wtu.xpehh", header = FALSE)
colnames(xpehh_data) <- c("pos", "gpos", "p1", "ihh1", "p2", "ihh2", "xpehh_raw", "xpehh_norm")

# 根据第8列（xpehh_norm）从大到小排序
sorted_data <- xpehh_data[order(-xpehh_data$xpehh_norm), ]

# 计算5%的阈值位置
total_count <- nrow(sorted_data)
top_5_percent_count <- ceiling(total_count * 0.05)
bottom_5_percent_count <- ceiling(total_count * 0.05)

# 获取最大的前5%
top_5_percent <- head(sorted_data, top_5_percent_count)

# 获取最小的后5%
bottom_5_percent <- tail(sorted_data, bottom_5_percent_count)

# 获取top_5_percent_count第八列的最后一位数
threshold_upper <- tail(top_5_percent$xpehh_norm, 1)

# 获取bottom_5_percent_count第八列的第一位数
threshold_lower <- head(bottom_5_percent$xpehh_norm, 1)

# 2. 提取染色体信息
extract_chr_info <- function(pos_string) {
  parts <- strsplit(pos_string, ":")[[1]]
  if (length(parts) == 2) {
    return(list(chr = parts[1], pos = as.numeric(parts[2])))
  }
  return(list(chr = NA, pos = NA))
}

chr_info <- lapply(xpehh_data$pos, extract_chr_info)
xpehh_data$chr <- sapply(chr_info, function(x) x$chr)
xpehh_data$physical_pos <- sapply(chr_info, function(x) x$pos)

# 3. 创建染色体排序（包含Z染色体）
all_chromosomes <- sort(unique(xpehh_data$chr))
numeric_chr <- all_chromosomes[grepl("^[0-9]+$", all_chromosomes)]
character_chr <- all_chromosomes[!grepl("^[0-9]+$", all_chromosomes)]
ordered_chr <- c(sort(as.numeric(numeric_chr)), sort(character_chr))

chr_to_num <- setNames(seq_along(ordered_chr), ordered_chr)
xpehh_data$chr_num <- chr_to_num[xpehh_data$chr]

# 4. 准备曼哈顿图数据
chr_max <- aggregate(physical_pos ~ chr, data = xpehh_data, max, na.rm = TRUE)
chr_max <- chr_max[order(chr_to_num[chr_max$chr]), ]
chr_max$bp_add <- c(0, cumsum(chr_max$physical_pos)[-nrow(chr_max)])

plot_data <- merge(xpehh_data, chr_max[, c("chr", "bp_add")], by = "chr")
plot_data$pos_abs <- plot_data$physical_pos + plot_data$bp_add

# 5. 计算坐标轴中心（包含所有染色体）
axis_centers <- aggregate(pos_abs ~ chr, data = plot_data, mean)
axis_centers <- axis_centers[order(chr_to_num[axis_centers$chr]), ]

# 6. 设置要显示标签的染色体
show_chr_labels <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "11", "13", "16", "19", "22", "27", "33", "Z")
axis_centers_display <- axis_centers[axis_centers$chr %in% show_chr_labels, ]

# 创建所有染色体的breaks，但只显示指定标签
all_breaks <- axis_centers$pos_abs
all_labels <- ifelse(axis_centers$chr %in% show_chr_labels, axis_centers$chr, "")

# 7. 设置8种颜色循环
color_palette <- c("#e54d38", "#51b9d2", "#db8933", "#f39c81", 
                   "#828fb2", "#97cfc4", "#3b5586", "#ae9c84")
plot_data$chr_index <- as.numeric(factor(plot_data$chr, levels = ordered_chr))
plot_data$color_group <- (plot_data$chr_index - 1) %% 8 + 1
plot_data$point_color <- color_palette[plot_data$color_group]

# 8. 绘制曼哈顿图
library(ggplot2)

signif_threshold <- 4

manhattan_plot <- ggplot(plot_data, 
                         aes(x = pos_abs, y = xpehh_norm)) +
  
  # 数据点 - 使用8种颜色循环
  geom_point(aes(color = as.factor(color_group)), alpha = 0.7, size = 0.5) +
  
  # 使用指定的8种颜色
  scale_color_manual(values = color_palette) +
  
  # 绘制动态阈值线
  geom_hline(yintercept = c(threshold_upper, threshold_lower), 
             linetype = "dashed", 
             color = "black", 
             linewidth = 0.8,
             alpha = 0.8) +
  
  # 坐标轴设置：所有染色体都有刻度线，但只显示指定标签
  scale_x_continuous(
    breaks = all_breaks,
    labels = all_labels,
    expand = expansion(mult = 0.02)
  ) +
  
  scale_y_continuous(
    expand = expansion(mult = 0.05)
  ) +
  
  labs(
    x = "Chr",
    y = "XP-EHH"
  ) +
  
  # 主题设置
  theme_classic() +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 12, face = "bold", margin = margin(t = 10)),
    axis.title.y = element_text(size = 12, face = "bold", margin = margin(r = 10)),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )


# 10. 保存结果
print(manhattan_plot)
ggsave("xpehh_manhattan_final.png", width = 14, height = 7, dpi = 300, bg = "white")
ggsave("xpehh_manhattan_final.pdf", width = 14, height = 7, bg = "white")

# 11. 输出信息
cat("最终版曼哈顿图已保存！\n")
cat("修改内容:\n")
cat("- 保留所有染色体（包括Z染色体）的刻度线\n")
cat("- 只在指定染色体下方显示标签: 1,2,3,4,5,6,7,8,9,11,13,16,19,22,27,33,Z\n")
cat("- 去掉上下阈值线，只保留中间的0值黑色虚线\n")
cat("- 阈值线置于最上层显示\n")
cat("- Y轴范围设置为-2.5到2.5\n")



# 根据第8列（xpehh_norm）从大到小排序
sorted_data <- xpehh_data[order(-xpehh_data$xpehh_norm), ]

# 计算5%的阈值位置
total_count <- nrow(sorted_data)
top_5_percent_count <- ceiling(total_count * 0.05)
bottom_5_percent_count <- ceiling(total_count * 0.05)

# 获取最大的前5%
top_5_percent <- head(sorted_data, top_5_percent_count)

# 获取最小的后5%
bottom_5_percent <- tail(sorted_data, bottom_5_percent_count)
