# 1. 加载必要的库
library(ggplot2)
library(dplyr)

# 2. 读取数据 (使用 fill=TRUE 兼容瑕疵行)
df <- read.table("aixiao-WXP-CLR.txt", 
                 header = TRUE, 
                 stringsAsFactors = FALSE, 
                 fill = TRUE,      
                 quote = "",       
                 comment.char = "") 

# 3. 准备绘图数据
plot_data <- df %>%
  mutate(chrom = factor(chrom, levels = c(1:33, "Z"))) %>%
  group_by(chrom) %>%
  summarise(chr_len = max(stop)) %>%
  mutate(tot = cumsum(as.numeric(chr_len)) - chr_len) %>%
  select(-chr_len) %>%
  left_join(df, ., by = c("chrom" = "chrom")) %>%
  arrange(chrom, start) %>%
  mutate(bp_cum = start + tot)

# --- X轴坐标显示保持不变 ---
keep_axis <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", 
               "12", "14", "17", "21", "33", "Z")

axis_df <- plot_data %>%
  group_by(chrom) %>%
  summarize(center = (max(bp_cum) + min(bp_cum)) / 2) %>%
  mutate(display_label = ifelse(chrom %in% keep_axis, as.character(chrom), ""))

# --- 4. 修复颜色的关键代码：锁定颜色映射 ---
color_palette <- c("#e54d38", "#51b9d2", "#db8933", "#f39c81", 
                   "#828fb2", "#97cfc4", "#3b5586", "#ae9c84")

# 建立一个包含所有可能染色体的列表
all_levels <- c(as.character(1:33), "Z")

# 为每一个染色体编号“点名”分配固定颜色
final_colors <- rep(color_palette, length.out = length(all_levels))
names(final_colors) <- all_levels

# --- 重点修改：手动指定 Z 染色体为红色 ---
# 将 Z 的颜色强制设置为 color_palette 的第一个颜色（红色 #e54d38）
final_colors["Z"] <- color_palette[1]

# 5. 开始绘图
# 注意：color 映射改为 as.character(chrom) 以便匹配字典名字
p <- ggplot(plot_data, aes(x = bp_cum, y = xpclr_norm, color = as.character(chrom))) +
  geom_point(alpha = 0.8, size = 1) +
  # 使用锁定好的颜色映射表
  scale_color_manual(values = final_colors) +
  # 保持你原来的 X 轴设置
  scale_x_continuous(label = axis_df$display_label, breaks = axis_df$center) +
  geom_hline(yintercept = quantile(plot_data$xpclr_norm, 0.99, na.rm=TRUE), 
             linetype = "dashed", color = "grey40") +
  labs(title = "XP-CLR Manhattan Plot",
       x = "Chromosome",
       y = "Normalized XP-CLR Score") +
  theme_minimal() +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14), 
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 0, size = 9, vjust = 0.5)
  )

# 6. 显示并保存
print(p)

ggsave("XPCLR_Manhattan_Plot_Clean.png", 
       plot = p, width = 12, height = 6, dpi = 300, bg = "white")

# 7. 提取 Top 5% 数据
threshold_5 <- quantile(df$xpclr_norm, 0.95, na.rm = TRUE)
df_top5 <- df %>%
  filter(xpclr_norm >= threshold_5) %>%
  select(2, 3, 4) %>%
  arrange(chrom, start)

write.table(df_top5, "XP-CLR_top5_R_output.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)

message("修复完成：颜色已锁定映射，X 轴保持精简显示。")