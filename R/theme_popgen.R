#### ============================================================
#### theme_popgen.R —— 全仓库唯一的样式出口
####
#### 作用：把「配色 / 字号 / 点大小 / 导出尺寸 / dpi」集中在这一次定义。
####      其它作图函数**只允许**从这里取值，不许自己写死颜色和数字。
####
#### 这是"换任何 AI 作图，出来都一个样"的唯一手段：
####  AI 改的是画法，样式永远从本文件读。
####
#### 修改样式只改这一个文件，全仓库所有图同步变。
#### ============================================================

suppressPackageStartupMessages(library(ggplot2))

#### ---------- 1. 配色 ----------
#### 8 色固定顺序。前三位刻意贴近原始图的 red / blue / green，
#### 这样历史图和新图能对上眼，但全套图共用同一份色号。
PALETTE_POPGEN <- c(
  "#D62728",  # 1 红
  "#1F77B4",  # 2 蓝
  "#2CA02C",  # 3 绿
  "#FF7F0E",  # 4 橙
  "#9467BD",  # 5 紫
  "#8C564B",  # 6 棕
  "#17BECF",  # 7 青
  "#E377C2"   # 8 粉
)

#### 取色：n 个群体的前 n 个颜色；n 超过 8 时循环取用
popgen_colors <- function(n) {
  if (missing(n) || is.null(n)) return(PALETTE_POPGEN)
  rep(PALETTE_POPGEN, length.out = n)
}

#### 要 n 个可区分色（如 ADMIXTURE 的 K 个祖先成分）：
#### n <= 8 直接取色板前 n 个（相邻色对比强）；
#### n >  8 才在色板上插值，保证任何 K 都出自同一套色系
popgen_colors_ramp <- function(n) {
  if (n <= length(PALETTE_POPGEN)) return(PALETTE_POPGEN[seq_len(n)])
  grDevices::colorRampPalette(PALETTE_POPGEN)(n)
}

#### 曼哈顿图专用：染色体交替色（CMplot 的 col 参数按染色体轮换）
#### 取值为原 Fst 脚本中实际使用的 6 色，未做改动
PALETTE_CHR <- c("#fe0000", "#ffd600", "#0000fe",
                 "#7effaa", "#6294eb", "#2d8956")

#### 多染色体循环色（8 色）：XP-CLR / XP-EHH 这类有 30+ 条染色体时用，
#### 取值为原手写脚本中的 8 色，未做改动
PALETTE_CHR8 <- c("#e54d38", "#51b9d2", "#db8933", "#f39c81",
                  "#828fb2", "#97cfc4", "#3b5586", "#ae9c84")

#### 韦恩图专用：Venny 2.1 的默认四色（红、黄、绿、蓝），
#### 半透明叠加后是它那个标志性的观感
PALETTE_VENN <- c("#FF0000", "#FFCC00", "#00CC00", "#0000FF")

#### ---------- 2. 尺寸与字号 ----------
POPGEN_POINT_SIZE <- 2      # 散点大小（与原始 PCA 图的 size = 2 一致）
POPGEN_BASE_SIZE  <- 11     # theme 基准字号
POPGEN_LINE_WIDTH <- 0.5    # 坐标轴线宽
POPGEN_FONT       <- ""     # 字体族；留空 = 系统默认。Windows 中文可设 "Microsoft YaHei"

#### 导出参数：单图默认尺寸（与原始 PCA 图 ggsave 参数一致）
POPGEN_WIDTH  <- 10
POPGEN_HEIGHT <- 4.5
POPGEN_DPI    <- 500

#### ---------- 3. 主题 ----------
#### 以 theme_light() 为底 —— 原始 PCA 图用的就是它，保持不变
theme_popgen <- function(base_size = POPGEN_BASE_SIZE,
                         base_family = POPGEN_FONT) {
  theme_light(base_size = base_size, base_family = base_family) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text        = element_text(colour = "black"),
      legend.title     = element_text(face = "bold")
    )
}

#### 无边框无网格的变体：LD 衰减这类曲线图原本用的就是 theme_minimal()，
#### 保留这个观感，但字号仍走统一常量
theme_popgen_minimal <- function(base_size = POPGEN_BASE_SIZE,
                                 base_family = POPGEN_FONT) {
  theme_minimal(base_size = base_size, base_family = base_family) +
    theme(
      axis.text    = element_text(colour = "black"),
      legend.title = element_text(face = "bold"),
      panel.border = element_blank()
    )
}

#### 白底黑轴变体：XP-EHH 这类曼哈顿图原本用 theme_classic()，
#### 保留坐标轴线 + 无网格的观感，字号仍走统一常量
theme_popgen_classic <- function(base_size = POPGEN_BASE_SIZE,
                                 base_family = POPGEN_FONT) {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.text           = element_text(colour = "black"),
      legend.title        = element_text(face = "bold"),
      panel.background    = element_rect(fill = "white", colour = NA),
      plot.background     = element_rect(fill = "white", colour = NA),
      panel.grid.major    = element_blank(),
      panel.grid.minor    = element_blank()
    )
}

#### ---------- 4. 统一导出 ----------
save_popgen <- function(plot, file,
                        width  = POPGEN_WIDTH,
                        height = POPGEN_HEIGHT,
                        dpi    = POPGEN_DPI) {
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  ggsave(filename = file, plot = plot,
         width = width, height = height, dpi = dpi)
  invisible(file)
}

#### ---------- 5. 矢量/位图双份导出（可选调用） ----------
## 需要 PDF + PNG 两版时用这个，保证两版尺寸一致
save_popgen_both <- function(plot, prefix,
                            width  = POPGEN_WIDTH,
                            height = POPGEN_HEIGHT,
                            dpi    = POPGEN_DPI) {
  save_popgen(plot, paste0(prefix, ".pdf"), width, height, dpi)
  save_popgen(plot, paste0(prefix, ".png"), width, height, dpi)
  invisible(prefix)
}

#### ---------- 6. 批量加载 ----------
## 用法：source("R/load_all.R") 即可载入全部函数
popgen_version <- "0.1.0"
