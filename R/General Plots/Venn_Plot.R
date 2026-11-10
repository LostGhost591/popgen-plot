#### ============================================================
#### Venn_Plot.R —— 韦恩图（复刻 Venny 2.1 的效果）
####
#### 对标 https://bioinfogp.cnb.csic.es/tools/venny/index.html
####   · 2 / 3 / 4 个集合，每集合一色，半透明叠加
####   · 每个区域显示元素个数
####   · 3 集合是标准三角布局，4 集合是标准对称四椭圆（不是按面积缩放的欧拉图）
####   · 线宽、字号、字体族可调
####
#### 实现用的是 VennDiagram 包，关键是 euler.d = FALSE + scaled = FALSE：
####   默认 TRUE 会按集合大小画成欧拉图，形状会变形，就不是 Venny 那个样子了。
####
#### 用法：
####   Venn_Plot(list(`List 1` = genesA, `List 2` = genesB, `List 3` = genesC),
####            out_file = "out/venn.png")
####
#### 注意：返回的是 grid grob，不是 ggplot 对象，所以用不了 save_popgen()。
#### ============================================================

## 内部工具：二集合韦恩图，自己用 grid 画
## 不用 VennDiagram 的 draw.pairwise.venn，因为它在 cat.dist 传向量时
## 会把 area1 画到右边、fill 顺序跟着错位，两个集合名还会叠在正上方。
.venn2Grob <- function(sets, colors, alpha, cex, cat_cex, lwd, font_family) {
  nA  <- length(sets[[1]])
  nB  <- length(sets[[2]])
  nAB <- length(intersect(sets[[1]], sets[[2]]))

  r  <- 0.30
  cx <- c(0.5 - 0.17, 0.5 + 0.17)

  circle <- function(col, x) {
    grid::circleGrob(x = grid::unit(x, "npc"), y = grid::unit(0.5, "npc"),
                     r = grid::unit(r, "snpc"),
                     gp = grid::gpar(fill = col, alpha = alpha,
                                     col = "black", lwd = lwd))
  }
  txt <- function(lbl, x, y, size, face = "plain") {
    grid::textGrob(lbl, x = grid::unit(x, "npc"), y = grid::unit(y, "npc"),
                   gp = grid::gpar(cex = size, fontfamily = font_family,
                                   fontface = face))
  }

  grid::grobTree(
    circle(colors[1], cx[1]),
    circle(colors[2], cx[2]),
    txt(nA - nAB, cx[1] - 0.13, 0.5, cex),
    txt(nB - nAB, cx[2] + 0.13, 0.5, cex),
    txt(nAB,      0.5,          0.5, cex),
    txt(names(sets)[1], cx[1] - 0.02, 0.5 + r + 0.05, cat_cex),
    txt(names(sets)[2], cx[2] + 0.02, 0.5 + r + 0.05, cat_cex)
  )
}

## 内部工具：枚举所有非空区域，返回「区域 -> 元素」的 list
## 区域名形如 "List 1 & List 2"，表示恰好只在这几个集合里出现的元素
vennRegions <- function(sets) {
  n   <- length(sets)
  nms <- names(sets)
  out <- list()
  for (mask in seq_len(2^n - 1)) {
    in_idx  <- which(bitwAnd(mask, bitwShiftL(1L, 0:(n - 1L))) > 0)
    out_idx <- setdiff(seq_len(n), in_idx)
    common  <- Reduce(intersect, sets[in_idx])
    if (length(out_idx)) {
      common <- setdiff(common, unlist(sets[out_idx], use.names = FALSE))
    }
    out[[paste(nms[in_idx], collapse = " & ")]] <- sort(unique(common))
  }
  out
}

Venn_Plot <- function(sets,                     out_file    = NULL,
                     colors      = NULL,       # NULL = PALETTE_VENN（红黄绿蓝）
                     alpha       = 0.4,        # 填充透明度，Venny 默认偏淡
                     title       = NULL,
                     cex         = 1.6,        # 区域数字字号
                     cat_cex     = 1.6,        # 集合名字号
                     line_width  = 2,          # 圆/椭圆边框线宽
                     font_family = "sans",     # "sans" / "serif" / "mono"
                     detail_out  = NULL,       # 可选：把每个区域的元素明细写成 tsv
                     width  = 8,
                     height = 8,
                     dpi    = 300) {

  if (!is.list(sets) || is.null(names(sets)) || any(!nzchar(names(sets)))) {
    stop("sets 必须是带名字的 list，例如：\n",
         "  list(`List 1` = c(\"a\",\"b\"), `List 2` = c(\"b\",\"c\"))")
  }
  n <- length(sets)
  if (n < 2 || n > 4) {
    stop("只支持 2~4 个集合（Venny 的上限也是 4），当前 ", n, " 个")
  }
  if (!requireNamespace("VennDiagram", quietly = TRUE)) {
    stop("缺少包 VennDiagram，请先 install.packages('VennDiagram')")
  }

  sets <- lapply(sets, function(x) unique(as.character(x)))

  #### ---- 配色 ----
  if (is.null(colors)) colors <- PALETTE_VENN[seq_len(n)]
  if (length(colors) < n) stop("colors 至少需要 ", n, " 个颜色")
  colors <- colors[seq_len(n)]

  #### ---- 交集辅助 ----
  I <- function(idx) Reduce(intersect, sets[idx])

  #### ---- 画图 ----
  ## 一律 euler.d = FALSE, scaled = FALSE：这是 Venny 那种标准韦恩布局。
  ## 默认 TRUE 会按集合大小画成欧拉图，圆会变形、位置会错乱。
  ##
  ## 这里刻意不用 do.call + c() 拼参数 —— 那样传进去后 2 集合的
  ## fill / category 顺序会错位、cat.pos 也不生效。

  if (n == 2) {
    g <- .venn2Grob(sets, colors, alpha, cex, cat_cex, line_width, font_family)

  } else if (n == 3) {
    g <- VennDiagram::draw.triple.venn(
      area1    = length(sets[[1]]),
      area2    = length(sets[[2]]),
      area3    = length(sets[[3]]),
      n12      = length(I(c(1, 2))),
      n23      = length(I(c(2, 3))),
      n13      = length(I(c(1, 3))),
      n123     = length(I(1:3)),
      category = names(sets),
      fill     = colors,
      alpha    = rep(alpha, 3),
      cex      = cex,
      cat.cex  = cat_cex,
      fontfamily = font_family,
      cat.fontfamily = font_family,
      lwd      = line_width,
      euler.d  = FALSE, scaled = FALSE, ind = FALSE,
      ## 默认 0.05 时正下方那个集合名会贴住圆边，稍微拉开一点
      cat.dist = c(0.07, 0.07, 0.07),
      margin   = 0.1
    )

  } else {
    g <- VennDiagram::draw.quad.venn(
      area1    = length(sets[[1]]),
      area2    = length(sets[[2]]),
      area3    = length(sets[[3]]),
      area4    = length(sets[[4]]),
      n12      = length(I(c(1, 2))),
      n13      = length(I(c(1, 3))),
      n14      = length(I(c(1, 4))),
      n23      = length(I(c(2, 3))),
      n24      = length(I(c(2, 4))),
      n34      = length(I(c(3, 4))),
      n123     = length(I(c(1, 2, 3))),
      n124     = length(I(c(1, 2, 4))),
      n134     = length(I(c(1, 3, 4))),
      n234     = length(I(c(2, 3, 4))),
      n1234    = length(I(1:4)),
      category = names(sets),
      fill     = colors,
      alpha    = rep(alpha, 4),
      cex      = cex,
      cat.cex  = cat_cex,
      fontfamily = font_family,
      cat.fontfamily = font_family,
      lwd      = line_width,
      euler.d  = FALSE, scaled = FALSE, ind = FALSE,
      margin   = 0.12
    )
  }

  #### ---- 标题（VennDiagram 本身不支持，用 grid 叠一层）----
  if (!is.null(title)) {
    g <- grid::grobTree(
      g,
      grid::textGrob(title,
                     y   = grid::unit(0.97, "npc"),
                     gp  = grid::gpar(cex = 1.5, fontface = "bold",
                                      fontfamily = font_family))
    )
  }

  #### ---- 导出 ----
  if (!is.null(out_file)) {
    dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
    if (grepl("\\.pdf$", out_file, ignore.case = TRUE)) {
      grDevices::pdf(out_file, width = width, height = height)
    } else {
      grDevices::png(out_file, width = width, height = height,
                     units = "in", res = dpi)
    }
    grid::grid.newpage()
    grid::grid.draw(g)
    grDevices::dev.off()
  }

  #### ---- 区域明细 ----
  regions <- vennRegions(sets)
  counts  <- vapply(regions, length, integer(1))

  if (!is.null(detail_out)) {
    dir.create(dirname(detail_out), showWarnings = FALSE, recursive = TRUE)
    df <- data.frame(
      region  = rep(names(regions), counts),
      element = unlist(regions, use.names = FALSE),
      stringsAsFactors = FALSE
    )
    utils::write.table(df, detail_out, sep = "\t",
                       quote = FALSE, row.names = FALSE)
  }

  invisible(list(plot = g, counts = counts, regions = regions))
}
