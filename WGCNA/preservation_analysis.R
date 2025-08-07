setwd("/Users/yf358/Desktop/PR1/round3/06_wgcna/v1")
load("wgcna.rdata")
library(WGCNA)

# preservation analysis: human network as reference
color_pig <- netwk_pig$colors
color_human <- netwk_human$colors
color_mouse <- netwk_mouse$colors
color_marmoset <- netwk_marmoset$colors
color_cattle <- netwk_cattle$colors
color_rat <- netwk_rat$colors
colorList = list(color_pig, color_human, color_mouse, 
                 color_marmoset, color_cattle, color_rat)
names(colorList) <- setLabels
mp = modulePreservation(multiExpr, colorList,
                        referenceNetworks = 2,
                        loadPermutedStatistics = FALSE,
                        nPermutations = 200,
                        networkType = "signed", 
                        verbose = 3)

# generate preservation z-score plot
ref = 1 # Select the human data as reference
plot_list <- list()  # 保存所有的图

for (test in c(1, 3, 4, 5, 6)) {
  species <- setLabels[test]
  
  # 提取统计数据
  statsObs <- cbind(mp$quality$observed[[ref]][[test]][, -1], 
                    mp$preservation$observed[[ref]][[test]][, -1])
  statsZ <- cbind(mp$quality$Z[[ref]][[test]], 
                  mp$preservation$Z[[ref]][[test]][, -1])
  
  # 交叉制表统计
  overlap <- cbind(mp$accuracy$observedCounts[[ref]][[test]], 
                   mp$accuracy$observedFisherPvalues[[ref]][[test]])
  write.csv(overlap, paste0(species, "_preservation.csv"))
  
  # 绘制图表
  input1 <- data.frame(
    moduleSize = statsZ$moduleSize,
    Zsummary.pres = signif(statsZ[, "Zsummary.pres", drop = FALSE], 2)
  )
  input1$module <- rownames(statsZ)
  input1 <- input1[-2,]
  
  # 创建颜色数据帧
  colors_df <- data.frame(
    module = netwk_human$colors,
    mergedColors_human
  )
  
  # 计算每种组合的数量
  color_combination_counts <- colors_df %>%
    group_by(module, mergedColors_human) %>%
    summarise(Count = n(), .groups = 'drop') %>%
    filter(Count > 0) %>%
    arrange(desc(Count))
  
  # 合并数据帧
  input1 <- merge(input1, color_combination_counts, by = "module")
  
  # 创建颜色映射
  color_mapping <- setNames(unique(input1$mergedColors_human), unique(input1$mergedColors_human))
  
  # 绘制图形
  dataset1Plot <- ggplot(input1, aes(x = Count, y = Zsummary.pres, fill = mergedColors_human)) +
    geom_point(shape = 21, size = 6) +
    theme_classic() +
    scale_fill_manual(values = color_mapping) +
    labs(y = "Preservation (z-summary)", x = "Number of genes in each module",
         title = paste0(tools::toTitleCase(species), " (Test) vs. Human (Ref)")) + 
    theme(plot.title = element_text(hjust = 0.01, size = 20, face = "bold"),
          plot.margin = unit(c(1, 1, 1, 1), "cm"),
          legend.position = "none",
          axis.title.x = element_text(size = 18),  # 增加 x 轴标题字体大小
          axis.title.y = element_text(size = 18),  # 增加 y 轴标题字体大小
          axis.text.x = element_text(size = 14),   # 增加 x 轴刻度标签字体大小
          axis.text.y = element_text(size = 14)) +    # 增加 y 轴刻度标签字体大小) + # 调整边距
    geom_hline(yintercept = 2, linetype = "dashed", color = "red", size = 1) +
    geom_hline(yintercept = 10, linetype = "dashed", color = "black", size = 1)
  
  plot_list[[length(plot_list) + 1]] <- dataset1Plot  # 保存图形
}

# 创建保存 pdf，并将图表排列在一起
pdf(paste0("preservation_plot_v2.pdf"), width = 16, height = 12)
grid.arrange(grobs = plot_list, nrow = 2, ncol = 3)  # 3列2行的布局
dev.off()
```

# Sankey Plot showing preservation: PIG
### Create circular sankey plot for each pair ###
ref = 1
# PIG
test = 1 #test in c(1, 3, 4, 5, 6) # 1pig, 3mouse, 4marmoset, 5cattle, 6rat
overlap_pig <- cbind(mp$accuracy$observedCounts[[ref]][[test]], 
                     mp$accuracy$observedFisherPvalues[[ref]][[test]])
overlap_pig_p <- overlap_pig[,c(39:ncol(overlap_pig))]
rownames(overlap_pig_p) <- paste0("P", "_",labels2colors(0:(nrow(overlap_pig_p)-1)))
colnames(overlap_pig_p) <- paste0("H", "_",labels2colors(0:(ncol(overlap_pig_p)-1)))
mat <- overlap_pig_p
mat <- -log(mat)
mat <- t(mat)

# 创建分组变量，扇形标签作为命名向量的名称
nm <- unique(unlist(dimnames(mat)))
group <- structure(gsub("^(.).*$", "\\1", nm), names = nm)
table(group)

grid.col <- structure(
  c(labels2colors(0:(nrow(mat)-1)),labels2colors(0:(ncol(mat)-1))),
  names = c(rownames(mat), colnames(mat))
)

# draw sankey plot
circos.clear()
# 最外层添加一个空白轨迹
png("pig_human.png", width = 3000, height = 3000, res = 300)
par(mai = c(0.5, 0.5, 0.5, 0.5))
circos.par(start.degree = 0)
chordDiagram(
  mat, group = group, 
  grid.col = grid.col, 
  annotationTrack = c("grid"), big.gap = 0.5, small.gap = 1.2,
  preAllocateTracks = list(
    track.height = mm_h(6),
    track.margin = c(mm_h(4), 0)
  )
)

# 定义一个计数器，用于控制漂移方向
shift_direction <- 1
circos.track(
  track.index = 2, 
  panel.fun = function(x, y) {
    sector.index = get.cell.meta.data("sector.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    
    # 根据计数器设置y轴的漂移量
    y_offset = ifelse(shift_direction == 1, 1.5, -1.5)  # 向上漂移3或向下漂移3
    x_offset = 0  # 保持x轴位置不变，或根据需要设置x轴漂移量
    
    # circos.text(
    #   mean(xlim) + x_offset,  # x轴位置
    #   mean(ylim) + y_offset,  # y轴位置加上漂移量
    #   sector.index, cex = 0.6,
    #   niceFacing = TRUE
    # )
    
    # 切换漂移方向
    shift_direction <<- -shift_direction
  }, 
  bg.border = NA,
)

# 高亮最外层的扇形格子
highlight.sector(
  rownames(mat), track.index = 1, col = "#b2df8a", 
  text = "HUMAN", cex = 0.9, text.col = "white", 
  niceFacing = TRUE
)
highlight.sector(
  colnames(mat), track.index = 1, col = "#9a96c6", 
  text = "PIG", cex = 0.9, text.col = "white", 
  niceFacing = TRUE
)
dev.off()
circos.clear()
```

# Sankey Plot showing preservation: MOUSE
```{r}
# MOUSE
test = 3 #test in c(1, 3, 4, 5, 6) # 1pig, 3mouse, 4marmoset, 5cattle, 6rat
overlap_mouse <- cbind(mp$accuracy$observedCounts[[ref]][[test]], 
                       mp$accuracy$observedFisherPvalues[[ref]][[test]])
overlap_mouse_p <- overlap_mouse[,c(39:ncol(overlap_mouse))]
rownames(overlap_mouse_p) <- paste0("M", "_",labels2colors(0:(nrow(overlap_mouse_p)-1)))
colnames(overlap_mouse_p) <- paste0("H", "_",labels2colors(0:(ncol(overlap_mouse_p)-1)))
mat <- overlap_mouse_p
mat <- -log(mat)
mat <- t(mat)


# 创建分组变量，扇形标签作为命名向量的名称
nm <- unique(unlist(dimnames(mat)))
group <- structure(gsub("^(.).*$", "\\1", nm), names = nm)
table(group)

grid.col <- structure(
  c(labels2colors(0:(nrow(mat)-1)),labels2colors(0:(ncol(mat)-1))),
  names = c(rownames(mat), colnames(mat))
)

library(circlize)
circos.clear()

# 最外层添加一个空白轨迹
png("mouse_human.png", width = 3000, height = 3000, res = 300)
par(mai = c(0.5, 0.5, 0.5, 0.5))
circos.par(start.degree = 0)
chordDiagram(
  mat, group = group, 
  grid.col = grid.col, 
  annotationTrack = c("grid"), big.gap = 0.5, small.gap = 1.2,
  preAllocateTracks = list(
    track.height = mm_h(6),
    track.margin = c(mm_h(4), 0)
  )
)

# 定义一个计数器，用于控制漂移方向
shift_direction <- 1
circos.track(
  track.index = 2, 
  panel.fun = function(x, y) {
    sector.index = get.cell.meta.data("sector.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    
    # 根据计数器设置y轴的漂移量
    y_offset = ifelse(shift_direction == 1, 1.5, -1.5)  # 向上漂移3或向下漂移3
    x_offset = 0  # 保持x轴位置不变，或根据需要设置x轴漂移量
    
    # circos.text(
    #   mean(xlim) + x_offset,  # x轴位置
    #   mean(ylim) + y_offset,  # y轴位置加上漂移量
    #   sector.index, cex = 0.6,
    #   niceFacing = TRUE
    # )
    
    # 切换漂移方向
    shift_direction <<- -shift_direction
  }, 
  bg.border = NA,
)

# 高亮最外层的扇形格子
highlight.sector(
  rownames(mat), track.index = 1, col = "#b2df8a", 
  text = "HUMAN", cex = 0.9, text.col = "white", 
  niceFacing = TRUE
)
highlight.sector(
  colnames(mat), track.index = 1, col = "#9a96c6", 
  text = "MOUSE", cex = 0.9, text.col = "white", 
  niceFacing = TRUE
)
dev.off()
circos.clear()
```

# Sankey Plot showing preservation: MARMOSET
```{r}
test = 4 #test in c(1, 3, 4, 5, 6) # 1pig, 3mouse, 4marmoset, 5cattle, 6rat
overlap_marmoset <- cbind(mp$accuracy$observedCounts[[ref]][[test]], 
                          mp$accuracy$observedFisherPvalues[[ref]][[test]])
overlap_marmoset_p <- overlap_marmoset[,c(39:ncol(overlap_marmoset))]
rownames(overlap_marmoset_p) <- paste0("M", "_",labels2colors(0:(nrow(overlap_marmoset_p)-1)))
colnames(overlap_marmoset_p) <- paste0("H", "_",labels2colors(0:(ncol(overlap_marmoset_p)-1)))
mat <- overlap_marmoset_p
mat <- -log(mat)
mat <- t(mat)


# 创建分组变量，扇形标签作为命名向量的名称
nm <- unique(unlist(dimnames(mat)))
group <- structure(gsub("^(.).*$", "\\1", nm), names = nm)
table(group)

grid.col <- structure(
  c(labels2colors(0:(nrow(mat)-1)),labels2colors(0:(ncol(mat)-1))),
  names = c(rownames(mat), colnames(mat))
)

library(circlize)
circos.clear()

# 最外层添加一个空白轨迹
png("marmoset_human.png", width = 3000, height = 3000, res = 300)
par(mai = c(0.5, 0.5, 0.5, 0.5))
circos.par(start.degree = 0)
chordDiagram(
  mat, group = group, 
  grid.col = grid.col, 
  annotationTrack = c("grid"), big.gap = 0.5, small.gap = 1.2,
  preAllocateTracks = list(
    track.height = mm_h(6),
    track.margin = c(mm_h(4), 0)
  )
)

# 定义一个计数器，用于控制漂移方向
shift_direction <- 1
circos.track(
  track.index = 2, 
  panel.fun = function(x, y) {
    sector.index = get.cell.meta.data("sector.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    
    # 根据计数器设置y轴的漂移量
    y_offset = ifelse(shift_direction == 1, 1.5, -1.5)  # 向上漂移3或向下漂移3
    x_offset = 0  # 保持x轴位置不变，或根据需要设置x轴漂移量
    
    # circos.text(
    #   mean(xlim) + x_offset,  # x轴位置
    #   mean(ylim) + y_offset,  # y轴位置加上漂移量
    #   sector.index, cex = 0.6,
    #   niceFacing = TRUE
    # )
    
    # 切换漂移方向
    shift_direction <<- -shift_direction
  }, 
  bg.border = NA,
)

# 高亮最外层的扇形格子
highlight.sector(
  rownames(mat), track.index = 1, col = "#b2df8a", 
  text = "HUMAN", cex = 0.9, text.col = "white", 
  niceFacing = TRUE
)
highlight.sector(
  colnames(mat), track.index = 1, col = "#9a96c6", 
  text = "MARMOSET", cex = 0.9, text.col = "white", 
  niceFacing = TRUE
)
dev.off()
circos.clear()

```

# Sankey Plot showing preservation: CATTLE
```{r}
test = 5 #test in c(1, 3, 4, 5, 6) # 1pig, 3mouse, 4marmoset, 5cattle, 6rat
overlap_cattle <- cbind(mp$accuracy$observedCounts[[ref]][[test]], 
                        mp$accuracy$observedFisherPvalues[[ref]][[test]])
overlap_cattle_p <- overlap_cattle[,c(39:ncol(overlap_cattle))]
rownames(overlap_cattle_p) <- paste0("C", "_",labels2colors(0:(nrow(overlap_cattle_p)-1)))
colnames(overlap_cattle_p) <- paste0("H", "_",labels2colors(0:(ncol(overlap_cattle_p)-1)))
mat <- overlap_cattle_p
mat <- -log(mat)
mat <- t(mat)


# 创建分组变量，扇形标签作为命名向量的名称
nm <- unique(unlist(dimnames(mat)))
group <- structure(gsub("^(.).*$", "\\1", nm), names = nm)
table(group)

grid.col <- structure(
  c(labels2colors(0:(nrow(mat)-1)),labels2colors(0:(ncol(mat)-1))),
  names = c(rownames(mat), colnames(mat))
)

circos.clear()
# 最外层添加一个空白轨迹
png("cattle_human.png", width = 3000, height = 3000, res = 300)
par(mai = c(0.5, 0.5, 0.5, 0.5))
circos.par(start.degree = 180)
chordDiagram(
  mat, group = group, 
  grid.col = grid.col, 
  annotationTrack = c("grid"), big.gap = 0.5, small.gap = 1.2,
  preAllocateTracks = list(
    track.height = mm_h(6),
    track.margin = c(mm_h(4), 0)
  )
)

# 定义一个计数器，用于控制漂移方向
shift_direction <- 1
circos.track(
  track.index = 2, 
  panel.fun = function(x, y) {
    sector.index = get.cell.meta.data("sector.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    
    # 根据计数器设置y轴的漂移量
    y_offset = ifelse(shift_direction == 1, 1.5, -1.5)  # 向上漂移3或向下漂移3
    x_offset = 0  # 保持x轴位置不变，或根据需要设置x轴漂移量
    
    # circos.text(
    #   mean(xlim) + x_offset,  # x轴位置
    #   mean(ylim) + y_offset,  # y轴位置加上漂移量
    #   sector.index, cex = 0.6,
    #   niceFacing = TRUE
    # )
    
    # 切换漂移方向
    shift_direction <<- -shift_direction
  }, 
  bg.border = NA,
)

# 高亮最外层的扇形格子
highlight.sector(
  rownames(mat), track.index = 1, col = "#b2df8a", 
  text = "HUMAN", cex = 0.9, text.col = "white", 
  niceFacing = TRUE
)
highlight.sector(
  colnames(mat), track.index = 1, col = "#9a96c6", 
  text = "CATTLE", cex = 0.9, text.col = "white", 
  niceFacing = TRUE
)

dev.off()
circos.clear()

```
# Sankey Plot showing preservation: RAT
```{r}
test = 6 #test in c(1, 3, 4, 5, 6) # 1pig, 3mouse, 4marmoset, 5cattle, 6rat
overlap_rat <- cbind(mp$accuracy$observedCounts[[ref]][[test]], 
                     mp$accuracy$observedFisherPvalues[[ref]][[test]])
overlap_rat_p <- overlap_rat[,c(39:ncol(overlap_rat))]
rownames(overlap_rat_p) <- paste0("R", "_",labels2colors(0:(nrow(overlap_rat_p)-1)))
colnames(overlap_rat_p) <- paste0("H", "_",labels2colors(0:(ncol(overlap_rat_p)-1)))
mat <- overlap_rat_p
mat <- -log(mat)
mat <- t(mat)


# 创建分组变量，扇形标签作为命名向量的名称
nm <- unique(unlist(dimnames(mat)))
group <- structure(gsub("^(.).*$", "\\1", nm), names = nm)
table(group)

grid.col <- structure(
  c(labels2colors(0:(nrow(mat)-1)),labels2colors(0:(ncol(mat)-1))),
  names = c(rownames(mat), colnames(mat))
)
grid.col

library(circlize)
circos.clear()

# 最外层添加一个空白轨迹
png("rat_human.png", width = 3000, height = 3000, res = 300)
par(mai = c(0.5, 0.5, 0.5, 0.5))
circos.par(start.degree = 0)
chordDiagram(
  mat, group = group, 
  grid.col = grid.col, 
  annotationTrack = c("grid"), big.gap = 0.5, small.gap = 1.2,
  preAllocateTracks = list(
    track.height = mm_h(6),
    track.margin = c(mm_h(4), 0)
  )
)

# 定义一个计数器，用于控制漂移方向
shift_direction <- 1
circos.track(
  track.index = 2, 
  panel.fun = function(x, y) {
    sector.index = get.cell.meta.data("sector.index")
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    
    # 根据计数器设置y轴的漂移量
    y_offset = ifelse(shift_direction == 1, 1.5, -1.5)  # 向上漂移3或向下漂移3
    x_offset = 0  # 保持x轴位置不变，或根据需要设置x轴漂移量
    
    # circos.text(
    #   mean(xlim) + x_offset,  # x轴位置
    #   mean(ylim) + y_offset,  # y轴位置加上漂移量
    #   sector.index, cex = 0.6,
    #   niceFacing = TRUE
    # )
    
    # 切换漂移方向
    shift_direction <<- -shift_direction
  }, 
  bg.border = NA,
)

# 高亮最外层的扇形格子
highlight.sector(
  rownames(mat), track.index = 1, col = "#b2df8a", 
  text = "HUMAN", cex = 0.9, text.col = "white", 
  niceFacing = TRUE
)
highlight.sector(
  colnames(mat), track.index = 1, col = "#9a96c6", 
  text = "RAT", cex = 0.9, text.col = "white", 
  niceFacing = TRUE
)

dev.off()
circos.clear()
