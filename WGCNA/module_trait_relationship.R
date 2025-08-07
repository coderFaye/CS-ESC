####### updated module-trait relationship plots ######## 

# PIG
setwd("/Users/yf358/Desktop/PR1/round3/06_wgcna/v1")
load("wgcna.rdata")
library(WGCNA)
set <- 1
species <- "pig"
expr_normalized <- multiExpr[[set]]$data
# create plot of module-trait correlations
mergedColors = labels2colors(netwk_pig$colors)
moduleColors=mergedColors


if(T){
  nGenes = ncol(expr_normalized)  # Number of genes
  nSamples = nrow(expr_normalized)  # Number of samples
  
  # Calculate module eigengenes
  MES0 <- moduleEigengenes(expr_normalized, moduleColors)$eigengenes
  MEs = orderMEs(MES0)  # Sort module eigengenes
  
  # Compute module-trait correlations and p-values
  moduleTraitCor <- cor(MEs, Traits[[set]]$data, use = "p")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
  
  # Create a formatted text matrix for displaying correlations and p-values
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
  dim(textMatrix) = dim(moduleTraitCor)  # Adjust dimensions of the text matrix
  
  # Set the width and height dynamically based on the size of the heatmap
  heatmap_width <- 2 + 1.7 * ncol(moduleTraitCor)  # Customize as needed based on number of traits
  heatmap_height <- 2 + 0.5 * nrow(moduleTraitCor)  # Customize based on number of modules
  
  # Save the plot to a PDF file
  pdf("V2_module-trait-relationships_pig_adjusted_font.pdf", 
      width = heatmap_width, height = heatmap_height)
  
  # Adjust margins for a cleaner layout
  par(mar = c(10, 12, 4, 2))
  
  # Use medium red and blue contrast
  labeledHeatmap(
    Matrix = moduleTraitCor, 
    xLabels = names(Traits[[set]]$data),
    yLabels = names(MEs),
    ySymbols = names(MEs),
    colorLabels = FALSE,
    colors = colorRampPalette(c("#87CEEB", "#FFFFFF", "#FF4C4C"))(50),  # Medium blue to soft red contrast
    textMatrix = textMatrix,
    setStdMargins = FALSE,
    cex.text = 0.9,  # Adjust text size for the correlation values
    cex.lab = 1.0,  # Larger font size for module labels (you can increase this further if needed)
    zlim = c(-1, 1),  # Correlation range from -1 to 1
    main = paste("Module-Trait Relationships in", species),
    cex.main = 1.0,  # Larger title font size
    borderColor = "grey85",  # Light gray border for separation
    frame = FALSE,  # Remove heavy outer borders for a sleek look
    grid = TRUE  # Add gridlines for clarity
  )
  
  # Close the file
  dev.off()
}

# Set the width and height dynamically based on the size of the heatmap
heatmap_width <- 2 + 1.5 * ncol(moduleTraitCor)  # Based on number of traits
heatmap_height <- 2 + 0.5 * nrow(moduleTraitCor)  # Based on number of modules

# Save the plot to a PNG file with the same aspect ratio
png("V2_module-trait-relationships_pig_adjusted_font.png", 
    width = heatmap_width, height = heatmap_height, units = "in", res = 300)  # High resolution

# Adjust margins for a cleaner layout
par(mar = c(10, 12, 4, 2))

# Use medium red and blue contrast in labeledHeatmap
labeledHeatmap(
  Matrix = moduleTraitCor, 
  xLabels = names(Traits[[set]]$data),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = colorRampPalette(c("#87CEEB", "#FFFFFF", "#FF4C4C"))(50),  # Medium blue to soft red contrast
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.9,  # Adjust text size for the correlation values
  cex.lab = 1.0,  # Font size for module labels
  zlim = c(-1, 1),  # Correlation range from -1 to 1
  main = paste("Module-Trait Relationships in", species),
  cex.main = 1.0,  # Title font size
  borderColor = "grey85",  # Light gray border
  frame = FALSE,  # Remove outer borders
  grid = TRUE  # Add gridlines
)

# Close the PNG device
dev.off()


# setLabels = c("pig", "human", "mouse", "marmoset", "cattle","rat")
set <- 2
species <- "human"
expr_normalized <- multiExpr[[set]]$data
# create plot of module-trait correlations
mergedColors = labels2colors(netwk_human$colors)
moduleColors=mergedColors

if(T){
  nGenes = ncol(expr_normalized)  # Number of genes
  nSamples = nrow(expr_normalized)  # Number of samples
  
  # Calculate module eigengenes
  MES0 <- moduleEigengenes(expr_normalized, moduleColors)$eigengenes
  MEs = orderMEs(MES0)  # Sort module eigengenes
  
  # Compute module-trait correlations and p-values
  moduleTraitCor <- cor(MEs, Traits[[set]]$data, use = "p")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
  
  # Create a formatted text matrix for displaying correlations and p-values
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
  dim(textMatrix) = dim(moduleTraitCor)  # Adjust dimensions of the text matrix
  
  # Set the width and height dynamically based on the size of the heatmap
  heatmap_width <- 2 + 2.5 * ncol(moduleTraitCor)  # Customize as needed based on number of traits
  heatmap_height <- 2 + 0.5 * nrow(moduleTraitCor)  # Customize based on number of modules
  
  # Save the plot to a PDF file
  pdf("V2_module-trait-relationships_human_adjusted_font.pdf", 
      width = heatmap_width, height = heatmap_height)
  
  # Adjust margins for a cleaner layout
  par(mar = c(12, 12, 8, 6))
  
  # Use medium red and blue contrast
  labeledHeatmap(
    Matrix = moduleTraitCor, 
    xLabels = names(Traits[[set]]$data),
    yLabels = names(MEs),
    ySymbols = names(MEs),
    colorLabels = FALSE,
    colors = colorRampPalette(c("#87CEEB", "#FFFFFF", "#FF4C4C"))(50),  # Medium blue to soft red contrast
    textMatrix = textMatrix,
    setStdMargins = FALSE,
    cex.text = 1.2,  # Adjust text size for the correlation values
    cex.lab = 1.5,  # Larger font size for module labels (you can increase this further if needed)
    zlim = c(-1, 1),  # Correlation range from -1 to 1
    main = paste("Module-Trait Relationships in", species),
    cex.main = 1.5,  # Larger title font size
    borderColor = "grey85",  # Light gray border for separation
    frame = FALSE,  # Remove heavy outer borders for a sleek look
    grid = TRUE  # Add gridlines for clarity
  )
  
  # Close the file
  dev.off()
}

# Set the width and height dynamically based on the size of the heatmap
heatmap_width <- 2 + 2.5 * ncol(moduleTraitCor)  # Adjust based on number of traits
heatmap_height <- 2 + 0.5 * nrow(moduleTraitCor)  # Adjust based on number of modules

# Save the plot to a PNG file with the same aspect ratio
png("V2_module-trait-relationships_human_adjusted_font.png", 
    width = heatmap_width, height = heatmap_height, units = "in", res = 300)  # High resolution

# Adjust margins for a cleaner layout
par(mar = c(12, 12, 8, 6))

# Use medium red and blue contrast in labeledHeatmap
labeledHeatmap(
  Matrix = moduleTraitCor, 
  xLabels = names(Traits[[set]]$data),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = colorRampPalette(c("#87CEEB", "#FFFFFF", "#FF4C4C"))(50),  # Medium blue to soft red contrast
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 1.2,  # Adjust text size for the correlation values
  cex.lab = 1.5,  # Larger font size for module labels
  zlim = c(-1, 1),  # Correlation range from -1 to 1
  main = paste("Module-Trait Relationships in", species),
  cex.main = 1.5,  # Title font size
  borderColor = "grey85",  # Light gray border
  frame = FALSE,  # Remove outer borders
  grid = TRUE  # Add gridlines
)

# Close the PNG device
dev.off()


# setLabels = c("pig", "human", "mouse", "marmoset", "cattle","rat")
set <- 3
species <- "mouse"
expr_normalized <- multiExpr[[set]]$data
# create plot of module-trait correlations
mergedColors = labels2colors(netwk_mouse$colors)
moduleColors=mergedColors

if(T){
  nGenes = ncol(expr_normalized)  # Number of genes
  nSamples = nrow(expr_normalized)  # Number of samples
  
  # Calculate module eigengenes
  MES0 <- moduleEigengenes(expr_normalized, moduleColors)$eigengenes
  MEs = orderMEs(MES0)  # Sort module eigengenes
  
  # Compute module-trait correlations and p-values
  moduleTraitCor <- cor(MEs, Traits[[set]]$data, use = "p")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
  
  # Create a formatted text matrix for displaying correlations and p-values
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
  dim(textMatrix) = dim(moduleTraitCor)  # Adjust dimensions of the text matrix
  
  # Set the width and height dynamically based on the size of the heatmap
  heatmap_width <- 2 + 2.5 * ncol(moduleTraitCor)  # Customize as needed based on number of traits
  heatmap_height <- 2 + 0.8 * nrow(moduleTraitCor)  # Customize based on number of modules
  
  # Save the plot to a PDF file
  pdf("V2_module-trait-relationships_mouse_adjusted_font.pdf", 
      width = heatmap_width, height = heatmap_height)
  
  # Adjust margins for a cleaner layout
  par(mar = c(12, 12, 8, 6))
  
  # Use medium red and blue contrast
  labeledHeatmap(
    Matrix = moduleTraitCor, 
    xLabels = names(Traits[[set]]$data),
    yLabels = names(MEs),
    ySymbols = names(MEs),
    colorLabels = FALSE,
    colors = colorRampPalette(c("#87CEEB", "#FFFFFF", "#FF4C4C"))(50),  # Medium blue to soft red contrast
    textMatrix = textMatrix,
    setStdMargins = FALSE,
    cex.text = 1.5,  # Adjust text size for the correlation values
    cex.lab = 1.5,  # Larger font size for module labels (you can increase this further if needed)
    zlim = c(-1, 1),  # Correlation range from -1 to 1
    main = paste("Module-Trait Relationships in", species),
    cex.main = 1.5,  # Larger title font size
    borderColor = "grey85",  # Light gray border for separation
    frame = FALSE,  # Remove heavy outer borders for a sleek look
    grid = TRUE  # Add gridlines for clarity
  )
  
  # Close the file
  dev.off()
}

# Set the width and height dynamically based on the size of the heatmap
heatmap_width <- 2 + 2.8 * ncol(moduleTraitCor)  # Based on number of traits
heatmap_height <- 2 + 0.8 * nrow(moduleTraitCor)  # Based on number of modules

# Save the plot to a PNG file with the same aspect ratio
png("V2_module-trait-relationships_mouse_adjusted_font.png", 
    width = heatmap_width, height = heatmap_height, units = "in", res = 300)  # High resolution

# Adjust margins for a cleaner layout
par(mar = c(12, 12, 8, 6))

# Use medium red and blue contrast in labeledHeatmap
labeledHeatmap(
  Matrix = moduleTraitCor, 
  xLabels = names(Traits[[set]]$data),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = colorRampPalette(c("#87CEEB", "#FFFFFF", "#FF4C4C"))(50),  # Medium blue to soft red contrast
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 1.5,  # Adjust text size for the correlation values
  cex.lab = 1.5,  # Larger font size for module labels
  zlim = c(-1, 1),  # Correlation range from -1 to 1
  main = paste("Module-Trait Relationships in", species),
  cex.main = 1.5,  # Title font size
  borderColor = "grey85",  # Light gray border
  frame = FALSE,  # Remove outer borders
  grid = TRUE  # Add gridlines
)

# Close the PNG device
dev.off()



# setLabels = c("pig", "human", "mouse", "marmoset", "cattle","rat")
set <- 4
species <- "marmoset"
expr_normalized <- multiExpr[[set]]$data
# create plot of module-trait correlations
mergedColors = labels2colors(netwk_marmoset$colors)
moduleColors=mergedColors


if(T){
  nGenes = ncol(expr_normalized)  # Number of genes
  nSamples = nrow(expr_normalized)  # Number of samples
  
  # Calculate module eigengenes
  MES0 <- moduleEigengenes(expr_normalized, moduleColors)$eigengenes
  MEs = orderMEs(MES0)  # Sort module eigengenes
  
  
  # Compute module-trait correlations and p-values
  moduleTraitCor <- cor(MEs, Traits[[set]]$data[, colSums(Traits[[set]]$data != 0) > 0], use = "p")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
  
  # Create a formatted text matrix for displaying correlations and p-values
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
  dim(textMatrix) = dim(moduleTraitCor)  # Adjust dimensions of the text matrix
  
  # Set the width and height dynamically based on the size of the heatmap
  heatmap_width <- 2 + 2.5 * ncol(moduleTraitCor)  # Customize as needed based on number of traits
  heatmap_height <- 2 + 0.8 * nrow(moduleTraitCor)  # Customize based on number of modules
  
  # Save the plot to a PDF file
  pdf("V2_module-trait-relationships_marmoset_adjusted_font.pdf", 
      width = heatmap_width, height = heatmap_height)
  
  # Adjust margins for a cleaner layout
  par(mar = c(12, 12, 8, 6))
  
  # Use medium red and blue contrast
  labeledHeatmap(
    Matrix = moduleTraitCor, 
    xLabels = names(Traits[[set]]$data[, colSums(Traits[[set]]$data != 0) > 0]),
    yLabels = names(MEs),
    ySymbols = names(MEs),
    colorLabels = FALSE,
    colors = colorRampPalette(c("#87CEEB", "#FFFFFF", "#FF4C4C"))(50),  # Medium blue to soft red contrast
    textMatrix = textMatrix,
    setStdMargins = FALSE,
    cex.text = 1.5,  # Adjust text size for the correlation values
    cex.lab = 1.5,  # Larger font size for module labels (you can increase this further if needed)
    zlim = c(-1, 1),  # Correlation range from -1 to 1
    main = paste("Module-Trait Relationships in", species),
    cex.main = 1.5,  # Larger title font size
    borderColor = "grey85",  # Light gray border for separation
    frame = FALSE,  # Remove heavy outer borders for a sleek look
    grid = TRUE  # Add gridlines for clarity
  )
  
  # Close the file
  dev.off()
}

# Set the width and height dynamically based on the size of the heatmap
heatmap_width <- 2 + 3.2 * ncol(moduleTraitCor)  # Adjust based on number of traits
heatmap_height <- 2 + 0.8 * nrow(moduleTraitCor)  # Adjust based on number of modules

# Save the plot to a PNG file with the same aspect ratio
png("V2_module-trait-relationships_marmoset_adjusted_font.png", 
    width = heatmap_width, height = heatmap_height, units = "in", res = 300)  # High resolution

# Adjust margins for a cleaner layout
par(mar = c(12, 12, 8, 6))

# Use medium red and blue contrast in labeledHeatmap
labeledHeatmap(
  Matrix = moduleTraitCor, 
  xLabels = names(Traits[[set]]$data[, colSums(Traits[[set]]$data != 0) > 0]),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = colorRampPalette(c("#87CEEB", "#FFFFFF", "#FF4C4C"))(50),  # Medium blue to soft red contrast
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 1.5,  # Adjust text size for the correlation values
  cex.lab = 1.5,  # Larger font size for module labels
  zlim = c(-1, 1),  # Correlation range from -1 to 1
  main = paste("Module-Trait Relationships in", species),
  cex.main = 1.5,  # Larger title font size
  borderColor = "grey85",  # Light gray border
  frame = FALSE,  # Remove outer borders
  grid = TRUE  # Add gridlines
)

# Close the PNG device
dev.off()



# setLabels = c("pig", "human", "mouse", "marmoset", "cattle","rat")
set <- 5
species <- "cattle"
expr_normalized <- multiExpr[[set]]$data
# create plot of module-trait correlations
mergedColors = labels2colors(netwk_cattle$colors)
moduleColors=mergedColors


if(T){
  nGenes = ncol(expr_normalized)  # Number of genes
  nSamples = nrow(expr_normalized)  # Number of samples
  
  # Calculate module eigengenes
  MES0 <- moduleEigengenes(expr_normalized, moduleColors)$eigengenes
  MEs = orderMEs(MES0)  # Sort module eigengenes
  
  
  # Compute module-trait correlations and p-values
  moduleTraitCor <- cor(MEs, Traits[[set]]$data[, colSums(Traits[[set]]$data != 0) > 0], use = "p")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
  
  # Create a formatted text matrix for displaying correlations and p-values
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
  dim(textMatrix) = dim(moduleTraitCor)  # Adjust dimensions of the text matrix
  
  # Set the width and height dynamically based on the size of the heatmap
  heatmap_width <- 2 + 3.5 * ncol(moduleTraitCor)  # Customize as needed based on number of traits
  heatmap_height <- 2 + 0.8 * nrow(moduleTraitCor)  # Customize based on number of modules
  
  # Save the plot to a PDF file
  pdf("V2_module-trait-relationships_cattle_adjusted_font.pdf", 
      width = heatmap_width, height = heatmap_height)
  
  # Adjust margins for a cleaner layout
  par(mar = c(12, 12, 8, 6))
  
  # Use medium red and blue contrast
  labeledHeatmap(
    Matrix = moduleTraitCor, 
    xLabels = names(Traits[[set]]$data[, colSums(Traits[[set]]$data != 0) > 0]),
    yLabels = names(MEs),
    ySymbols = names(MEs),
    colorLabels = FALSE,
    colors = colorRampPalette(c("#87CEEB", "#FFFFFF", "#FF4C4C"))(50),  # Medium blue to soft red contrast
    textMatrix = textMatrix,
    setStdMargins = FALSE,
    cex.text = 1.5,  # Adjust text size for the correlation values
    cex.lab = 1.5,  # Larger font size for module labels (you can increase this further if needed)
    zlim = c(-1, 1),  # Correlation range from -1 to 1
    main = paste("Module-Trait Relationships in", species),
    cex.main = 1.5,  # Larger title font size
    borderColor = "grey85",  # Light gray border for separation
    frame = FALSE,  # Remove heavy outer borders for a sleek look
    grid = TRUE  # Add gridlines for clarity
  )
  
  # Close the file
  dev.off()
}

# Set the width and height dynamically based on the size of the heatmap
heatmap_width <- 2 + 3.5 * ncol(moduleTraitCor)  # Based on number of traits
heatmap_height <- 2 + 0.8 * nrow(moduleTraitCor)  # Based on number of modules
# Save the plot to a PNG file with the same aspect ratio
png("V2_module-trait-relationships_cattle_adjusted_font.png", 
    width = heatmap_width, height = heatmap_height, units = "in", res = 300)  # High resolution

# Adjust margins for a cleaner layout
par(mar = c(12, 12, 8, 6))
# Use medium red and blue contrast in labeledHeatmap
labeledHeatmap(
  Matrix = moduleTraitCor, 
  xLabels = names(Traits[[set]]$data[, colSums(Traits[[set]]$data != 0) > 0]),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = colorRampPalette(c("#87CEEB", "#FFFFFF", "#FF4C4C"))(50),  # Medium blue to soft red contrast
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 1.5,  # Adjust text size for the correlation values
  cex.lab = 1.5,  # Larger font size for module labels
  zlim = c(-1, 1),  # Correlation range from -1 to 1
  main = paste("Module-Trait Relationships in", species),
  cex.main = 1.5,  # Larger title font size
  borderColor = "grey85",  # Light gray border
  frame = FALSE,  # Remove outer borders
  grid = TRUE  # Add gridlines
)
# Close the PNG device
dev.off()

# setLabels = c("pig", "human", "mouse", "marmoset", "cattle","rat")
set <- 6
species <- "rat"
expr_normalized <- multiExpr[[set]]$data
# create plot of module-trait correlations
mergedColors = labels2colors(netwk_rat$colors)
moduleColors=mergedColors


if(T){
  nGenes = ncol(expr_normalized)  # Number of genes
  nSamples = nrow(expr_normalized)  # Number of samples
  
  # Calculate module eigengenes
  MES0 <- moduleEigengenes(expr_normalized, moduleColors)$eigengenes
  MEs = orderMEs(MES0)  # Sort module eigengenes
  
  
  # Compute module-trait correlations and p-values
  moduleTraitCor <- cor(MEs, Traits[[set]]$data[, colSums(Traits[[set]]$data != 0) > 0], use = "p")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
  
  # Create a formatted text matrix for displaying correlations and p-values
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
  dim(textMatrix) = dim(moduleTraitCor)  # Adjust dimensions of the text matrix
  
  # Set the width and height dynamically based on the size of the heatmap
  heatmap_width <- 2 + 3.5 * ncol(moduleTraitCor)  # Customize as needed based on number of traits
  heatmap_height <- 2 + 2.5 * nrow(moduleTraitCor)  # Customize based on number of modules
  
  # Save the plot to a PDF file
  pdf("V2_module-trait-relationships_rat_adjusted_font.pdf", 
      width = heatmap_width, height = heatmap_height)
  
  # Adjust margins for a cleaner layout
  par(mar = c(12, 12, 8, 6))
  
  # Use medium red and blue contrast
  labeledHeatmap(
    Matrix = moduleTraitCor, 
    xLabels = names(Traits[[set]]$data[, colSums(Traits[[set]]$data != 0) > 0]),
    yLabels = names(MEs),
    ySymbols = names(MEs),
    colorLabels = FALSE,
    colors = colorRampPalette(c("#87CEEB", "#FFFFFF", "#FF4C4C"))(50),  # Medium blue to soft red contrast
    textMatrix = textMatrix,
    setStdMargins = FALSE,
    cex.text = 1.5,  # Adjust text size for the correlation values
    cex.lab = 1.5,  # Larger font size for module labels (you can increase this further if needed)
    zlim = c(-1, 1),  # Correlation range from -1 to 1
    main = paste("Module-Trait Relationships in", species),
    cex.main = 1.5,  # Larger title font size
    borderColor = "grey85",  # Light gray border for separation
    frame = FALSE,  # Remove heavy outer borders for a sleek look
    grid = TRUE  # Add gridlines for clarity
  )
  
  # Close the file
  dev.off()
}

heatmap_width <- 2 + 3.5 * ncol(moduleTraitCor)  # Adjust based on number of traits
heatmap_height <- 2 + 2.5 * nrow(moduleTraitCor)  # Adjust based on number of modules
# Set resolution and aspect ratio with png()
png("V2_module-trait-relationships_rat_adjusted_font.png",
    width = heatmap_width, height = heatmap_height, units = "in", res = 300)  # Set high resolution

# Adjust margins for a cleaner layout
par(mar = c(12, 12, 8, 6))

# Use medium red and blue contrast in labeledHeatmap
labeledHeatmap(
  Matrix = moduleTraitCor, 
  xLabels = names(Traits[[set]]$data[, colSums(Traits[[set]]$data != 0) > 0]),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = colorRampPalette(c("#87CEEB", "#FFFFFF", "#FF4C4C"))(50),  # Medium blue to soft red contrast
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 1.5,  # Adjust text size for the correlation values
  cex.lab = 1.5,  # Larger font size for module labels
  zlim = c(-1, 1),  # Correlation range from -1 to 1
  main = paste("Module-Trait Relationships in", species),
  cex.main = 1.5,  # Larger title font size
  borderColor = "grey85",  # Light gray border
  frame = FALSE,  # Remove outer borders
  grid = TRUE  # Add gridlines
)

# Close the PNG device
dev.off()

```

```{r}
# summarize sonserved human modules across all the species
ref = 1
test = 6 #(1,3,4,5,6)
species <- setLabels[test]

# 提取统计数据
statsObs <- cbind(mp$quality$observed[[ref]][[test]][, -1], 
                  mp$preservation$observed[[ref]][[test]][, -1])
statsZ <- cbind(mp$quality$Z[[ref]][[test]], 
                mp$preservation$Z[[ref]][[test]][, -1])


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

input1_filtered <- input1[input1$Zsummary.pres>2,]
input3_filtered <- input1[input1$Zsummary.pres>2,]
input4_filtered <- input1[input1$Zsummary.pres>2,]
input5_filtered <- input1[input1$Zsummary.pres>2,]
input6_filtered <- input1[input1$Zsummary.pres>2,]

all_intersections <- Reduce(intersect, list(
  input1_filtered$mergedColors_human,
  input3_filtered$mergedColors_human,
  input4_filtered$mergedColors_human,
  input5_filtered$mergedColors_human,
  input6_filtered$mergedColors_human
))

all_intersections