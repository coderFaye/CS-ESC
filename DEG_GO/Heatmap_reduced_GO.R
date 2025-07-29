setwd("/Users/yf358/Desktop/PR1/round4/Figure2")
load("n_p_deg.rds")
# merging go terms (down)
human_go <- data.frame(ID = go_down_human$term_id, 
                       Description = go_down_human$term_name,
                       human_pajdut = go_down_human$p_value)
go_term_all <- human_go

mouse_go <- data.frame(ID = go_down_mouse$term_id, mouse_pajdut = go_down_mouse$p_value)
go_term_all <- merge(go_term_all, mouse_go, by = "ID", all.x= T )

pig_go <- data.frame(ID = go_down_pig$term_id, pig_pajdut = go_down_pig$p_value)
go_term_all <- merge(go_term_all, pig_go, by = "ID", all.x= T )

marmoset_go <- data.frame(ID = go_down_marmoset$term_id, marmoset_pajdut = go_down_marmoset$p_value)
go_term_all <- merge(go_term_all, marmoset_go, by = "ID", all.x= T )

crab_macaque_go <- data.frame(ID = go_down_crab_macaque$term_id, crab_macaque_pajdut = go_down_crab_macaque$p_value)
go_term_all <- merge(go_term_all, crab_macaque_go, by = "ID", all.x= T )

go_term_all[is.na(go_term_all)] <- 0.99

selected_rows <- rowSums(go_term_all[, 3:7] < 0.05) >= 1
filtered_go_terms_down <- go_term_all[selected_rows, ]
write.csv(filtered_go_terms_down, "n_p_go_down.csv")

# merging go terms (up)
human_go <- data.frame(ID = go_up_human$term_id, 
                       Description = go_up_human$term_name,
                       human_pajdut = go_up_human$p_value)
go_term_all <- human_go

mouse_go <- data.frame(ID = go_up_mouse$term_id, mouse_pajdut = go_up_mouse$p_value)
go_term_all <- merge(go_term_all, mouse_go, by = "ID", all.x= T )

pig_go <- data.frame(ID = go_up_pig$term_id, pig_pajdut = go_up_pig$p_value)
go_term_all <- merge(go_term_all, pig_go, by = "ID", all.x= T )

marmoset_go <- data.frame(ID = go_up_marmoset$term_id, marmoset_pajdut = go_up_marmoset$p_value)
go_term_all <- merge(go_term_all, marmoset_go, by = "ID", all.x= T )

crab_macaque_go <- data.frame(ID = go_up_crab_macaque$term_id, crab_macaque_pajdut = go_up_crab_macaque$p_value)
go_term_all <- merge(go_term_all, crab_macaque_go, by = "ID", all.x= T )

go_term_all[is.na(go_term_all)] <- 0.99

selected_rows <- rowSums(go_term_all[, 3:7] < 0.05) >= 1
filtered_go_terms_up <- go_term_all[selected_rows, ]
write.csv(filtered_go_terms_up, "n_p_go_up.csv")

# Heatmaps
# create heatmap showing go terms
# Up
up_select <- read.table("n_p_up_go_select.txt")
colnames(up_select) <- "ID"
filtered_go_terms_up1 <- merge(up_select, filtered_go_terms_up, by.x = "ID")

if (all(filtered_go_terms_up1[, 3:7] > 0)) {
  filtered_go_terms_up1[, 3:7] <- -log(filtered_go_terms_up1[, 3:7])
} else {
  stop("All values must be positive to take the logarithm.")
}
filtered_go_terms_up1$Description
filtered_go_terms_up1 <- filtered_go_terms_up1[-c(7,8,9,16,18,19,20,21,47), ]
# draw heatmap
heatmap_df <- as.matrix(filtered_go_terms_up1[, 3:7])
colnames(heatmap_df) <- c("human","mouse","pig","marmoset","C.macaque")
col_order <- c("C.macaque", "human", "marmoset", "mouse", "pig")
heatmap_df <- heatmap_df[, col_order]
rownames(heatmap_df) <- filtered_go_terms_up1$Description
rownames(heatmap_df) <- gsub("signaling pathway", "sig. path.", rownames(heatmap_df))
rownames(heatmap_df) <- gsub("regulation of", "reg. of", rownames(heatmap_df))

# Define a new color gradient: from dark blue (low), through white (neutral), to dark red (high)
breaks <- c(seq(0, 0.99, length.out = 50), seq(1, 10, length.out = 50))
colors <- colorRampPalette(c("#313695", "white", "#A50026"))(length(breaks))  # Smooth blue-white-red gradient

# Draw the heatmap with refined aesthetics
p <- pheatmap(heatmap_df, 
              color = colors, 
              breaks = breaks, 
              labels_row = rownames(heatmap_df), 
              name = "-log10(P)",  # Keep a meaningful label for the legend
              treeheight_col = 20,  # Adjust tree height for a balanced look
              treeheight_row = 40,  # Make the row tree slightly larger for better differentiation
              fontsize = 14,  # Increase font size for better readability in high-res figures
              legend = TRUE,  # Keep the legend for understanding the scale
              border_color = NA,  # Remove borders for a clean, modern look
              cluster_cols = F,  # Cluster the columns
              cluster_rows = T,  # Cluster the rows
              cellwidth = 20,  # Adjust cell width for better readability in publications
              cellheight = 20)  # Adjust cell height
p

png("go_heatmap_n_p_up_all_species_v4.png",res = 300, width=2800, height=4600)
print(p)
dev.off()

svg("go_heatmap_n_p_up_all_species_v4.svg", width = 9.33, height = 15.33)  
print(p)
dev.off()

# Down
down_select <- read.table("n_p_down_go_select.txt")
colnames(down_select) <- "ID"
filtered_go_terms_down1 <- merge(down_select, filtered_go_terms_down, by.x = "ID")

if (all(filtered_go_terms_down1[, 3:7] > 0)) {
  filtered_go_terms_down1[, 3:7] <- -log(filtered_go_terms_down1[, 3:7])
} else {
  stop("All values must be positive to take the logarithm.")
}

filtered_go_terms_down1$Description
# draw heatmap
heatmap_df <- as.matrix(filtered_go_terms_down1[, 3:7])
colnames(heatmap_df) <- c("human","mouse","pig","marmoset","C.macaque")
col_order <- c("C.macaque", "human", "marmoset", "mouse", "pig")
heatmap_df <- heatmap_df[, col_order]
rownames(heatmap_df) <- filtered_go_terms_down1$Description

rownames(heatmap_df) <- gsub("signaling pathway", "sig. path.", rownames(heatmap_df))
rownames(heatmap_df) <- gsub("regulation of", "reg. of", rownames(heatmap_df))

# Define a new color gradient: from dark blue (low), through white (neutral), to dark red (high)
breaks <- c(seq(0, 0.99, length.out = 50), seq(1, 10, length.out = 50))
colors <- colorRampPalette(c("#313695", "white", "#A50026"))(length(breaks))  # Smooth blue-white-red gradient

# Draw the heatmap with refined aesthetics
p <- pheatmap(heatmap_df, 
              color = colors, 
              breaks = breaks, 
              labels_row = rownames(heatmap_df), 
              name = "-log10(P)",  # Keep a meaningful label for the legend
              treeheight_col = 20,  # Adjust tree height for a balanced look
              treeheight_row = 40,  # Make the row tree slightly larger for better differentiation
              fontsize = 14,  # Increase font size for better readability in high-res figures
              legend = TRUE,  # Keep the legend for understanding the scale
              border_color = NA,  # Remove borders for a clean, modern look
              cluster_cols = F,  # Cluster the columns
              cluster_rows = T,  # Cluster the rows
              cellwidth = 20,  # Adjust cell width for better readability in publications
              cellheight = 20)  # Adjust cell height
p

png("go_heatmap_n_p_down_all_species_v4.png",res = 300, width=2800, height=3500)
print(p)
dev.off()

svg("go_heatmap_n_p_down_all_species_v4.svg", width = 12.33, height = 15.33)
print(p)
dev.off()


# Reduce go terms
# BiocManager::install("rrvgo")
# BiocManager::install("org.Hs.eg.db")

library(rrvgo)
filtered_go_terms_down_filtered <- filtered_go_terms_down[
  rowSums(filtered_go_terms_down[, 3:7] < 0.05) >= 4, 
]

filtered_go_terms_up_filtered <- filtered_go_terms_up[
  rowSums(filtered_go_terms_up[, 3:7] < 0.05) >= 4, 
]

# List of GO terms you want to visualize
go_terms <- filtered_go_terms_up_filtered$ID 
simMatrix <- calculateSimMatrix(go_terms,
                                orgdb="org.Hs.eg.db",
                                ont="BP",
                                method="Rel")
#scores <- setNames(-log10(go_analysis$qvalue), go_analysis$ID)
reducedTerms <- reduceSimMatrix(simMatrix,
                                threshold=0.8,
                                orgdb="org.Hs.eg.db")
plot <- scatterPlot(simMatrix, reducedTerms,labelSize = 5)
plot

# Now, apply your custom color scale and other customizations
custom_colors <- c("#313695", "#A3C1DA", "#264653", "#A50026", 
                   "#FF6F61", "#1f77b4", "#ff7f0e", "#2ca02c", 
                   "#d62728", "#9467bd", "#8c564b", "#e377c2", 
                   "#7f7f7f", "#bcbd22", "#17becf", "#fdae61")
expanded_colors <- colorRampPalette(custom_colors)(37)

# Apply the custom color palette
plot <- plot +
  theme_minimal() +  # Start with a minimal theme
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_rect(fill = "white", color = NA),  # Set white background
    plot.background = element_rect(fill = "white", color = NA),  # Set white plot background
    axis.line = element_blank(),  # Keep black axis lines (frame)
    axis.text = element_blank(),  # Remove axis text
    axis.ticks = element_blank(),  # Remove axis ticks
    axis.title = element_blank(),  # Remove axis titles
    legend.position = "right"  # Keep the legend on the right
  ) +
  scale_color_manual(values = expanded_colors)


png("simigo_up_v3.png",res = 300, width=5900, height=2500)
print(plot)
dev.off()

# Save the plot as an SVG file
svg("simigo_up_v3.svg", width = 25, height = 16)  # Adjust width and height as needed
print(plot)
dev.off()




# List of GO terms you want to visualize
go_terms <- filtered_go_terms_down_filtered$ID 
simMatrix <- calculateSimMatrix(go_terms,
                                orgdb="org.Hs.eg.db",
                                ont="BP",
                                method="Rel")
#scores <- setNames(-log10(go_analysis$qvalue), go_analysis$ID)
reducedTerms <- reduceSimMatrix(simMatrix,
                                threshold=0.8,
                                orgdb="org.Hs.eg.db")
plot <- scatterPlot(simMatrix, reducedTerms, labelSize = 5)
# Remove the existing color scale to avoid conflicts
plot <- plot + scale_color_discrete(NULL)  # This removes the existing color scale

# Now, apply your custom color scale and other customizations
custom_colors <- c("#313695", "#A3C1DA", "#264653", "#A50026", 
                   "#FF6F61", "#1f77b4", "#ff7f0e", "#2ca02c", 
                   "#d62728", "#9467bd", "#8c564b", "#e377c2", 
                   "#7f7f7f", "#bcbd22", "#17becf", "#fdae61")

# Apply the custom color palette
plot <- plot +
  theme_minimal() +  # Start with a minimal theme
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    panel.background = element_rect(fill = "white", color = NA),  # Set white background
    plot.background = element_rect(fill = "white", color = NA),  # Set white plot background
    axis.line = element_blank(),  # Keep black axis lines (frame)
    axis.text = element_blank(),  # Remove axis text
    axis.ticks = element_blank(),  # Remove axis ticks
    axis.title = element_blank(),  # Remove axis titles
    legend.position = "right"  # Keep the legend on the right
  ) +
  scale_color_manual(values = custom_colors)

png("simigo_down_v3.png",res = 300, width=3900, height=1400)
print(plot)
dev.off()

# Save the plot as an SVG file
svg("simigo_down_v3.svg", width = 16, height = 5)  # Adjust width and height as needed
print(plot)
dev.off()