setwd("/Users/yf358/Desktop/PR1/Final_version/transcriptome_amalgamation/")

library(ComplexHeatmap)
library(circlize)

# Load sample metadata and standardize pluripotency state labels
sampleinfo <- read.table("v1_sampleinfo.txt", header = TRUE)
sampleinfo$type[grepl("ex", sampleinfo$type)] <- "expanded/extended"
sampleinfo$type[grepl("na", sampleinfo$type)] <- "naïve"

# Prepare ordering of samples by species and pluripotency state
ordering_df <- data.frame(
  species = sampleinfo$species,
  type = sampleinfo$type,
  original_order = seq_along(sampleinfo$species)
)
ordering_df <- ordering_df[order(ordering_df$species, ordering_df$type), ]

# Color palettes for species and pluripotency states
cor_species <- c(
  "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", 
  "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22"
)
cor_type <- c(
  "#17becf", "#ff9896", "#c5b0d5", "#9edae5"
)

# Heatmap top annotation
top <- HeatmapAnnotation(
  border = FALSE,
  Pluripotency = ordering_df$type,
  Species = ordering_df$species,
  col = list(
    Pluripotency = c(
      "formative" = cor_type[1],
      "expanded/extended" = cor_type[2],
      "primed" = cor_type[3],
      "naïve" = cor_type[4]
    ),
    Species = c(
      "cattle" = cor_species[1], "human" = cor_species[2],
      "mouse" = cor_species[3], "pig" = cor_species[4],
      "sheep" = cor_species[5], "horse" = cor_species[6],
      "rat" = cor_species[7], "crab_macaque" = cor_species[8],
      "marmoset" = cor_species[9]
    )
  )
)

# Define heatmap color scale
myColor <- colorRamp2(c(0.6, 0.8, 1), c("#53A5FE", "white", "#FF6363"))

### --- Raw count correlation analysis --- ###
raw <- read.table("v1_rawcount.txt", header = TRUE)
raw <- raw[-c(1:17)]  # Remove metadata columns
raw <- raw[, ordering_df$original_order]
cor_raw <- round(cor(raw), 4)

# Plot and export heatmap
heatmap_raw <- Heatmap(
  cor_raw,
  name = "Correlation",
  col = myColor,
  top_annotation = top,
  show_row_names = FALSE,
  show_column_names = FALSE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(
    title = "Correlation",
    legend_direction = "horizontal",
    title_gp = gpar(fontsize = 12, fontface = "bold"),
    labels_gp = gpar(fontsize = 10)
  ),
  border = FALSE
)
png("raw_cor.png", res = 300, width = 3000, height = 2000)
draw(heatmap_raw)
dev.off()

### --- Log-TPM correlation analysis --- ###
tpm <- read.table("v1_tpm.txt", header = TRUE)
tpm <- log2(tpm[, -c(1:17)] + 1)  # Log transform TPM after removing metadata
tpm <- tpm[, ordering_df$original_order]
cor_tpm <- round(cor(tpm), 4)

heatmap_tpm <- Heatmap(
  cor_tpm,
  name = "Correlation",
  col = myColor,
  top_annotation = top,
  show_row_names = FALSE,
  show_column_names = FALSE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(
    title = "Correlation",
    legend_direction = "horizontal",
    title_gp = gpar(fontsize = 12, fontface = "bold"),
    labels_gp = gpar(fontsize = 10)
  ),
  border = FALSE
)
png("logtpm_cor.png", res = 300, width = 3000, height = 2000)
draw(heatmap_tpm)
dev.off()

### --- SVA-corrected log-TPM correlation analysis --- ###
tc <- read.table("v1_sva_log_tpm.txt", header = TRUE)
tc <- tc[, -c(1:17)]
tc <- tc[, ordering_df$original_order]
cor_tc <- round(cor(tc), 4)

heatmap_tc <- Heatmap(
  cor_tc,
  name = "Correlation",
  col = myColor,
  top_annotation = top,
  show_row_names = FALSE,
  show_column_names = FALSE,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_heatmap_legend = TRUE,
  heatmap_legend_param = list(
    title = "Correlation",
    legend_direction = "horizontal",
    title_gp = gpar(fontsize = 12, fontface = "bold"),
    labels_gp = gpar(fontsize = 10)
  ),
  border = FALSE
)
png("sva_tpm_heatmap_v2.png", res = 300, width = 3000, height = 2000)
draw(heatmap_tc)
dev.off()

# Save final R workspace (alternative: save each correlation matrix with saveRDS)
save.image("correlation.rds")