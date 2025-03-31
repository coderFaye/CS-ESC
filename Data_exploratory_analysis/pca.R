setwd("/Users/yf358/Desktop/PR1/Final_version/Data_exploratory_analysis")

# Load required libraries
library(ggplot2)
library(grid)  # For unit()
library(dplyr)

# --- Load & preprocess data --- #

# Load log-transformed, SVA-corrected TPM data and remove metadata columns
tpm <- read.table("/Users/yf358/Desktop/PR1/Final_version/transcriptome_amalgamation/v1_sva_log_tpm.txt", header = TRUE)
tpm <- tpm[, -c(1:17)]

# Load sample metadata and unify type labels
sampleinfo <- read.table("/Users/yf358/Desktop/PR1/Final_version/transcriptome_amalgamation/v1_sampleinfo.txt", header = TRUE)
sampleinfo$type[grepl("ex", sampleinfo$type)] <- "expanded/extended"

# Transpose expression data to samples × genes and remove zero-variance genes
data <- t(tpm)
data_filtered <- data[, apply(data, 2, var) != 0]

# --- PCA computation --- #

data.pca <- prcomp(data_filtered, scale. = TRUE)
pca_result <- predict(data.pca)

# Extract top 3 PCs and calculate variance explained
pca_df <- as.data.frame(pca_result[, 1:3])
variance_explained <- summary(data.pca)$importance[2, 1:2]

# Prepare labels for axes
xlab <- paste0("PC1: ", round(variance_explained[1] * 100, 2), "% variance")
ylab <- paste0("PC2: ", round(variance_explained[2] * 100, 2), "% variance")

# Combine PCA result with sample metadata
PCA <- data.frame(
  PC1 = pca_result[, 1],
  PC2 = pca_result[, 2],
  PC3 = pca_result[, 3],
  Species = sampleinfo$species,
  Pluripotency = sampleinfo$type
)

# Define plotting levels and color palettes
pluripotency_levels <- c("expanded/extended", "naïve", "formative", "primed")
species_colors <- c(
  "cattle" = "#A6CEE3",     # Soft Blue
  "human" = "#1F78B4",      # Muted Dark Blue
  "mouse" = "#B2DF8A",      # Light Green
  "pig" = "#33A02C",        # Dark Green
  "sheep" = "#FB9A99",      # Soft Red
  "horse" = "#E31A1C",      # Muted Red
  "rat" = "#FDBF6F",        # Light Orange
  "crab_macaque" = "#FF7F00",  # Bold Orange
  "marmoset" = "#CAB2D6"    # Lavender
)

# --- PCA plot --- #

p <- ggplot(PCA, aes(x = PC1, y = PC2)) +
  geom_point(aes(fill = Species, shape = Pluripotency), size = 3.5, alpha = 1) +
  scale_fill_manual(values = species_colors) +
  scale_shape_manual(values = c(21, 22, 23, 24), limits = pluripotency_levels) +
  labs(x = xlab, y = ylab, title = "PCA Plot") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10),
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  ) +
  guides(
    fill = guide_legend(title = "Species", override.aes = list(color = species_colors)),
    shape = guide_legend(title = "Pluripotency")
  )

# View plot
p

# Save plot in high-quality formats
ggsave("pca.png", plot = p, width = 7, height = 5, dpi = 300)
ggsave("pca.svg", plot = p, width = 7, height = 5, dpi = 300)