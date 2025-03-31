setwd("/Users/yf358/Desktop/PR1/round3/03_pca_batch/03_pca_batch_04.17")

# Load required libraries
library(ape)
library(ggtree)
library(ggalluvial)
library(ggprism)
library(reshape)
library(dplyr)  # Needed for group_by and summarise
library(ggplot2)

### --- Plot species tree --- ###
tree <- read.tree("merged.raxml.bestTree")

# Replace generic tip labels with actual species names
species_id <- c(
  "Rattus_norvegicus", "Mus_musculus", "Callithrix_jacchus", "Homo_sapiens", 
  "Macaca_fascicularis", "Ovis_aries", "Bos_taurus", "Sus_scrofa", "Equus_caballus"
)
tree$tip.label <- species_id

# Plot species tree using ggtree
p <- ggtree(tree) + 
  geom_tiplab(align = FALSE, size = 3)
p
ggsave("species_tree_v2.png", p, width = 6, height = 4, units = "in")
ggsave("species_tree_v2.svg", p, width = 6, height = 6, units = "in")

### --- Stacked barplot of sample numbers per species and state --- ###

# Load sample metadata and standardize type labels
sampleinfo <- read.table("v1_sampleinfo.txt", header = TRUE)
sampleinfo$type[grepl("ex", sampleinfo$type)] <- "expanded/extended"
sampleinfo$type[grepl("na", sampleinfo$type)] <- "naïve"

# Summarize sample count per species and type
sample_species <- sampleinfo %>% 
  group_by(species, type) %>% 
  summarise(count = n(), .groups = "drop")

# Convert to wide format and save for record-keeping
table <- xtabs(count ~ species + type, data = sample_species)
write.csv(table, "sample_group.csv", row.names = TRUE)

# Reload and reshape for plotting
sample_group <- read.csv("sample_group.csv")
colnames(sample_group)[1] <- "species"
df1 <- melt(sample_group, id.vars = "species")
names(df1)[2:3] <- c("type", "number")

# Factor species order for consistent x-axis arrangement
df1$species <- factor(df1$species, levels = c(
  "horse", "pig", "sheep", "cattle", "rat", "mouse",
  "marmoset", "human", "crab_macaque"
))

# Define custom color palette for pluripotency types
col <- c("#C2C1E0", "#264653", "#A3C1DA", "#696969")

# Generate the stacked bar plot
p1 <- ggplot(df1, aes(x = species, y = number, fill = type)) +
  geom_col(position = 'stack', width = 0.7) +
  scale_y_continuous(expand = c(0, 0), labels = scales::comma) +
  labs(
    x = "Species", 
    y = "Sample Numbers", 
    fill = "Pluripotency",
    title = "Sample Distribution Across Species"
  ) +
  guides(fill = guide_legend(
    keywidth = 1, keyheight = 1, 
    title.position = "top", title.hjust = 0.5
  )) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(
      color = "black", size = 12, family = "Helvetica",
      face = "plain", angle = 45, hjust = 1
    ),
    axis.text.y = element_text(color = "black", size = 12, family = "Helvetica"),
    axis.title = element_text(face = "bold", size = 14),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, face = "italic", size = 12),
    legend.position = "right",
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12)
  ) +
  scale_fill_manual(values = col)
p1

# Save plot
ggsave("stack_plot_col2_v2.png", p1, width = 7, height = 7)
ggsave("stack_plot_col2_v2.svg", p1, width = 7, height = 7)