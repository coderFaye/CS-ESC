setwd("/Users/yf358/Desktop/PR1/round3/05_deg")
#devtools::install_github("junjunlab/ClusterGVis")
library(DESeq2)
library(ClusterGVis)
library(dplyr)
library(IHW)

# deseq2-vst normalization
species <- "human"
sampleinfo_human <- read.table("./04.16_v1/v1_sampleinfo_human.txt", header = T)
sampleinfo_human$type[grepl("ex", sampleinfo_human$type)] <- "expanded/extended"
condition_human <- factor(paste0(sampleinfo_human$species, "_", sampleinfo_human$type))
condition_human

raw_human <- read.table("./04.16_v1/v1_raw_human.txt", header = T)
integer_counts_human <- apply(raw_human, c(1, 2), as.integer)
colData_human <- data.frame(row.names = colnames(integer_counts_human), condition_human)
dds_human <- DESeqDataSetFromMatrix(integer_counts_human, colData_human, 
                                    design = ~condition_human)
human_vst <- vst(dds_human)
human_mat <- as.data.frame(assay(human_vst))

# calculate the median value
gene_medians <- apply(human_mat, 1, function(gene_expr) {
  tapply(gene_expr, condition_human, median)
})
human_medians <- t(gene_medians)

# load the human degs
load("./mfuzz/v1/human_e_n.RData")
load("./mfuzz/v1/human_n_f.RData")
load("./mfuzz/v1/human_f_p.RData")
load("./mfuzz/v1/human_e_f.RData")
load("./mfuzz/v1/human_e_p.RData")
load("./mfuzz/v1/human_n_p.RData")

human_e_unique <- Reduce(intersect, list(
  human_e_n[human_e_n$log2FoldChange < 0, "ensembl_gene_id"], 
  human_e_f[human_e_f$log2FoldChange < 0, "ensembl_gene_id"], 
  human_e_p[human_e_p$log2FoldChange < 0, "ensembl_gene_id"]
))

human_n_unique <- Reduce(intersect, list(
  human_e_n[human_e_n$log2FoldChange > 0, "ensembl_gene_id"], 
  human_n_f[human_n_f$log2FoldChange < 0, "ensembl_gene_id"], 
  human_n_p[human_n_p$log2FoldChange < 0, "ensembl_gene_id"]
))

human_f_unique <- Reduce(intersect, list(
  human_e_f[human_e_f$log2FoldChange > 0, "ensembl_gene_id"], 
  human_n_f[human_n_f$log2FoldChange > 0, "ensembl_gene_id"], 
  human_f_p[human_f_p$log2FoldChange < 0, "ensembl_gene_id"]
))

human_p_unique <- Reduce(intersect, list(
  human_e_p[human_e_p$log2FoldChange > 0, "ensembl_gene_id"], 
  human_n_p[human_n_p$log2FoldChange > 0, "ensembl_gene_id"], 
  human_f_p[human_f_p$log2FoldChange > 0, "ensembl_gene_id"]
))          

human_deg <- unique(c(human_e_unique, human_n_unique, human_f_unique, human_p_unique))

# filter using 1-1 orthology
orthology <- read.csv("/Users/yf358/Desktop/PR1/round3/02_merge/orthologous_filter-9species.csv")
human_ortho <- orthology$Gene.stable.ID
human_deg_ortho <- intersect(human_deg, human_ortho)
human_deg_medians <- human_medians[human_deg_ortho,]
colnames(human_deg_medians)
human_deg_medians <- human_deg_medians[,c(1,3,2,4)]

# clustergvis
pdf('./mfuzz/v1/human_elbow.pdf',height = 6,width = 6)
getClusters(exp = human_deg_medians)
dev.off()
cm_human <- clusterData(exp = human_deg_medians,
                        cluster.method = "mfuzz",
                        cluster.num = 5)

pdf('./mfuzz/v1/human_mfuzz_10_6.pdf',height = 8,width = 15)
visCluster(object = cm_human,
            plot.type = "line")
dev.off()

pdf('./mfuzz/v1/human_cluster_10_6.pdf',height = 10,width = 8)
visCluster(object = cm_human,
           plot.type = "both",                   # Show both line and heatmap
           add.box = TRUE,                       # Add boxplot for cluster annotations
           boxcol = ggsci::pal_lancet()(4),      # Use a professional color palette for the boxplot
           ht.col.list = list(col_range = c(-2, 0, 2), col_color = c("#2C7BB6", "#FFFFFF", "#D7191C")), # Define a custom color range for heatmap
           border = TRUE,                        # Add borders to heatmap
           line.size = 0.5,                      # Adjust line size for line plots
           mline.size = 2,                       # Set median line size
           mline.col = "#CC3333",                # Set median line color to red
           add.mline = TRUE,                     # Add median lines to the plot
           show_row_names = FALSE,               # Don't show row names for cleaner heatmap
           cluster_columns = F,               # Cluster the columns for better visualization
           add.sampleanno = TRUE,                # Add sample annotation
           sample.col = ggsci::pal_lancet()(4),
           panel.arg = c(3,0.4,4,"grey90",NA),
           add.line = F,
           show_row_dend = F
)
dev.off()



######### mouse ##########
# deseq2-vst normalization
species <- "mouse"
sampleinfo_mouse <- read.table("./04.16_v1/v1_sampleinfo_mouse.txt", header = T)
sampleinfo_mouse$type[grepl("ex", sampleinfo_mouse$type)] <- "expanded/extended"
condition_mouse <- factor(paste0(sampleinfo_mouse$species, "_", sampleinfo_mouse$type))
condition_mouse

raw_mouse <- read.table("./04.16_v1/v1_raw_mouse.txt", header = T)
integer_counts_mouse <- apply(raw_mouse, c(1, 2), as.integer)
colData_mouse <- data.frame(row.names = colnames(integer_counts_mouse), condition_mouse)
dds_mouse <- DESeqDataSetFromMatrix(integer_counts_mouse, colData_mouse, 
                                    design = ~condition_mouse)
mouse_vst <- vst(dds_mouse)
mouse_mat <- as.data.frame(assay(mouse_vst))

# calculate the median value
gene_medians <- apply(mouse_mat, 1, function(gene_expr) {
  tapply(gene_expr, condition_mouse, median)
})
mouse_medians <- t(gene_medians)

# load the mouse degs
load("./mfuzz/v1/mouse_e_n.RData")
load("./mfuzz/v1/mouse_n_f.RData")
load("./mfuzz/v1/mouse_f_p.RData")
load("./mfuzz/v1/mouse_e_f.RData")
load("./mfuzz/v1/mouse_e_p.RData")
load("./mfuzz/v1/mouse_n_p.RData")

mouse_e_unique <- Reduce(intersect, list(
  mouse_e_n[mouse_e_n$log2FoldChange < 0, "ensembl_gene_id"], 
  mouse_e_f[mouse_e_f$log2FoldChange < 0, "ensembl_gene_id"], 
  mouse_e_p[mouse_e_p$log2FoldChange < 0, "ensembl_gene_id"]
))

mouse_n_unique <- Reduce(intersect, list(
  mouse_e_n[mouse_e_n$log2FoldChange > 0, "ensembl_gene_id"], 
  mouse_n_f[mouse_n_f$log2FoldChange < 0, "ensembl_gene_id"], 
  mouse_n_p[mouse_n_p$log2FoldChange < 0, "ensembl_gene_id"]
))

mouse_f_unique <- Reduce(intersect, list(
  mouse_e_f[mouse_e_f$log2FoldChange > 0, "ensembl_gene_id"], 
  mouse_n_f[mouse_n_f$log2FoldChange > 0, "ensembl_gene_id"], 
  mouse_f_p[mouse_f_p$log2FoldChange < 0, "ensembl_gene_id"]
))

mouse_p_unique <- Reduce(intersect, list(
  mouse_e_p[mouse_e_p$log2FoldChange > 0, "ensembl_gene_id"], 
  mouse_n_p[mouse_n_p$log2FoldChange > 0, "ensembl_gene_id"], 
  mouse_f_p[mouse_f_p$log2FoldChange > 0, "ensembl_gene_id"]
))          

mouse_deg <- unique(c(mouse_e_unique, mouse_n_unique, mouse_f_unique, mouse_p_unique))

# filter using 1-1 orthology
orthology <- read.csv("/Users/yf358/Desktop/PR1/round3/02_merge/orthologous_filter-9species.csv")
mouse_ortho <- orthology$Mouse.gene.stable.ID
mouse_deg_ortho <- intersect(mouse_deg, mouse_ortho)
mouse_deg_medians <- mouse_medians[mouse_deg_ortho,]

# clustergvis
pdf('./mfuzz/v1/mouse_elbow.pdf',height = 6,width = 6)
getClusters(exp = mouse_deg_medians)
dev.off()

colnames(mouse_deg_medians)
mouse_deg_medians <- mouse_deg_medians[,c(1,3,2,4)]
cm_mouse <- clusterData(exp = mouse_deg_medians,
                        cluster.method = "mfuzz",
                        cluster.num = 7)

pdf('./mfuzz/v1/mouse_mfuzz_10_6.pdf',height = 8,width = 14)
visCluster(object = cm_mouse,
            plot.type = "line")
dev.off()

pdf('./mfuzz/v1/mouse_cluster_10_6.pdf',height = 10,width = 8)
visCluster(object = cm_mouse,
           plot.type = "both",                   # Show both line and heatmap
           add.box = TRUE,                       # Add boxplot for cluster annotations
           boxcol = ggsci::pal_lancet()(4),      # Use a professional color palette for the boxplot
           ht.col.list = list(col_range = c(-2, 0, 2), col_color = c("#2C7BB6", "#FFFFFF", "#D7191C")), # Define a custom color range for heatmap
           border = TRUE,                        # Add borders to heatmap
           line.size = 0.5,                      # Adjust line size for line plots
           mline.size = 2,                       # Set median line size
           mline.col = "#CC3333",                # Set median line color to red
           add.mline = TRUE,                     # Add median lines to the plot
           show_row_names = FALSE,               # Don't show row names for cleaner heatmap
           cluster_columns = F,               # Cluster the columns for better visualization
           add.sampleanno = TRUE,                # Add sample annotation
           sample.col = ggsci::pal_lancet()(4),
           panel.arg = c(3,0.4,4,"grey90",NA),
           add.line = F,
           show_row_dend = F
)
dev.off()




######### pig ##########
# deseq2-vst normalization
species <- "pig"
sampleinfo_pig <- read.table("./04.16_v1/v1_sampleinfo_pig.txt", header = T)
sampleinfo_pig$type[grepl("ex", sampleinfo_pig$type)] <- "expanded/extended"
condition_pig <- factor(paste0(sampleinfo_pig$species, "_", sampleinfo_pig$type))
condition_pig

raw_pig <- read.table("./04.16_v1/v1_raw_pig.txt", header = T)
integer_counts_pig <- apply(raw_pig, c(1, 2), as.integer)
colData_pig <- data.frame(row.names = colnames(integer_counts_pig), condition_pig)
dds_pig <- DESeqDataSetFromMatrix(integer_counts_pig, colData_pig, 
                                    design = ~condition_pig)
pig_vst <- vst(dds_pig)
pig_mat <- as.data.frame(assay(pig_vst))

# calculate the median value
gene_medians <- apply(pig_mat, 1, function(gene_expr) {
  tapply(gene_expr, condition_pig, median)
})
pig_medians <- t(gene_medians)

# load the pig degs
load("./mfuzz/v1/pig_e_n.RData")
load("./mfuzz/v1/pig_n_f.RData")
load("./mfuzz/v1/pig_f_p.RData")
load("./mfuzz/v1/pig_e_f.RData")
load("./mfuzz/v1/pig_e_p.RData")
load("./mfuzz/v1/pig_n_p.RData")

pig_e_unique <- Reduce(intersect, list(
  pig_e_n[pig_e_n$log2FoldChange < 0, "ensembl_gene_id"], 
  pig_e_f[pig_e_f$log2FoldChange < 0, "ensembl_gene_id"], 
  pig_e_p[pig_e_p$log2FoldChange < 0, "ensembl_gene_id"]
))

pig_n_unique <- Reduce(intersect, list(
  pig_e_n[pig_e_n$log2FoldChange > 0, "ensembl_gene_id"], 
  pig_n_f[pig_n_f$log2FoldChange < 0, "ensembl_gene_id"], 
  pig_n_p[pig_n_p$log2FoldChange < 0, "ensembl_gene_id"]
))

pig_f_unique <- Reduce(intersect, list(
  pig_e_f[pig_e_f$log2FoldChange > 0, "ensembl_gene_id"], 
  pig_n_f[pig_n_f$log2FoldChange > 0, "ensembl_gene_id"], 
  pig_f_p[pig_f_p$log2FoldChange < 0, "ensembl_gene_id"]
))

pig_p_unique <- Reduce(intersect, list(
  pig_e_p[pig_e_p$log2FoldChange > 0, "ensembl_gene_id"], 
  pig_n_p[pig_n_p$log2FoldChange > 0, "ensembl_gene_id"], 
  pig_f_p[pig_f_p$log2FoldChange > 0, "ensembl_gene_id"]
))          

pig_deg <- unique(c(pig_e_unique, pig_n_unique, pig_f_unique, pig_p_unique))

# filter using 1-1 orthology
orthology <- read.csv("/Users/yf358/Desktop/PR1/round3/02_merge/orthologous_filter-9species.csv")
pig_ortho <- orthology$Pig.gene.stable.ID
pig_deg_ortho <- intersect(pig_deg, pig_ortho)
pig_deg_medians <- pig_medians[pig_deg_ortho,]
colnames(pig_deg_medians)
pig_deg_medians <- pig_deg_medians[,c(1,3,2,4)]

# clustergvis
pdf('./mfuzz/v1/pig_elbow.pdf',height = 6,width = 6)
getClusters(exp = pig_deg_medians)
dev.off()
cm_pig <- clusterData(exp = pig_deg_medians,
                      cluster.method = "mfuzz",
                      cluster.num = 5)
pdf('./mfuzz/v1/pig_mfuzz_10_6.pdf',height = 8,width = 14)
visCluster(object = cm_pig,
           plot.type = "line")
dev.off()

pdf('./mfuzz/v1/pig_cluster_10_6.pdf',height = 8,width = 7)
visCluster(object = cm_pig,
           plot.type = "both",                   # Show both line and heatmap
           add.box = TRUE,                       # Add boxplot for cluster annotations
           boxcol = ggsci::pal_lancet()(4),      # Use a professional color palette for the boxplot
           ht.col.list = list(col_range = c(-2, 0, 2), col_color = c("#2C7BB6", "#FFFFFF", "#D7191C")), # Define a custom color range for heatmap
           border = TRUE,                        # Add borders to heatmap
           line.size = 0.5,                      # Adjust line size for line plots
           mline.size = 2,                       # Set median line size
           mline.col = "#CC3333",                # Set median line color to red
           add.mline = TRUE,                     # Add median lines to the plot
           show_row_names = FALSE,               # Don't show row names for cleaner heatmap
           cluster_columns = F,               # Cluster the columns for better visualization
           add.sampleanno = TRUE,                # Add sample annotation
           sample.col = ggsci::pal_lancet()(4),
           panel.arg = c(3,0.4,4,"grey90",NA),
           add.line = F,
           show_row_dend = F
)
dev.off()
#save.image("mfuzz_v1.rds")


human_cluster <- data.frame(
  cluster = c(
    rep("expanded/extended", length(human_e_unique)),
    rep("naive", length(human_n_unique)),
    rep("formative", length(human_f_unique)),
    rep("primed", length(human_p_unique))
  ),
  ensembl_gene_id = c(human_e_unique, human_n_unique, human_f_unique, human_p_unique)
)
human_cluster_ortho <- human_cluster[human_cluster$ensembl_gene_id %in% human_deg_ortho, ]
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
symbol  <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                 filters = 'ensembl_gene_id',
                 values = human_cluster_ortho$ensembl_gene_id,
                 mart = ensembl)
human_cluster_ortho <- merge(human_cluster_ortho, symbol, by = "ensembl_gene_id")
human_cluster_ortho <- human_cluster_ortho[ ,-1]

mouse_cluster <- data.frame(
  cluster = c(
    rep("expanded/extended", length(mouse_e_unique)),
    rep("naive", length(mouse_n_unique)),
    rep("formative", length(mouse_f_unique)),
    rep("primed", length(mouse_p_unique))
  ),
  ensembl_gene_id = c(mouse_e_unique, mouse_n_unique, mouse_f_unique, mouse_p_unique)
)
mouse_cluster_ortho <- mouse_cluster[mouse_cluster$ensembl_gene_id %in% mouse_deg_ortho, ]
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
symbol  <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                 filters = 'ensembl_gene_id',
                 values = mouse_cluster_ortho$ensembl_gene_id,
                 mart = ensembl)
mouse_cluster_ortho <- merge(mouse_cluster_ortho, symbol, by = "ensembl_gene_id")
mouse_cluster_ortho <- mouse_cluster_ortho[ ,-1]
mouse_cluster_ortho$external_gene_name <- toupper(mouse_cluster_ortho$external_gene_name)

pig_cluster <- data.frame(
  cluster = c(
    rep("expanded/extended", length(pig_e_unique)),
    rep("naive", length(pig_n_unique)),
    rep("formative", length(pig_f_unique)),
    rep("primed", length(pig_p_unique))
  ),
  ensembl_gene_id = c(pig_e_unique, pig_n_unique, pig_f_unique, pig_p_unique)
)
pig_cluster_ortho <- pig_cluster[pig_cluster$ensembl_gene_id %in% pig_deg_ortho, ]
ensembl <- useMart("ensembl", dataset = "sscrofa_gene_ensembl")
symbol  <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                 filters = 'ensembl_gene_id',
                 values = pig_cluster_ortho$ensembl_gene_id,
                 mart = ensembl)
pig_cluster_ortho <- merge(pig_cluster_ortho, symbol, by = "ensembl_gene_id")
pig_cluster_ortho <- pig_cluster_ortho[ ,-1]


#### River Plot ####
library(networkD3)
# Merge human and mouse clusters by gene
mouse_human <- full_join(mouse_cluster_ortho, human_cluster_ortho, by = "external_gene_name", suffix = c("_mouse", "_human"))
mouse_human <- as.data.frame(mouse_human)

# Count the number of overlapping genes between human and mouse clusters
mouse_human_flow <- mouse_human %>%
  count(cluster_mouse,cluster_human) %>%
  filter(!is.na(cluster_human) & !is.na(cluster_mouse))
# Add species prefix to the cluster names for human and mouse
mouse_human_flow <- mouse_human_flow %>%
  mutate(cluster_human = paste("human", cluster_human, sep = "_"),
         cluster_mouse = paste("mouse", cluster_mouse, sep = "_"))
# Create nodes from unique clusters, now with species prefixes
nodes <- data.frame(name = unique(c(mouse_human_flow$cluster_mouse, mouse_human_flow$cluster_human)))
# Create links (flows) between human and mouse clusters with counts
library(dplyr)
links <- mouse_human_flow %>%
  mutate(source = match(cluster_mouse, nodes$name) - 1,  # source cluster
         target = match(cluster_human, nodes$name) - 1,  # target cluster
         value = n) %>%
  dplyr::select(source, target, value)
# Check structure
head(links)
head(nodes)
links$group <- as.factor(human_mouse_flow$cluster_mouse)

# Create a color palette for the nodes
node_colors <- 'd3.scaleOrdinal()
                .domain(["human_formative", "human_naive", "human_expanded/extended", "human_primed",
                         "mouse_formative", "mouse_naive", "mouse_expanded/extended", "mouse_primed"])
                .range(["#1f78b4", "#33a02c", "#e31a1c", "#a6cee3",
                        "#ff7f00", "#6a3d9a", "#b15928", "#cab2d6"])' 

# Create the Sankey plot with custom node colors
sankey <- sankeyNetwork(Links = links, 
              Nodes = nodes, 
              Source = "source", 
              Target = "target", 
              Value = "value", 
              NodeID = "name",  # Node labels
              LinkGroup = "group",  # Group colors for links
              sinksRight = FALSE,  # Position the nodes
              nodeWidth = 40,  # Node width
              fontSize = 0,  # Node font size
              nodePadding = 4,  # Padding between nodes
              colourScale = node_colors)  # Apply the custom color scale
sankey
saveNetwork(sankey, file = './mfuzz/v1/mouse_human_flow.html')


# Merge human and pig clusters by gene
human_pig <- full_join(human_cluster_ortho, pig_cluster_ortho, by = "external_gene_name", suffix = c("_human", "_pig"))
human_pig <- as.data.frame(human_pig)

# Count the number of overlapping genes between human and pig clusters
human_pig_flow <- human_pig %>%
  count(cluster_human, cluster_pig) %>%
  filter(!is.na(cluster_human) & !is.na(cluster_pig))
# Add species prefix to the cluster names for human and pig
human_pig_flow <- human_pig_flow %>%
  mutate(cluster_human = paste("human", cluster_human, sep = "_"),
         cluster_pig = paste("pig", cluster_pig, sep = "_"))
# Create nodes from unique clusters, now with species prefixes
nodes <- data.frame(name = unique(c(human_pig_flow$cluster_human, human_pig_flow$cluster_pig)))
# Create links (flows) between human and pig clusters with counts
links <- human_pig_flow %>%
  mutate(source = match(cluster_human, nodes$name) - 1,  # source cluster
         target = match(cluster_pig, nodes$name) - 1,  # target cluster
         value = n) %>%
  select(source, target, value)
# Check structure
head(links)
head(nodes)
links$group <- as.factor(human_pig_flow$cluster_human)

# Create a color palette for the nodes
node_colors <- 'd3.scaleOrdinal()
                .domain(["human_formative", "human_naive", "human_expanded/extended", "human_primed",
                         "pig_formative", "pig_naive", "pig_expanded/extended", "pig_primed"])
                .range(["#1f78b4", "#33a02c", "#e31a1c", "#a6cee3",
                        "#ff7f00", "#6a3d9a", "#b15928", "#cab2d6"])' 

# Create the Sankey plot with custom node colors
sankey <- sankeyNetwork(Links = links, 
              Nodes = nodes, 
              Source = "source", 
              Target = "target", 
              Value = "value", 
              NodeID = "name",  # Node labels
              LinkGroup = "group",  # Group colors for links
              sinksRight = T,  # Position the nodes
              nodeWidth = 20,  # Node width
              fontSize = 0,  # Node font size
              nodePadding = 4,  # Padding between nodes
              colourScale = node_colors)  # Apply the custom color scale
saveNetwork(sankey, file = './mfuzz/v1/human_pig_flow.html')

