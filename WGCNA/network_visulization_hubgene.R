setwd("/Users/yf358/Desktop/PR1/round3/06_wgcna/v1")
load("wgcna.rdata")

# The human blue module
### Part1. draw the correlation plot for the genes within one module with the trait ###
set = 2
trait = Traits[[set]]$data$primed
names(trait) <- "primed"
expr_normalized <- multiExpr[[set]]$data
moduleColors = labels2colors(netwk_human$colors)

MES0 <- moduleEigengenes(expr_normalized,moduleColors)$eigengenes  # 计算模块特征向量
MEs = orderMEs(MES0)
modNames = substring(names(MEs), 3) # names (colors) of the modules

# calculate the correlation between the normalized gene expression data (expr_normalized) and the module eigengenes (MEs)
geneModuleMembership = as.data.frame(cor(expr_normalized, MEs, use = "p"))
nSamples = nrow(expr_normalized)
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")
# calculate gene-trait significance 
geneTraitSignificance = as.data.frame(cor(expr_normalized, trait, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

module = "blue" #########################putting the color below the plot
column = match(module, modNames)
moduleGenes = moduleColors==module

#sizeGrWindow(7, 7)
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for primed state",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)



### Part2. Exporting the network to a cytoscape format ###
#Recalculating topological overlap, if necessary
TOM = TOMsimilarityFromExpr(expr_normalized, power = 10);
module = "blue"; #chose modules that u want to export
#Select the gene modules
genes = colnames(expr_normalized)
inModule = (moduleColors==module)
modGenes = genes[inModule]

#Select the corresponding topologic overlap 
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modGenes, modGenes)

#Export the network in list files os n edges that cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = "CytoscapeEdgeFile.txt",
                               nodeFile = "CytoscapeNodeFile.txt",
                               weighted = TRUE,
                               threshold = 0.02,  # can adjust
                               nodeNames = modGenes,
                               nodeAttr = moduleColors[inModule])


### Part3. Identify top hub genes ###
library(igraph)
# Create an igraph object from the adjacency matrix (TOM)
modGraph <- graph.adjacency(modTOM, mode = "undirected", weighted = TRUE, diag = FALSE)
modGraph <- simplify(modGraph)
hub_score <- hub_score(modGraph, scale = TRUE, weights = NULL, options = arpack_defaults)
top_hub_genes <- names(sort(hub_score[["vector"]], decreasing = TRUE))[1:10]
top_hub_genes


### Part4.Visulize the network - interactive mode ###
library(networkD3)
# Convert the TOM matrix to a forceNetwork object
edges <- as.data.frame(as.table(modTOM))
edges <- edges[edges$Freq > 0.27, ]  # Apply a threshold #
colnames(edges) <- c("source", "target", "weight")

node_names <- unique(c(edges$source, edges$target))
nodes <- data.frame(name = node_names)
nodes$group <- ifelse(nodes$name %in% top_hub_genes, 1, 2)

# Convert source and target to zero-based indices
edges$source <- match(edges$source, node_names) - 1
edges$target <- match(edges$target, node_names) - 1

network <- forceNetwork(Links = edges,
                        Nodes = nodes,
                        Source = "source",
                        Target = "target",
                        Value = "weight",
                        NodeID = "name",
                        Group = "group",
                        opacity = 0.8,
                        zoom = T)

# save the network
saveNetwork(network, "gene_interaction_network.html")

# static visulization
library("geomnet")
edges <- as.data.frame(as.table(modTOM))
edges <- edges[edges$Freq > 0.20, ]  # Apply a threshold
colnames(edges) <- c("source", "target", "weight")
edges$source <- as.character(edges$source)
edges$target <- as.character(edges$target)
# Remove self-loops
edges <- edges[edges$source != edges$target, ]

node_names <- unique(c(edges$source, edges$target))
nodes <- data.frame(name = node_names)
hub_names <- names(hub_score[["vector"]])
nodes$hub_score <- hub_score[["vector"]][match(nodes$name, hub_names)]

netall <- fortify(as.edgedf(edges), nodes)
head(netall)

pgeomnet <- 
  ggplot(data = netall, aes(from_id = from_id, to_id = to_id)) +
  geom_net(layout.alg = 'fruchtermanreingold',
           aes(color = hub_score),  # 根据权重调整边的颜色
           linewidth = 0.8,      # 增加线条宽度
           size = 7,             # 增大节点大小
           alpha = 0.85,
           #labelon = T
  ) +       
  scale_color_viridis_c(option = "C", end = 0.5) + 
  theme_minimal(base_size = 15) +  # 使用简洁主题并调整基础字体大小
  theme(
    legend.position = "bottom",  # 将图例置于底部
    legend.title = element_text(size = 14, face = "bold"),  # 调整图例标题字体
    legend.text = element_text(size = 8),  # 调整图例字体
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # 居中标题
    axis.title = element_blank(),  # 去除轴标题
    axis.text = element_blank(),   # 去除轴标签
    panel.grid = element_blank(),  # 去除网格线
    plot.margin = margin(10, 10, 10, 10)  # 调整图表边距
  ) +
  labs(color = "Hub_score", title = "Gene Interaction Network")  # 添加图例和标题
pgeomnet

png("Gene_Interaction_Network_label.png", width = 15, height = 12, units = "in", res = 300)
print(pgeomnet)
dev.off()


### For the black module ###
#Recalculating topological overlap, if necessary
TOM = TOMsimilarityFromExpr(expr_normalized, power = 10);
module = c("black"); #chose modules that u want to export
#Select the gene modules
genes = colnames(expr_normalized)
inModule = (moduleColors==module)
modGenes = genes[inModule]

#Select the corresponding topologic overlap 
modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modGenes, modGenes)

#Export the network in list files os n edges that cytoscape can read
# cyt = exportNetworkToCytoscape(modTOM,
#                                edgeFile = "CytoscapeEdgeFile.txt",
#                                nodeFile = "CytoscapeNodeFile.txt",
#                                weighted = TRUE,
#                                threshold = 0.02,  # can adjust
#                                nodeNames = modGenes,
#                                nodeAttr = moduleColors[inModule])


### Part3. Identify top hub genes ###
library(igraph)
# Create an igraph object from the adjacency matrix (TOM)
modGraph <- graph.adjacency(modTOM, mode = "undirected", weighted = TRUE, diag = FALSE)
modGraph <- simplify(modGraph)
hub_score <- hub_score(modGraph, scale = TRUE, weights = NULL, options = arpack_defaults)
top_hub_genes <- names(sort(hub_score[["vector"]], decreasing = TRUE))[1:10]
top_hub_genes


### Part4.Visulize the network - interactive mode ###
library(networkD3)
# Convert the TOM matrix to a forceNetwork object
edges <- as.data.frame(as.table(modTOM))
edges <- edges[edges$Freq > 0.18, ]  # Apply a threshold #
colnames(edges) <- c("source", "target", "weight")

node_names <- unique(c(edges$source, edges$target))
nodes <- data.frame(name = node_names)
nodes$group <- ifelse(nodes$name %in% top_hub_genes, 1, 2)

# Convert source and target to zero-based indices
edges$source <- match(edges$source, node_names) - 1
edges$target <- match(edges$target, node_names) - 1

network <- forceNetwork(Links = edges,
                        Nodes = nodes,
                        Source = "source",
                        Target = "target",
                        Value = "weight",
                        NodeID = "name",
                        Group = "group",
                        opacity = 0.8,
                        zoom = T)

# # save the network
saveNetwork(network, "gene_interaction_network_black.html")