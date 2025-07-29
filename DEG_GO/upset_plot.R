setwd("/Users/yf358/Desktop/PR1/round4/Figure2")
# Upset plots
#upset figure: up
gene_list_human <- up_human[up_human != ""]
gene_list_mouse <- up_mouse[up_mouse != ""]
gene_list_pig <- up_pig[up_pig != ""]
gene_list_crab_macaque <- up_crab_macaque[up_crab_macaque != ""]
gene_list_marmoset <- up_marmoset[up_marmoset != ""]
gene_df <- list(human = gene_list_human, mouse = toupper(gene_list_mouse), 
                pig = gene_list_pig, crab_macaque = gene_list_crab_macaque, 
                marmoset = gene_list_marmoset)
species_order <- c("pig", "mouse", "marmoset", "human", "crab_macaque")
# Generate the UpSet plot
p <- upset(
  fromList(gene_df),
  nsets = 6,
  keep.order = TRUE,  # Keep the specified order in `sets`
  sets = species_order,  # Explicitly set the species order here
  order.by = "degree",  # Order sets by degree
  mainbar.y.label = "Up-regulated Genes in Primed State",  # Clearer label for the main bar
  main.bar.color = "#A3C1DA",  # Use a more striking color for the main bars
  sets.bar.color = "#264653",  # Contrasting color for set size bars
  point.size = 2.0,  # Increase point size for better visibility
  sets.x.label = "Species",  # Capitalize label for consistency
  line.size = 1.2,  # Thicker lines for clarity
  text.scale = 1.6,  # Larger text for readability
  nintersects = NA,
  number.angles = 360,
  queries = list(list(query = intersects, 
                      params = list("human", "mouse", "marmoset", "pig", "crab_macaque"),
                      active = TRUE, color="#00A7A7"))
)
p
png("upsetR_n_p_up.png", res = 300, width = 3400, height = 1500)
print(p)
dev.off()

# Save the plot as SVG
svg("upsetR_n_p_UP.svg", width = 11.33, height = 5)  # Adjust width and height as needed
print(p)
dev.off()

common_genes <- Reduce(intersect, list(
  gene_list_human,
  gene_list_mouse,
  gene_list_pig,
  gene_list_crab_macaque,
  gene_list_marmoset
))
common_genes

#upset figure: down
gene_list_human <- down_human[down_human != ""]
gene_list_mouse <- down_mouse[down_mouse != ""]
gene_list_pig <- down_pig[down_pig != ""]
gene_list_crab_macaque <- down_crab_macaque[down_crab_macaque != ""]
gene_list_marmoset <- down_marmoset[down_marmoset != ""]
gene_df <- list(human = gene_list_human, mouse = toupper(gene_list_mouse), 
                pig = gene_list_pig, crab_macaque = gene_list_crab_macaque, 
                marmoset = gene_list_marmoset)
species_order <- c("pig", "mouse", "marmoset", "human", "crab_macaque")
p <- upset(
  fromList(gene_df),
  nsets = 6,
  keep.order = T,
  sets = species_order,
  order.by = "degree",  # Order sets by degree
  mainbar.y.label = "Up-regulated Genes in Naïve State",  # Clearer label for the main bar
  main.bar.color = "#A3C1DA",  # Use a more striking color for the main bars (Coral)
  sets.bar.color = "#264653",  # Use a contrasting color for set size bars (Forest Green)
  point.size = 2.0,  # Increase point size for better visibility
  sets.x.label = "Species",  # Capitalize the label for consistency
  line.size = 1.2,  # Slightly thicker lines for better clarity
  text.scale = 1.6,  # Larger text for improved readability in publications
  nintersects = NA,
  number.angles = 360,
  queries = list(list(query = intersects, 
                      params = list("human","mouse","marmoset","pig","crab_macaque"),
                      active = T,color="#00A7A7"))
)
p

# Display the plot
png("upsetR_n_p_down.png", res = 300, width = 3400, height = 1500)
print(p)
dev.off()

# Save the plot as SVG
svg("upsetR_n_p_down.svg", width = 11.33, height = 5)  # Adjust width and height as needed
print(p)
dev.off()

common_genes <- Reduce(intersect, list(
  gene_list_human,
  gene_list_mouse,
  gene_list_pig,
  gene_list_crab_macaque,
  gene_list_marmoset
))
common_genes
