# Set working directory
setwd("/Users/yf358/Desktop/PR1/Final_version/transcriptome_amalgamation/")

# Load orthologous gene data
orthology <- read.csv("orthologous_filter-9species.csv", header = TRUE, row.names = 1)

# Function to merge data
merge_species_data <- function(orthology, species_file, id_col) {
  species_data <- read.table(species_file, header = TRUE, row.names = 1)
  species_data[[id_col]] <- rownames(species_data)
  merge(orthology, species_data, by = id_col, all.x = TRUE)
}

# File paths and gene ID columns for each species
species_files <- list(
  cow = list(file = "./raw/v1_raw_cattle.txt", id_col = "Cow.gene.stable.ID"),
  human = list(file = "./raw/v1_raw_human.txt", id_col = "Gene.stable.ID"),
  mouse = list(file = "./raw/v1_raw_mouse.txt", id_col = "Mouse.gene.stable.ID"),
  pig = list(file = "./raw/v1_raw_pig.txt", id_col = "Pig.gene.stable.ID"),
  sheep = list(file = "./raw/v1_raw_sheep.txt", id_col = "Sheep.gene.stable.ID"),
  horse = list(file = "./raw/v1_raw_horse.txt", id_col = "Horse.gene.stable.ID"),
  rat = list(file = "./raw/v1_raw_rat.txt", id_col = "Rat.gene.stable.ID"),
  macaque = list(file = "./raw/v1_raw_crab_macaque.txt", id_col = "Crab.eating.macaque.gene.stable.ID"),
  marmoset = list(file = "./raw/v1_raw_marmoset.txt", id_col = "White.tufted.ear.marmoset.gene.stable.ID")
)

# Merge raw count data across species
raw_count_data <- orthology
for (species in names(species_files)) {
  raw_count_data <- merge_species_data(
    raw_count_data,
    species_files[[species]]$file,
    species_files[[species]]$id_col
  )
}
# Convert all numeric columns to integers
raw_count_data[, -(1:17)] <- apply(raw_count_data[, -(1:17)], c(1, 2), as.integer)
raw_count_data <- na.omit(raw_count_data)
# Save raw count data again
write.table(raw_count_data, "v1_rawcount.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)

# Merge TPM data
tpm_data <- orthology
species_files <- list(
  cow = list(file = "./raw/v1_tpm_cattle.txt", id_col = "Cow.gene.stable.ID"),
  human = list(file = "./raw/v1_tpm_human.txt", id_col = "Gene.stable.ID"),
  mouse = list(file = "./raw/v1_tpm_mouse.txt", id_col = "Mouse.gene.stable.ID"),
  pig = list(file = "./raw/v1_tpm_pig.txt", id_col = "Pig.gene.stable.ID"),
  sheep = list(file = "./raw/v1_tpm_sheep.txt", id_col = "Sheep.gene.stable.ID"),
  horse = list(file = "./raw/v1_tpm_horse.txt", id_col = "Horse.gene.stable.ID"),
  rat = list(file = "./raw/v1_tpm_rat.txt", id_col = "Rat.gene.stable.ID"),
  macaque = list(file = "./raw/v1_tpm_crab_macaque.txt", id_col = "Crab.eating.macaque.gene.stable.ID"),
  marmoset = list(file = "./raw/v1_tpm_marmoset.txt", id_col = "White.tufted.ear.marmoset.gene.stable.ID")
)
for (species in names(species_files)) {
  tpm_data <- merge_species_data(
    tpm_data,
    species_files[[species]]$file,
    species_files[[species]]$id_col
  )
}
tpm_data <- na.omit(tpm_data)
# Save TPM data
write.table(tpm_data, "v1_tpm.txt", row.names = FALSE, col.names = TRUE, quote = FALSE)


# SVA-log-TPM normalization
sampleinfo <- read.table("v1_sampleinfo.txt", header = T)
sampleinfo$type[grepl("ex", sampleinfo$type)] <- "expanded/extended"
sampleinfo$type[grepl("na", sampleinfo$type)] <- "naïve"
table(sampleinfo$type)

tpm <- read.table("v1_tpm.txt", header = T)
rownames(tpm) <- tpm$Gene.stable.ID
colnames(tpm)
gene_id <- tpm[,c(1:17)]
tpm <- tpm[,-c(1:17)]
tpm <- log2(tpm+1)

library(sva)
cleanY = function(y, mod, svs) {
  X = cbind(mod, svs)
  Hat = solve(t(X) %*% X) %*% t(X)
  beta = (Hat %*% t(y))
  P = ncol(mod)
  return(y - t(as.matrix(X[,-c(1:P)]) %*% beta[-c(1:P),]))
}
colnames(tpm)
rownames(sampleinfo) <- colnames(tpm)
mod = model.matrix(~ as.factor(species) + as.factor(type), data=sampleinfo)
mod0 = model.matrix(~ 1,data=sampleinfo)

# filter out the genes that are 0 in all samples, otherwise there will be an error
exp_mean <- apply(tpm,1,mean)
exp_mean2 <- tpm[which(exp_mean!=0),]

sva1 = sva(dat=as.matrix(exp_mean2), mod=mod, mod0=mod0, B=10)
tc = cleanY(y=tpm, mod=mod, svs=sva1$sv)
tc <- cbind(gene_id, tc)
write.table(tc,"v1_sva_log_tpm.txt", col.names = T, row.names = T, quote = F)
