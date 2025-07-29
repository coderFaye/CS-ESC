setwd("/Users/yf358/Desktop/PR1/round4/Figure2")
# Load necessary library
library(DESeq2)
#BiocManager::install("IHW")
library(IHW)
library(biomaRt)
#BiocManager::install("gprofiler2")
library(gprofiler2)
library(RColorBrewer)
library(pheatmap)
#install.packages("UpSetR")
library(UpSetR)

# identify plruipotency types for comparisons
com_types <- c("naive", "primed")


### Human ###
species <- "human"
sampleinfo_human <- read.table("./04.16_v1/v1_sampleinfo_human.txt", header = T)
select_indices <- sampleinfo_human$type %in% com_types
condition_human <- factor(paste0(sampleinfo_human$species, "_", sampleinfo_human$type)[select_indices])
condition_human

raw_human <- read.table("./04.16_v1/v1_raw_human.txt", header = T)
raw_human <- raw_human[,select_indices]
colnames(raw_human)
integer_counts_human <- apply(raw_human, c(1, 2), as.integer)
colData_human <- data.frame(row.names = colnames(integer_counts_human), condition_human)
dds_human <- DESeqDataSetFromMatrix(integer_counts_human, colData_human, 
                                    design = ~condition_human)

# pre-filtering
condition_counts <- table(condition_human)
smallestGroupSize <- min(condition_counts)
smallestGroupSize
keep <- rowSums(counts(dds_human) >= 10) >= smallestGroupSize
dds_human <- dds_human[keep,]
dds_norm_human <- DESeq(dds_human)
png("MA_plot_before_shrinkage_n_p_human.png", res = 300, width = 1400, height = 1400)
plotMA(dds_norm_human)
dev.off()

# get result
contrast = c("condition_human", "human_primed", "human_naive")
resIHW_human <- results(dds_norm_human, contrast = contrast, filterFun=ihw)
summary(resIHW_human)
ihw_df_human <- as.data.frame(resIHW_human)
write.csv(ihw_df_human, "DEG_results_n_p_human.csv")

# get significant degs
p_cutoff <- 0.05
logFC_cutoff <- 1
ihw_df_human <- na.omit(ihw_df_human)
significant_degs <- ihw_df_human[ihw_df_human$padj < p_cutoff & abs(ihw_df_human$log2FoldChange) > logFC_cutoff, ]
human_n_p <- significant_degs[order(significant_degs$log2FoldChange, decreasing = TRUE), ]
head(human_n_p)
write.csv(human_n_p, "significant_DEG_n_p_human.csv")
sum(significant_degs$log2FoldChange > 0)
sum(significant_degs$log2FoldChange < 0)
human_n_p$cluster <- rep(1,nrow(human_n_p))

# lfc-shrinkage
resultsNames(dds_norm_human)
resLFC_human <- lfcShrink(dds_norm_human, coef=2, type="apeglm", svalue = F, lfcThreshold = 0)
resLFC_df_human <- as.data.frame(resLFC_human)
write.csv(resLFC_df_human, "lfc_n_p_human.csv")
png("MA_plot_shrinkage_n_p_human.png", res = 300, width = 1400, height = 1400)
plotMA(resLFC_human)
dev.off()

# go analysis
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
symbol  <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                 filters = 'ensembl_gene_id',
                 values = rownames(human_n_p),
                 mart = ensembl)
human_n_p$ensembl_gene_id <- rownames(human_n_p)
human_n_p <- merge(human_n_p, symbol, by = "ensembl_gene_id")

gene_list_human <- human_n_p$external_gene_name
head(gene_list_human)
up_human <- human_n_p[which(human_n_p$log2FoldChange > 0), "external_gene_name"]
down_human <- human_n_p[which(human_n_p$log2FoldChange < 0), "external_gene_name"]

# enrich upregulated genes
pathway_up_human <- gost(up_human,
                         organism = "hsapiens",
                         ordered_query = FALSE,
                         multi_query = FALSE,
                         significant = FALSE,
                         exclude_iea = FALSE,
                         measure_underrepresentation = FALSE,
                         evcodes = FALSE,
                         user_threshold = 0.05,
                         correction_method = "g_SCS",
                         domain_scope = "annotated",
                         custom_bg = NULL,
                         numeric_ns = "",
                         sources = NULL,
                         as_short_link = FALSE,
                         highlight = FALSE
)
gostplot(pathway_up_human, capped = FALSE, interactive = FALSE)
png("pathway_up_human_n_p.png", width = 2900, height = 1500, res = 300)
gostplot(pathway_up_human, capped = FALSE, interactive = FALSE)
dev.off()

pathway_up_human <- as.data.frame(pathway_up_human[["result"]])
pathway_up_human$parents <- sapply(pathway_up_human$parents, function(x) paste(x, collapse = ", "))
write.csv(pathway_up_human, "human_pathway_n_p_up.csv")
go_up_human <- pathway_up_human[which(pathway_up_human$source=="GO:BP"),]
write.csv(go_up_human, "human_go_n_p_up.csv")

# enrich down regulated genes
pathway_down_human <- gost(down_human,
                           organism = "hsapiens",
                           ordered_query = FALSE,
                           multi_query = FALSE,
                           significant = FALSE,
                           exclude_iea = FALSE,
                           measure_underrepresentation = FALSE,
                           evcodes = FALSE,
                           user_threshold = 0.05,
                           correction_method = "g_SCS",
                           domain_scope = "annotated",
                           custom_bg = NULL,
                           numeric_ns = "",
                           sources = NULL,
                           as_short_link = FALSE,
                           highlight = FALSE
)
gostplot(pathway_down_human, capped = FALSE, interactive = FALSE)
png("pathway_down_human_n_p.png", width = 2900, height = 1500, res = 300)
gostplot(pathway_down_human, capped = FALSE, interactive = FALSE)
dev.off()
pathway_down_human <- as.data.frame(pathway_down_human[["result"]])
pathway_down_human$parents <- sapply(pathway_down_human$parents, function(x) paste(x, collapse = ", "))
write.csv(pathway_down_human, "human_pathway_n_p_down.csv")
go_down_human <- pathway_down_human[which(pathway_down_human$source=="GO:BP"),]
write.csv(go_down_human, "human_go_n_p_down.csv")

       
### Mouse ###
species <- "mouse"
sampleinfo_mouse <- read.table("./04.16_v1/v1_sampleinfo_mouse.txt", header = T)
select_indices <- sampleinfo_mouse$type %in% com_types
condition_mouse <- factor(paste0(sampleinfo_mouse$species, "_", sampleinfo_mouse$type)[select_indices])
condition_mouse

raw_mouse <- read.table("./04.16_v1/v1_raw_mouse.txt", header = T)
raw_mouse <- raw_mouse[,select_indices]
colnames(raw_mouse)
integer_counts_mouse <- apply(raw_mouse, c(1, 2), as.integer)
colData_mouse <- data.frame(row.names = colnames(integer_counts_mouse), condition_mouse)
dds_mouse <- DESeqDataSetFromMatrix(integer_counts_mouse, colData_mouse, 
                                    design = ~condition_mouse)
#pre-filtering
condition_counts <- table(condition_mouse)
smallestGroupSize <- min(condition_counts)
smallestGroupSize
keep <- rowSums(counts(dds_mouse) >= 10) >= smallestGroupSize
dds_mouse <- dds_mouse[keep,]
dds_norm_mouse <- DESeq(dds_mouse)
png("MA_plot_before_shrinkage_n_p_mouse.png", res = 300, width = 1400, height = 1400)
plotMA(dds_norm_mouse)
dev.off()

# get result
contrast = c("condition_mouse", "mouse_primed", "mouse_naive")
resIHW_mouse <- results(dds_norm_mouse, contrast = contrast, filterFun=ihw)
summary(resIHW_mouse)
ihw_df_mouse <- as.data.frame(resIHW_mouse)
write.csv(ihw_df_mouse, "DEG_results_n_p_mouse.csv")

# get significant degs
p_cutoff <- 0.05
logFC_cutoff <- 1
ihw_df_mouse <- na.omit(ihw_df_mouse)
significant_degs <- ihw_df_mouse[ihw_df_mouse$padj < p_cutoff & abs(ihw_df_mouse$log2FoldChange) > logFC_cutoff, ]
sum(significant_degs$log2FoldChange > 0)
sum(significant_degs$log2FoldChange < 0)
mouse_n_p <- significant_degs[order(significant_degs$log2FoldChange, decreasing = TRUE), ]
write.csv(mouse_n_p, "significant_DEG_n_p_mouse.csv")
mouse_n_p$cluster <- rep(1,nrow(mouse_n_p))

#lfc-shrinkage
resultsNames(dds_norm_mouse)
resLFC_mouse <- lfcShrink(dds_norm_mouse, coef=2, type="apeglm", svalue = F, lfcThreshold = 0)
summary(resLFC_mouse)
resLFC_df_mouse <- as.data.frame(resLFC_mouse)
write.csv(resLFC_df_mouse, "lfc_n_p_mouse.csv")
png("MA_plot_shrinkage_n_p_mouse.png", res = 300, width = 1400, height = 1400)
plotMA(resLFC_mouse)
dev.off()
plotMA(resLFC_mouse)

# go analysis
ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
symbol  <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                 filters = 'ensembl_gene_id',
                 values = rownames(mouse_n_p),
                 mart = ensembl)
mouse_n_p$ensembl_gene_id <- rownames(mouse_n_p)
mouse_n_p <- merge(mouse_n_p, symbol, by = "ensembl_gene_id")
gene_list_mouse <- mouse_n_p$external_gene_name
head(gene_list_mouse)
up_mouse <- mouse_n_p[which(mouse_n_p$log2FoldChange > 0), "external_gene_name"]
down_mouse <- mouse_n_p[which(mouse_n_p$log2FoldChange < 0), "external_gene_name"]

pathway_up_mouse <- gost(up_mouse,
                         organism = "mmusculus",
                         ordered_query = FALSE,
                         multi_query = FALSE,
                         significant = FALSE,
                         exclude_iea = FALSE,
                         measure_underrepresentation = FALSE,
                         evcodes = FALSE,
                         user_threshold = 0.05,
                         correction_method = "g_SCS",
                         domain_scope = "annotated",
                         custom_bg = NULL,
                         numeric_ns = "",
                         sources = NULL,
                         as_short_link = FALSE,
                         highlight = FALSE
)
gostplot(pathway_up_mouse, capped = FALSE, interactive = FALSE)
png("pathway_up_mouse_n_p.png", width = 2900, height = 1500, res = 300)
gostplot(pathway_up_mouse, capped = FALSE, interactive = FALSE)
dev.off()
pathway_up_mouse <- as.data.frame(pathway_up_mouse[["result"]])
pathway_up_mouse$parents <- sapply(pathway_up_mouse$parents, function(x) paste(x, collapse = ", "))
write.csv(pathway_up_mouse, "mouse_pathway_n_p_up.csv")
go_up_mouse <- pathway_up_mouse[which(pathway_up_mouse$source=="GO:BP"),]
write.csv(go_up_mouse, "mouse_go_n_p_up.csv")

pathway_down_mouse <- gost(down_mouse,
                           organism = "mmusculus",
                           ordered_query = FALSE,
                           multi_query = FALSE,
                           significant = FALSE,
                           exclude_iea = FALSE,
                           measure_underrepresentation = FALSE,
                           evcodes = FALSE,
                           user_threshold = 0.05,
                           correction_method = "g_SCS",
                           domain_scope = "annotated",
                           custom_bg = NULL,
                           numeric_ns = "",
                           sources = NULL,
                           as_short_link = FALSE,
                           highlight = FALSE
)
gostplot(pathway_down_mouse, capped = FALSE, interactive = FALSE)
png("pathway_down_mouse_n_p.png", width = 2900, height = 1500, res = 300)
gostplot(pathway_down_mouse, capped = FALSE, interactive = FALSE)
dev.off()
pathway_down_mouse <- as.data.frame(pathway_down_mouse[["result"]])
pathway_down_mouse$parents <- sapply(pathway_down_mouse$parents, function(x) paste(x, collapse = ", "))
write.csv(pathway_down_mouse, "mouse_pathway_n_p_down.csv")
go_down_mouse <- pathway_down_mouse[which(pathway_down_mouse$source=="GO:BP"),]
write.csv(go_down_mouse, "mouse_go_n_p_down.csv")

                                     
### Pig ###
species <- "pig"
sampleinfo_pig <- read.table("./04.16_v1/v1_sampleinfo_pig.txt", header = T)
select_indices <- sampleinfo_pig$type %in% com_types
condition_pig <- factor(paste0(sampleinfo_pig$species, "_", sampleinfo_pig$type)[select_indices])
condition_pig

raw_pig <- read.table("./04.16_v1/v1_raw_pig.txt", header = T)
raw_pig <- raw_pig[,select_indices]
colnames(raw_pig)
integer_counts_pig <- apply(raw_pig, c(1, 2), as.integer)
colData_pig <- data.frame(row.names = colnames(integer_counts_pig), condition_pig)
dds_pig <- DESeqDataSetFromMatrix(integer_counts_pig, colData_pig, 
                                  design = ~condition_pig)

#pre-filtering
condition_counts <- table(condition_pig)
smallestGroupSize <- min(condition_counts)
smallestGroupSize
keep <- rowSums(counts(dds_pig) >= 10) >= smallestGroupSize
dds_pig <- dds_pig[keep,]
dds_norm_pig <- DESeq(dds_pig)
png("MA_plot_before_shrinkage_n_p_pig.png", res = 300, width = 1400, height = 1400)
plotMA(dds_norm_pig)
dev.off()

# get result
contrast = c("condition_pig", "pig_primed", "pig_naive")
resIHW_pig <- results(dds_norm_pig, contrast = contrast, filterFun=ihw)
summary(resIHW_pig)
ihw_df_pig <- as.data.frame(resIHW_pig)
write.csv(ihw_df_pig, "DEG_results_n_p_pig.csv")

# get significant degs
p_cutoff <- 0.05
logFC_cutoff <- 1
ihw_df_pig <- na.omit(ihw_df_pig)
significant_degs <- ihw_df_pig[ihw_df_pig$padj < p_cutoff & abs(ihw_df_pig$log2FoldChange) > logFC_cutoff, ]
sum(significant_degs$log2FoldChange > 0)
sum(significant_degs$log2FoldChange < 0)
pig_n_p <- significant_degs[order(significant_degs$log2FoldChange, decreasing = TRUE), ]
write.csv(pig_n_p, "significant_DEG_n_p_pig.csv")
pig_n_p$cluster <- rep(1,nrow(pig_n_p))

# lfc-shrinkage
resultsNames(dds_norm_pig)
resLFC_pig <- lfcShrink(dds_norm_pig, coef=2, type="apeglm", svalue = F, lfcThreshold = 0)
png("MA_plot_shrinkage_n_p_pig.png", res = 300, width = 1400, height = 1400)
plotMA(resLFC_pig) 
dev.off()
resLFC_df_pig <- as.data.frame(resLFC_pig)
write.csv(resLFC_df_pig, "lfc_n_p_pig.csv")
plotMA(resLFC_pig)

ensembl <- useMart("ensembl", dataset = "sscrofa_gene_ensembl")
symbol  <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                 filters = 'ensembl_gene_id',
                 values = rownames(pig_n_p),
                 mart = ensembl)
pig_n_p$ensembl_gene_id <- rownames(pig_n_p)
pig_n_p <- merge(pig_n_p, symbol, by = "ensembl_gene_id")
gene_list_pig <- pig_n_p$external_gene_name
up_pig <- pig_n_p[which(pig_n_p$log2FoldChange > 0), "external_gene_name"]
down_pig <- pig_n_p[which(pig_n_p$log2FoldChange < 0), "external_gene_name"]
head(gene_list_pig)

pathway_up_pig <- gost(up_pig,
                       organism = "sscrofa",
                       ordered_query = FALSE,
                       multi_query = FALSE,
                       significant = FALSE,
                       exclude_iea = FALSE,
                       measure_underrepresentation = FALSE,
                       evcodes = FALSE,
                       user_threshold = 0.05,
                       correction_method = "g_SCS",
                       domain_scope = "annotated",
                       custom_bg = NULL,
                       numeric_ns = "",
                       sources = NULL,
                       as_short_link = FALSE,
                       highlight = FALSE
)
gostplot(pathway_up_pig, capped = FALSE, interactive = FALSE)
png("pathway_up_pig_n_p.png", width = 2900, height = 1500, res = 300)
gostplot(pathway_up_pig, capped = FALSE, interactive = FALSE)
dev.off()
pathway_up_pig <- as.data.frame(pathway_up_pig[["result"]])
pathway_up_pig$parents <- sapply(pathway_up_pig$parents, function(x) paste(x, collapse = ", "))
write.csv(pathway_up_pig, "pig_pathway_n_p_up.csv")
go_up_pig <- pathway_up_pig[which(pathway_up_pig$source=="GO:BP"),]
write.csv(go_up_pig, "pig_go_n_p_up.csv")

pathway_down_pig <- gost(down_pig,
                         organism = "sscrofa",
                         ordered_query = FALSE,
                         multi_query = FALSE,
                         significant = FALSE,
                         exclude_iea = FALSE,
                         measure_underrepresentation = FALSE,
                         evcodes = FALSE,
                         user_threshold = 0.05,
                         correction_method = "g_SCS",
                         domain_scope = "annotated",
                         custom_bg = NULL,
                         numeric_ns = "",
                         sources = NULL,
                         as_short_link = FALSE,
                         highlight = FALSE
)
gostplot(pathway_down_pig, capped = FALSE, interactive = FALSE)
png("pathway_down_pig_n_p.png", width = 2900, height = 1500, res = 300)
gostplot(pathway_down_pig, capped = FALSE, interactive = FALSE)
dev.off()
pathway_down_pig <- as.data.frame(pathway_down_pig[["result"]])
pathway_down_pig$parents <- sapply(pathway_down_pig$parents, function(x) paste(x, collapse = ", "))
write.csv(pathway_down_pig, "pig_pathway_n_p_down.csv")
go_down_pig <- pathway_down_pig[which(pathway_down_pig$source=="GO:BP"),]
write.csv(go_down_pig, "pig_go_n_p_down.csv")

                                   
### Marmoset ###
species <- "marmoset"
sampleinfo_marmoset <- read.table("./04.16_v1/v1_sampleinfo_marmoset.txt", header = T)
sampleinfo_marmoset$type[grepl("na", sampleinfo_marmoset$type)] <- "naive"
select_indices <- sampleinfo_marmoset$type %in% com_types
condition_marmoset <- factor(paste0(sampleinfo_marmoset$species, "_", sampleinfo_marmoset$type)[select_indices])
condition_marmoset

raw_marmoset <- read.table("./04.16_v1/v1_raw_marmoset.txt", header = T)
raw_marmoset <- raw_marmoset[,select_indices]
colnames(raw_marmoset)
integer_counts_marmoset <- apply(raw_marmoset, c(1, 2), as.integer)
colData_marmoset <- data.frame(row.names = colnames(integer_counts_marmoset), 
                               condition_marmoset)
dds_marmoset <- DESeqDataSetFromMatrix(integer_counts_marmoset, colData_marmoset, 
                                       design = ~condition_marmoset)

#pre-filtering
condition_counts <- table(condition_marmoset)
smallestGroupSize <- min(condition_counts)
smallestGroupSize
keep <- rowSums(counts(dds_marmoset) >= 10) >= smallestGroupSize
dds_marmoset <- dds_marmoset[keep,]
dds_norm_marmoset <- DESeq(dds_marmoset)

png("MA_plot_before_shrinkage_n_p_marmoset.png", res = 300, width = 1400, height = 1400)
plotMA(dds_norm_marmoset)
dev.off()

# get result
contrast = c("condition_marmoset", "marmoset_primed", "marmoset_naive")
resIHW_marmoset <- results(dds_norm_marmoset, contrast = contrast, filterFun=ihw)
summary(resIHW_marmoset)
ihw_df_marmoset <- as.data.frame(resIHW_marmoset)
write.csv(ihw_df_marmoset, "DEG_results_n_p_marmoset.csv")

# get significant degs
p_cutoff <- 0.05
logFC_cutoff <- 1
ihw_df_marmoset <- na.omit(ihw_df_marmoset)
significant_degs <- ihw_df_marmoset[ihw_df_marmoset$padj < p_cutoff & abs(ihw_df_marmoset$log2FoldChange) > logFC_cutoff, ]
sum(significant_degs$log2FoldChange > 0)
sum(significant_degs$log2FoldChange < 0)
marmoset_n_p <- significant_degs[order(significant_degs$log2FoldChange, decreasing = TRUE), ]
write.csv(marmoset_n_p, "significant_DEG_n_p_marmoset.csv")
marmoset_n_p$cluster <- rep(1,nrow(marmoset_n_p))

#lfc-shrinkage
resultsNames(dds_norm_marmoset)
resLFC_marmoset <- lfcShrink(dds_norm_marmoset, coef=2, type="apeglm", svalue = F, lfcThreshold = 0)
resLFC_df_marmoset <- as.data.frame(resLFC_marmoset)
write.csv(resLFC_df_marmoset, "lfc_n_p_marmoset.csv")
png("MA_plot_shrinkage_n_p_marmoset.png", res = 300, width = 1400, height = 1400)
plotMA(resLFC_marmoset) 
dev.off()

# pathway enrichment
ensembl <- useMart("ensembl", dataset = "cjacchus_gene_ensembl")
symbol  <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                 filters = 'ensembl_gene_id',
                 values = rownames(marmoset_n_p),
                 mart = ensembl)
marmoset_n_p$ensembl_gene_id <- rownames(marmoset_n_p)
marmoset_n_p <- merge(marmoset_n_p, symbol, by = "ensembl_gene_id")
gene_list_marmoset <- marmoset_n_p$external_gene_name
up_marmoset <- marmoset_n_p[which(marmoset_n_p$log2FoldChange > 0), "external_gene_name"]
down_marmoset <- marmoset_n_p[which(marmoset_n_p$log2FoldChange < 0), "external_gene_name"]
head(gene_list_marmoset)

up_marmoset <- up_marmoset[up_marmoset != ""]
up_marmoset <- up_marmoset[up_marmoset != "Metazoa_SRP"]
pathway_up_marmoset <- gost(up_marmoset,
                            organism = "cjacchus",
                            ordered_query = FALSE,
                            multi_query = FALSE,
                            significant = FALSE,
                            exclude_iea = FALSE,
                            measure_underrepresentation = FALSE,
                            evcodes = FALSE,
                            user_threshold = 0.05,
                            correction_method = "g_SCS",
                            domain_scope = "annotated",
                            custom_bg = NULL,
                            numeric_ns = "",
                            sources = NULL,
                            as_short_link = FALSE,
                            highlight = FALSE
)
gostplot(pathway_up_marmoset, capped = FALSE, interactive = FALSE)
png("pathway_up_marmoset_n_p.png", width = 2900, height = 1500, res = 300)
gostplot(pathway_up_marmoset, capped = FALSE, interactive = FALSE)
dev.off()
pathway_up_marmoset <- as.data.frame(pathway_up_marmoset[["result"]])
pathway_up_marmoset$parents <- sapply(pathway_up_marmoset$parents, function(x) paste(x, collapse = ", "))
write.csv(pathway_up_marmoset, "marmoset_pathway_n_p_up.csv")
go_up_marmoset <- pathway_up_marmoset[which(pathway_up_marmoset$source=="GO:BP"),]
write.csv(go_up_marmoset, "marmoset_go_n_p_up.csv")

down_marmoset <- down_marmoset[down_marmoset != ""]
down_marmoset <- down_marmoset[down_marmoset != "Metazoa_SRP"]
pathway_down_marmoset <- gost(down_marmoset,
                              organism = "cjacchus",
                              ordered_query = FALSE,
                              multi_query = FALSE,
                              significant = FALSE,
                              exclude_iea = FALSE,
                              measure_underrepresentation = FALSE,
                              evcodes = FALSE,
                              user_threshold = 0.05,
                              correction_method = "g_SCS",
                              domain_scope = "annotated",
                              custom_bg = NULL,
                              numeric_ns = "",
                              sources = NULL,
                              as_short_link = FALSE,
                              highlight = FALSE
)
gostplot(pathway_down_marmoset, capped = FALSE, interactive = FALSE)
png("pathway_down_marmoset_n_p.png", width = 2900, height = 1500, res = 300)
gostplot(pathway_down_marmoset, capped = FALSE, interactive = FALSE)
dev.off()
pathway_down_marmoset <- as.data.frame(pathway_down_marmoset[["result"]])
pathway_down_marmoset$parents <- sapply(pathway_down_marmoset$parents, function(x) paste(x, collapse = ", "))
write.csv(pathway_down_marmoset, "marmoset_pathway_n_p_down.csv")
go_down_marmoset <- pathway_down_marmoset[which(pathway_down_marmoset$source=="GO:BP"),]
write.csv(go_down_marmoset, "marmoset_go_n_p_down.csv")


### Crab_macaque ###
species <- "crab_macaque"
sampleinfo_crab_macaque <- read.table("./04.16_v1/v1_sampleinfo_crab_macaque.txt", header = T)
sampleinfo_crab_macaque$type[grepl("na", sampleinfo_crab_macaque$type)] <- "naive"
select_indices <- sampleinfo_crab_macaque$type %in% com_types
condition_crab_macaque <- factor(paste0(sampleinfo_crab_macaque$species, "_", sampleinfo_crab_macaque$type)[select_indices])
condition_crab_macaque

raw_crab_macaque <- read.table("./04.16_v1/v1_raw_crab_macaque.txt", header = T)
raw_crab_macaque <- raw_crab_macaque[,select_indices]
colnames(raw_crab_macaque)
integer_counts_crab_macaque <- apply(raw_crab_macaque, c(1, 2), as.integer)
colData_crab_macaque <- data.frame(row.names = colnames(integer_counts_crab_macaque), 
                                   condition_crab_macaque)
dds_crab_macaque <- DESeqDataSetFromMatrix(integer_counts_crab_macaque, colData_crab_macaque, 
                                           design = ~condition_crab_macaque)

#pre-filtering
condition_counts <- table(condition_crab_macaque)
smallestGroupSize <- min(condition_counts)
smallestGroupSize
keep <- rowSums(counts(dds_crab_macaque) >= 10) >= smallestGroupSize
dds_crab_macaque <- dds_crab_macaque[keep,]
dds_norm_crab_macaque <- DESeq(dds_crab_macaque)

png("MA_plot_before_shrinkage_n_p_crab_macaque.png", res = 300, width = 1400, height = 1400)
plotMA(dds_norm_crab_macaque)
dev.off()

# get result
contrast = c("condition_crab_macaque", "crab_macaque_primed", "crab_macaque_naive")
resIHW_crab_macaque <- results(dds_norm_crab_macaque, contrast = contrast, filterFun=ihw)
summary(resIHW_crab_macaque)
ihw_df_crab_macaque <- as.data.frame(resIHW_crab_macaque)
write.csv(ihw_df_crab_macaque, "DEG_results_n_p_crab_macaque.csv")

# get significant degs
p_cutoff <- 0.05
logFC_cutoff <- 1
ihw_df_crab_macaque <- na.omit(ihw_df_crab_macaque)
significant_degs <- ihw_df_crab_macaque[ihw_df_crab_macaque$padj < p_cutoff & abs(ihw_df_crab_macaque$log2FoldChange) > logFC_cutoff, ]
sum(significant_degs$log2FoldChange > 0)
sum(significant_degs$log2FoldChange < 0)
crab_macaque_n_p <- significant_degs[order(significant_degs$log2FoldChange, decreasing = TRUE), ]
write.csv(crab_macaque_n_p, "significant_DEG_n_p_crab_macaque.csv")
crab_macaque_n_p$cluster <- rep(1,nrow(crab_macaque_n_p))

#lfc-shrinkage
resultsNames(dds_norm_crab_macaque)
resLFC_crab_macaque <- lfcShrink(dds_norm_crab_macaque, coef=2, type="apeglm", svalue = F, lfcThreshold = 0)
resLFC_df_crab_macaque <- as.data.frame(resLFC_crab_macaque)
write.csv(resLFC_df_crab_macaque, "lfc_n_p_crab_macaque.csv")
png("MA_plot_shrinkage_n_p_crab_macaque.png", res = 300, width = 1400, height = 1400)
plotMA(resLFC_crab_macaque) 
dev.off()

# pathway enrichment
ensembl <- useMart("ensembl", dataset = "mfascicularis_gene_ensembl")
symbol  <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'),
                 filters = 'ensembl_gene_id',
                 values = rownames(crab_macaque_n_p),
                 mart = ensembl)
crab_macaque_n_p$ensembl_gene_id <- rownames(crab_macaque_n_p)
crab_macaque_n_p <- merge(crab_macaque_n_p, symbol, by = "ensembl_gene_id")
gene_list_crab_macaque <- crab_macaque_n_p$external_gene_name
up_crab_macaque <- crab_macaque_n_p[which(crab_macaque_n_p$log2FoldChange > 0), "external_gene_name"]
down_crab_macaque <- crab_macaque_n_p[which(crab_macaque_n_p$log2FoldChange < 0), "external_gene_name"]
head(gene_list_crab_macaque)

pathway_up_crab_macaque <- gost(up_crab_macaque,
                                organism = "mfascicularis",
                                ordered_query = FALSE,
                                multi_query = FALSE,
                                significant = FALSE,
                                exclude_iea = FALSE,
                                measure_underrepresentation = FALSE,
                                evcodes = FALSE,
                                user_threshold = 0.05,
                                correction_method = "g_SCS",
                                domain_scope = "annotated",
                                custom_bg = NULL,
                                numeric_ns = "",
                                sources = NULL,
                                as_short_link = FALSE,
                                highlight = FALSE
)
gostplot(pathway_up_crab_macaque, capped = FALSE, interactive = FALSE)
png("pathway_up_crab_macaque_n_p.png", width = 2900, height = 1500, res = 300)
gostplot(pathway_up_crab_macaque, capped = FALSE, interactive = FALSE)
dev.off()
pathway_up_crab_macaque <- as.data.frame(pathway_up_crab_macaque[["result"]])
pathway_up_crab_macaque$parents <- sapply(pathway_up_crab_macaque$parents, function(x) paste(x, collapse = ", "))
write.csv(pathway_up_crab_macaque, "crab_macaque_pathway_n_p_up.csv")
go_up_crab_macaque <- pathway_up_crab_macaque[which(pathway_up_crab_macaque$source=="GO:BP"),]
write.csv(go_up_crab_macaque, "crab_macaque_go_n_p_up.csv")

pathway_down_crab_macaque <- gost(down_crab_macaque,
                                  organism = "mfascicularis",
                                  ordered_query = FALSE,
                                  multi_query = FALSE,
                                  significant = FALSE,
                                  exclude_iea = FALSE,
                                  measure_underrepresentation = FALSE,
                                  evcodes = FALSE,
                                  user_threshold = 0.05,
                                  correction_method = "g_SCS",
                                  domain_scope = "annotated",
                                  custom_bg = NULL,
                                  numeric_ns = "",
                                  sources = NULL,
                                  as_short_link = FALSE,
                                  highlight = FALSE
)
gostplot(pathway_down_crab_macaque, capped = FALSE, interactive = FALSE)
png("pathway_down_crab_macaque_n_p.png", width = 2900, height = 1500, res = 300)
gostplot(pathway_down_crab_macaque, capped = FALSE, interactive = FALSE)
dev.off()
pathway_down_crab_macaque <- as.data.frame(pathway_down_crab_macaque[["result"]])
pathway_down_crab_macaque$parents <- sapply(pathway_down_crab_macaque$parents, function(x) paste(x, collapse = ", "))
write.csv(pathway_down_crab_macaque, "crab_macaque_pathway_n_p_down.csv")
go_down_crab_macaque <- pathway_down_crab_macaque[which(pathway_down_crab_macaque$source=="GO:BP"),]
write.csv(go_down_crab_macaque, "crab_macaque_go_n_p_down.csv")
save.image("n_p_deg.rds")
