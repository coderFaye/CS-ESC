setwd("/Users/yf358/Desktop/PR1/Final_version/Transcriptome_amalgamation/Data_exploratory_analysis")

sampleinfo <- read.table("/Users/yf358/Desktop/PR1/Final_version/transcriptome_amalgamation/v1_sampleinfo.txt", header = T)
sampleinfo$type[grepl("ex", sampleinfo$type)] <- "expanded/extended"
sampleinfo$type[grepl("na", sampleinfo$type)] <- "naïve"
new_condition <- paste0(sampleinfo$species,"_", sampleinfo$type)

tpm <- read.table("/Users/yf358/Desktop/PR1/Final_version/transcriptome_amalgamation/v1_sva_log_tpm.txt", header = T)
rownames(tpm) <- tpm$Gene.stable.ID
colnames(tpm)
tpm <- tpm[,-c(1:17)]


library("variancePartition")
form <- ~ type + species
varPart <- fitExtractVarPartModel(tpm, form, sampleinfo)
colnames(varPart)
colnames(varPart) <- c("Pluripotency", "Species", "Residuals")
# figure s1a
plotVarPart(varPart)
png("varPart_plot.png", width=1200, height=1000, res = 300)
plotVarPart(varPart)
dev.off()

# use gene symbol as rownames
library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
symbol  <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), 
                 filters = 'ensembl_gene_id', 
                 values = rownames(varPart), 
                 mart = ensembl)
varPart$ensembl_gene_id <- rownames(varPart)
varPart <- merge(varPart, symbol, by = "ensembl_gene_id")
rownames(varPart) <- varPart$external_gene_name
colnames(varPart)
varPart <- varPart[,-c(1,5)]


# sort by the value of type
sortedVarPart <- varPart[order(varPart$Pluripotency, decreasing = TRUE), ]
plotPercentBars(sortedVarPart[1:30, ])

sortedVarPart_1 <- varPart[order(varPart$Species, decreasing = TRUE), ]
plotPercentBars(sortedVarPart_1[1:30, ])


# figure s1b
library(gghighlight)
library(tidyverse)
varPart$name <- rownames(varPart)
varPart <- as.data.frame(varPart)
p <- varPart %>%
  ggplot(aes(x = Pluripotency, y = Species)) +
  geom_point(color = alpha("#CD443F", 0.8), size = 2.5) +  
  gghighlight(Pluripotency > 0.45 ,
              unhighlighted_colour = alpha("#C0BC99", 0.4),
              use_direct_label = TRUE,  
              label_key = name) +
  labs(x = "Proportion of gene expression variance explained by pluripotency", 
       y = "Proportion of gene expression variance explained by species",
       title = "Variance Partitioning") +
  theme_minimal(base_size = 18) +  
  theme(
    plot.title = element_text(face = "bold", size = 16),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 16),
    legend.position = "none"  
  )
p
png(filename = "Variance_Partitioning_Highlighting.png", width = 10, height = 7, units = 'in', res = 300)
print(p)
dev.off()

