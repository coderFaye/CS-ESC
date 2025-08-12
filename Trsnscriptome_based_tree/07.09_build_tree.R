#### Remaining part: bootstrapping ####
setwd("/Users/yf358/Desktop/PR1/round3/04_tree")

# read in sampleinfo and tpm
sampleinfo <- read.table("v1_sampleinfo.txt", header = T)
sampleinfo$type[grepl("ex", sampleinfo$type)] <- "expanded"
sampleinfo$species[grepl("crab", sampleinfo$species)] <- "crabmacaque"
table(sampleinfo$type)

tpm <- read.table("v1_sva_log_tpm.txt", header = T)
rownames(tpm) <- tpm$Gene.stable.ID
colnames(tpm)
tpm <- tpm[,-c(1:17)]

### BUILD TREE ###
library(ape)
library(cowplot)
library(tidyverse)
library(dplyr)
library(tidyr)
library(ggtree)
library(scales)

# rename the samples
name <- paste0(sampleinfo$species, "_", sampleinfo$type)
combined_key <- paste(sampleinfo$species, sampleinfo$type, sep="_")
counter <- ave(rep(1, length(combined_key)), combined_key, FUN = seq_along)
name <- paste0(combined_key, ".TR_RP", counter)
print(name)
colnames(tpm) <- name
colnames(tpm)

# calculate the median for each group
expression_filtered_df <- as.data.frame(tpm)
medians <- matrix(nrow = nrow(expression_filtered_df), ncol = 24)
for (i in 1:nrow(expression_filtered_df)) {
  # This assumes that 'expression_filtered_df' is a dataframe where each row should be summarized
  median_value <- expression_filtered_df[i, ] %>%
    pivot_longer(cols = colnames(expression_filtered_df), names_to = "colname", values_to = "value") %>%
    separate(colname, into = c("species", "pluripotency", "type", "rep"), sep = "(\\.|_)", remove = FALSE) %>%
    group_by(species, pluripotency) %>%
    summarise(median = median(value, na.rm = TRUE)) %>%
    ungroup() %>%
    pull(median)  # Extract the median value
  # Assign the computed median to the 'medians' matrix or data frame
  medians[i, ] <- median_value
}
group <- sampleinfo %>% group_by(species, type) %>% summarise() %>% ungroup()
name_median <- paste0("median_", group$species,"_", group$type, ".TR")
colnames(medians) <- name_median
colnames(medians)
#save(medians, file = "07_09_sva_lot_tpm_median.RData")
#load("07_09_sva_lot_tpm_median.RData")

# combine the medians with the expression matrix
expression_combined <- cbind(tpm, medians)
expression_combined <- as.data.frame(expression_combined)

# estimate of within-species vriance and measurement error ####
error_estimate = tibble( colname = colnames(expression_combined) ) %>%
  filter( ! str_detect( colname, "median" ) ) %>%
  filter( colname != "Gene_ID" ) %>%
  separate( colname, c( "species", "tissue", "type", "rep" ), "(\\.|_)", remove=FALSE  ) %>%
  group_by( species, tissue, type ) %>%
  summarise( meanVar = mean( genefilter::rowVars( expression_combined[,colname,drop=FALSE] ), na.rm=TRUE ) ) 
error_estimate = data.frame(error_estimate)

# count the sample number of each group
group_number = tibble( colname = colnames(expression_combined) ) %>%
  filter( ! str_detect( colname, "median" ) ) %>%
  filter( colname != "Gene_ID" ) %>%
  separate( colname, c( "species", "type", "rep" ), "(\\.|_)", remove=FALSE  ) %>%
  group_by( species, type ) %>% summarise(n = n())

## divide by group number
error_estimate$dividedVar = error_estimate$meanVar/group_number$n

# estimate within group variance
variance_estimates_TR = tibble( colname = colnames(expression_combined) ) %>%
  filter( str_detect( colname, "median" ) ) %>%
  filter( colname != "Gene_ID" ) %>%
  separate( "colname", c( "median", "species", "tissue", "type" ),  "(\\.|_)", remove=FALSE ) %>%
  spread( type, colname ) %>%
  group_by_all() %>%
  summarise( TR_var = var( expression_combined[ , TR ], na.rm=TRUE ) ) %>%
  ungroup
variance_estimates = data.frame(variance_estimates_TR)


### set tissue: primed ###
tissue = "primed"
group_number[which(group_number$type=="primed"), 1]
organisms_list = c("cattle","human","marmoset","mouse","pig","crabmacaque","sheep") #primed

# get distance matrix
TR_distances_with_errors = crossing( sp1 = organisms_list, sp2 = organisms_list ) %>%
  filter( sp1 != sp2 ) %>%
  rowwise() %>%
  group_by_all() %>%
  reframe(var_divergence = 
            var(expression_combined[ , paste("median_",sp1,"_",tissue,".","TR",sep = "")] - 
                  expression_combined[ , paste("median_",sp2,"_",tissue,".","TR",sep = "")], na.rm=TRUE ) - 
            error_estimate[error_estimate$species == sp1 & error_estimate$tissue == tissue & error_estimate$type == "TR", "dividedVar"] -
            error_estimate[error_estimate$species == sp2 & error_estimate$tissue == tissue & error_estimate$type == "TR", "dividedVar"]) %>%
  #ungroup %>%
  spread( sp2, var_divergence, fill=0 ) %>%
  column_to_rownames( "sp1" ) %>%
  as.matrix
TR_distances_with_errors_normalized = TR_distances_with_errors
variance_sum = 0
error_sum = 0

# normalization
for(sp in organisms_list) {
  variance_sum = variance_sum + variance_estimates[variance_estimates$species == sp & variance_estimates$tissue == tissue, "TR_var"] 
  error_sum = error_sum + error_estimate[error_estimate$species == sp & error_estimate$type == "TR" & error_estimate$tissue == tissue, "dividedVar"]
}
denumerator = (variance_sum - error_sum)/length(organisms_list)

for(sp1 in organisms_list) {
  for(sp2 in organisms_list) {
    TR_distances_with_errors_normalized[sp1,sp2] = TR_distances_with_errors[sp1,sp2]/denumerator
  }
}
TR_distances_with_errors_normalized_primed <- TR_distances_with_errors_normalized

# build tree
TR_tree_primed = nj(TR_distances_with_errors_normalized_primed)
ggtree(TR_tree_primed)+geom_tiplab()
ggtree(TR_tree_primed, layout = 'daylight', branch.length='none')

# bootstrap
TR_tree_primed
bootstrap_support1 = c()
bootstrap_support2 = c()
bootstrap_support3 = c()
bootstrap_support4 = c()
bootstrap_support5 = c()
bootstrap_support6 = c()
control_tree <- TR_tree_primed

for(i in 1:1000) {
  expression_bootstrapped = expression_combined[sample(nrow(expression_combined),replace=T),]
  
  #### estimate of within-species vriance and measurement error ####
  
  error_estimate = tibble( colname = colnames(expression_bootstrapped) ) %>%
    filter( grepl("_R", colnames(expression_bootstrapped) )) %>%
    separate( colname, c( "species", "tissue", "type", "rep" ), "(\\.|_)", remove=FALSE  ) %>%
    group_by( species, tissue, type ) %>%
    summarise( medianVar = median( genefilter::rowVars( expression_bootstrapped[,colname,drop=FALSE] ), na.rm=TRUE ) )
  error_estimate = data.frame(error_estimate)
  
  ####  variance estimates ####
  
  variance_estimates_RP_TR = tibble( colname = colnames(expression_bootstrapped) ) %>%
    filter( str_detect( colname, "median" ) ) %>%
    filter( colname != "Gene_ID" ) %>%
    separate( "colname", c( "median", "species", "tissue", "type" ),  "(\\.|_)", remove=FALSE ) %>%
    #  filter( tissue == tissue) %>%
    spread( type, colname ) %>%
    group_by_all() %>%
    summarise( TR_var = var( expression_bootstrapped[ , TR ], na.rm=TRUE ) ) %>%
    ungroup
  variance_estimates = data.frame(variance_estimates_RP_TR)
  
  #### expression divergence calculations ####
  
  TR_distances_with_errors = crossing( sp1 = organisms_list, sp2 = organisms_list ) %>%
    filter( sp1 != sp2 ) %>%
    rowwise() %>%
    group_by_all() %>%
    summarise(var_divergence =
                var(expression_bootstrapped[ , paste("median_",sp1,"_",tissue,".","TR",sep = "")] -
                      expression_bootstrapped[ , paste("median_",sp2,"_",tissue,".","TR",sep = "")], na.rm=TRUE ) -
                error_estimate[error_estimate$species == sp1 & error_estimate$tissue == tissue & error_estimate$type == "TR", "medianVar"] -
                error_estimate[error_estimate$species == sp2 & error_estimate$tissue == tissue & error_estimate$type == "TR", "medianVar"]) %>%
    ungroup %>%
    spread( sp2, var_divergence, fill=0 ) %>%
    column_to_rownames( "sp1" ) %>%
    as.matrix
  
  TR_distances_with_errors_normalized = TR_distances_with_errors
  for(sp1 in organisms_list) {
    for(sp2 in organisms_list) {
      denominator = variance_estimates[variance_estimates$species == sp1 & variance_estimates$tissue == tissue, "TR_var"] +
        variance_estimates[variance_estimates$species == sp2 & variance_estimates$tissue == tissue, "TR_var"] -
        error_estimate[error_estimate$species == sp1 & error_estimate$type == "TR" & error_estimate$tissue == tissue, "medianVar"] -
        error_estimate[error_estimate$species == sp2 & error_estimate$type == "TR" & error_estimate$tissue == tissue, "medianVar"]
      TR_distances_with_errors_normalized[sp1,sp2] = TR_distances_with_errors[sp1,sp2]/(denominator/2)
    }
  }
  
  TR_tree = nj( TR_distances_with_errors_normalized )
  
  bootstrap_support1 = c(bootstrap_support1,prop.clades(control_tree,TR_tree)[1])
  bootstrap_support2 = c(bootstrap_support2,prop.clades(control_tree,TR_tree)[2])
  bootstrap_support3 = c(bootstrap_support3,prop.clades(control_tree,TR_tree)[3])
  bootstrap_support4 = c(bootstrap_support4,prop.clades(control_tree,TR_tree)[4])
  bootstrap_support5 = c(bootstrap_support5,prop.clades(control_tree,TR_tree)[5])
  bootstrap_support6 = c(bootstrap_support6,prop.clades(control_tree,TR_tree)[6])
  
}

prop_1 <- function(x) {
  sum(x == 1, na.rm = TRUE) / length(x)
}

boots_primed <- c(
  prop_1(bootstrap_support1),
  prop_1(bootstrap_support2),
  prop_1(bootstrap_support3),
  prop_1(bootstrap_support4),
  prop_1(bootstrap_support5),
  prop_1(bootstrap_support6)
)
plot(TR_tree_primed, main = "primed")
nodelabels(boots_primed)



## set tissue: naive
# estimate of within-species vriance and measurement error ####
error_estimate = tibble( colname = colnames(expression_combined) ) %>%
  filter( ! str_detect( colname, "median" ) ) %>%
  filter( colname != "Gene_ID" ) %>%
  separate( colname, c( "species", "tissue", "type", "rep" ), "(\\.|_)", remove=FALSE  ) %>%
  group_by( species, tissue, type ) %>%
  summarise( meanVar = mean( genefilter::rowVars( expression_combined[,colname,drop=FALSE] ), na.rm=TRUE ) ) 
error_estimate = data.frame(error_estimate)

# count the sample number of each group
group_number = tibble( colname = colnames(expression_combined) ) %>%
  filter( ! str_detect( colname, "median" ) ) %>%
  filter( colname != "Gene_ID" ) %>%
  separate( colname, c( "species", "type", "rep" ), "(\\.|_)", remove=FALSE  ) %>%
  group_by( species, type ) %>% summarise(n = n())

## divide by group number
error_estimate$dividedVar = error_estimate$meanVar/group_number$n

# estimate within group variance
variance_estimates_TR = tibble( colname = colnames(expression_combined) ) %>%
  filter( str_detect( colname, "median" ) ) %>%
  filter( colname != "Gene_ID" ) %>%
  separate( "colname", c( "median", "species", "tissue", "type" ),  "(\\.|_)", remove=FALSE ) %>%
  spread( type, colname ) %>%
  group_by_all() %>%
  summarise( TR_var = var( expression_combined[ , TR ], na.rm=TRUE ) ) %>%
  ungroup
variance_estimates = data.frame(variance_estimates_TR)

tissue = "naïve"
group_number[which(group_number$type=="naïve"), 1]
organisms_list = c("human","mouse","pig","crabmacaque","rat","marmoset") #naive

#### transcriptome layer divergence calculations, corrected ####
TR_distances_with_errors = crossing( sp1 = organisms_list, sp2 = organisms_list ) %>%
  filter( sp1 != sp2 ) %>%
  rowwise() %>%
  group_by_all() %>%
  reframe(var_divergence = 
            var(expression_combined[ , paste("median_",sp1,"_",tissue,".","TR",sep = "")] - 
                  expression_combined[ , paste("median_",sp2,"_",tissue,".","TR",sep = "")], na.rm=TRUE ) - 
            error_estimate[error_estimate$species == sp1 & error_estimate$tissue == tissue & error_estimate$type == "TR", "dividedVar"] -
            error_estimate[error_estimate$species == sp2 & error_estimate$tissue == tissue & error_estimate$type == "TR", "dividedVar"]) %>%
  #ungroup %>%
  spread( sp2, var_divergence, fill=0 ) %>%
  column_to_rownames( "sp1" ) %>%
  as.matrix
TR_distances_with_errors_normalized = TR_distances_with_errors
variance_sum = 0
error_sum = 0


for(sp in organisms_list) {
  variance_sum = variance_sum + variance_estimates[variance_estimates$species == sp & variance_estimates$tissue == tissue, "TR_var"] 
  error_sum = error_sum + error_estimate[error_estimate$species == sp & error_estimate$type == "TR" & error_estimate$tissue == tissue, "dividedVar"]
}
denumerator = (variance_sum - error_sum)/length(organisms_list)

for(sp1 in organisms_list) {
  for(sp2 in organisms_list) {
    TR_distances_with_errors_normalized[sp1,sp2] = TR_distances_with_errors[sp1,sp2]/denumerator
  }
}

TR_distances_with_errors_normalized_naive <- TR_distances_with_errors_normalized
TR_tree_naive = nj(TR_distances_with_errors_normalized_naive)
ggtree(TR_tree_naive)+geom_tiplab()

# bootstrap
TR_tree_naive
bootstrap_support1 = c()
bootstrap_support2 = c()
bootstrap_support3 = c()
bootstrap_support4 = c()

control_tree <- TR_tree_naive

for(i in 1:1000) {
  expression_bootstrapped = expression_combined[sample(nrow(expression_combined),replace=T),]
  
  #### estimate of within-species vriance and measurement error ####
  
  error_estimate = tibble( colname = colnames(expression_bootstrapped) ) %>%
    filter( grepl("_R", colnames(expression_bootstrapped) )) %>%
    separate( colname, c( "species", "tissue", "type", "rep" ), "(\\.|_)", remove=FALSE  ) %>%
    group_by( species, tissue, type ) %>%
    summarise( medianVar = median( genefilter::rowVars( expression_bootstrapped[,colname,drop=FALSE] ), na.rm=TRUE ) )
  error_estimate = data.frame(error_estimate)
  
  ####  variance estimates ####
  
  variance_estimates_RP_TR = tibble( colname = colnames(expression_bootstrapped) ) %>%
    filter( str_detect( colname, "median" ) ) %>%
    filter( colname != "Gene_ID" ) %>%
    separate( "colname", c( "median", "species", "tissue", "type" ),  "(\\.|_)", remove=FALSE ) %>%
    #  filter( tissue == tissue) %>%
    spread( type, colname ) %>%
    group_by_all() %>%
    summarise( TR_var = var( expression_bootstrapped[ , TR ], na.rm=TRUE ) ) %>%
    ungroup
  variance_estimates = data.frame(variance_estimates_RP_TR)
  
  #### expression divergence calculations ####
  
  TR_distances_with_errors = crossing( sp1 = organisms_list, sp2 = organisms_list ) %>%
    filter( sp1 != sp2 ) %>%
    rowwise() %>%
    group_by_all() %>%
    summarise(var_divergence =
                var(expression_bootstrapped[ , paste("median_",sp1,"_",tissue,".","TR",sep = "")] -
                      expression_bootstrapped[ , paste("median_",sp2,"_",tissue,".","TR",sep = "")], na.rm=TRUE ) -
                error_estimate[error_estimate$species == sp1 & error_estimate$tissue == tissue & error_estimate$type == "TR", "medianVar"] -
                error_estimate[error_estimate$species == sp2 & error_estimate$tissue == tissue & error_estimate$type == "TR", "medianVar"]) %>%
    ungroup %>%
    spread( sp2, var_divergence, fill=0 ) %>%
    column_to_rownames( "sp1" ) %>%
    as.matrix
  
  TR_distances_with_errors_normalized = TR_distances_with_errors
  for(sp1 in organisms_list) {
    for(sp2 in organisms_list) {
      denominator = variance_estimates[variance_estimates$species == sp1 & variance_estimates$tissue == tissue, "TR_var"] +
        variance_estimates[variance_estimates$species == sp2 & variance_estimates$tissue == tissue, "TR_var"] -
        error_estimate[error_estimate$species == sp1 & error_estimate$type == "TR" & error_estimate$tissue == tissue, "medianVar"] -
        error_estimate[error_estimate$species == sp2 & error_estimate$type == "TR" & error_estimate$tissue == tissue, "medianVar"]
      TR_distances_with_errors_normalized[sp1,sp2] = TR_distances_with_errors[sp1,sp2]/(denominator/2)
    }
  }
  
  TR_tree = nj( TR_distances_with_errors_normalized )
  
  bootstrap_support1 = c(bootstrap_support1,prop.clades(control_tree,TR_tree)[1])
  bootstrap_support2 = c(bootstrap_support2,prop.clades(control_tree,TR_tree)[2])
  bootstrap_support3 = c(bootstrap_support3,prop.clades(control_tree,TR_tree)[3])
  bootstrap_support4 = c(bootstrap_support4,prop.clades(control_tree,TR_tree)[4])
  
}

prop_1 <- function(x) {
  sum(x == 1, na.rm = TRUE) / length(x)
}

boots_naive <- c(
  prop_1(bootstrap_support1),
  prop_1(bootstrap_support2),
  prop_1(bootstrap_support3),
  prop_1(bootstrap_support4)
)
plot(TR_tree_naive, main = "naïve")
nodelabels(boots_naive)


## set tissue: formative
# estimate of within-species vriance and measurement error ####
error_estimate = tibble( colname = colnames(expression_combined) ) %>%
  filter( ! str_detect( colname, "median" ) ) %>%
  filter( colname != "Gene_ID" ) %>%
  separate( colname, c( "species", "tissue", "type", "rep" ), "(\\.|_)", remove=FALSE  ) %>%
  group_by( species, tissue, type ) %>%
  summarise( meanVar = mean( genefilter::rowVars( expression_combined[,colname,drop=FALSE] ), na.rm=TRUE ) ) 
error_estimate = data.frame(error_estimate)

# count the sample number of each group
group_number = tibble( colname = colnames(expression_combined) ) %>%
  filter( ! str_detect( colname, "median" ) ) %>%
  filter( colname != "Gene_ID" ) %>%
  separate( colname, c( "species", "type", "rep" ), "(\\.|_)", remove=FALSE  ) %>%
  group_by( species, type ) %>% summarise(n = n())

## divide by group number
error_estimate$dividedVar = error_estimate$meanVar/group_number$n

# estimate within group variance
variance_estimates_TR = tibble( colname = colnames(expression_combined) ) %>%
  filter( str_detect( colname, "median" ) ) %>%
  filter( colname != "Gene_ID" ) %>%
  separate( "colname", c( "median", "species", "tissue", "type" ),  "(\\.|_)", remove=FALSE ) %>%
  spread( type, colname ) %>%
  group_by_all() %>%
  summarise( TR_var = var( expression_combined[ , TR ], na.rm=TRUE ) ) %>%
  ungroup
variance_estimates = data.frame(variance_estimates_TR)

tissue = "formative"
group_number[which(group_number$type=="formative"), 1]
organisms_list = c("cattle","human","horse","mouse","pig","rat") #formative

#### transcriptome layer divergence calculations, corrected ####
TR_distances_with_errors = crossing( sp1 = organisms_list, sp2 = organisms_list ) %>%
  filter( sp1 != sp2 ) %>%
  rowwise() %>%
  group_by_all() %>%
  reframe(var_divergence = 
            var(expression_combined[ , paste("median_",sp1,"_",tissue,".","TR",sep = "")] - 
                  expression_combined[ , paste("median_",sp2,"_",tissue,".","TR",sep = "")], na.rm=TRUE ) - 
            error_estimate[error_estimate$species == sp1 & error_estimate$tissue == tissue & error_estimate$type == "TR", "dividedVar"] -
            error_estimate[error_estimate$species == sp2 & error_estimate$tissue == tissue & error_estimate$type == "TR", "dividedVar"]) %>%
  #ungroup %>%
  spread( sp2, var_divergence, fill=0 ) %>%
  column_to_rownames( "sp1" ) %>%
  as.matrix

TR_distances_with_errors_normalized = TR_distances_with_errors

variance_sum = 0
error_sum = 0

for(sp in organisms_list) {
  variance_sum = variance_sum + variance_estimates[variance_estimates$species == sp & variance_estimates$tissue == tissue, "TR_var"] 
  error_sum = error_sum + error_estimate[error_estimate$species == sp & error_estimate$type == "TR" & error_estimate$tissue == tissue, "dividedVar"]
}
denumerator = (variance_sum - error_sum)/length(organisms_list)

for(sp1 in organisms_list) {
  for(sp2 in organisms_list) {
    TR_distances_with_errors_normalized[sp1,sp2] = TR_distances_with_errors[sp1,sp2]/denumerator
  }
}

TR_distances_with_errors_normalized_formative <- TR_distances_with_errors_normalized
TR_tree_formative = nj(TR_distances_with_errors_normalized_formative)
ggtree(TR_tree_formative)+geom_tiplab()

# bootstrap
TR_tree_formative
bootstrap_support1 = c()
bootstrap_support2 = c()
bootstrap_support3 = c()

control_tree <- TR_tree_formative

for(i in 1:1000) {
  expression_bootstrapped = expression_combined[sample(nrow(expression_combined),replace=T),]
  
  #### estimate of within-species vriance and measurement error ####
  
  error_estimate = tibble( colname = colnames(expression_bootstrapped) ) %>%
    filter( grepl("_R", colnames(expression_bootstrapped) )) %>%
    separate( colname, c( "species", "tissue", "type", "rep" ), "(\\.|_)", remove=FALSE  ) %>%
    group_by( species, tissue, type ) %>%
    summarise( medianVar = median( genefilter::rowVars( expression_bootstrapped[,colname,drop=FALSE] ), na.rm=TRUE ) )
  error_estimate = data.frame(error_estimate)
  
  ####  variance estimates ####
  
  variance_estimates_RP_TR = tibble( colname = colnames(expression_bootstrapped) ) %>%
    filter( str_detect( colname, "median" ) ) %>%
    filter( colname != "Gene_ID" ) %>%
    separate( "colname", c( "median", "species", "tissue", "type" ),  "(\\.|_)", remove=FALSE ) %>%
    #  filter( tissue == tissue) %>%
    spread( type, colname ) %>%
    group_by_all() %>%
    summarise( TR_var = var( expression_bootstrapped[ , TR ], na.rm=TRUE ) ) %>%
    ungroup
  variance_estimates = data.frame(variance_estimates_RP_TR)
  
  #### expression divergence calculations ####
  
  TR_distances_with_errors = crossing( sp1 = organisms_list, sp2 = organisms_list ) %>%
    filter( sp1 != sp2 ) %>%
    rowwise() %>%
    group_by_all() %>%
    summarise(var_divergence =
                var(expression_bootstrapped[ , paste("median_",sp1,"_",tissue,".","TR",sep = "")] -
                      expression_bootstrapped[ , paste("median_",sp2,"_",tissue,".","TR",sep = "")], na.rm=TRUE ) -
                error_estimate[error_estimate$species == sp1 & error_estimate$tissue == tissue & error_estimate$type == "TR", "medianVar"] -
                error_estimate[error_estimate$species == sp2 & error_estimate$tissue == tissue & error_estimate$type == "TR", "medianVar"]) %>%
    ungroup %>%
    spread( sp2, var_divergence, fill=0 ) %>%
    column_to_rownames( "sp1" ) %>%
    as.matrix
  
  TR_distances_with_errors_normalized = TR_distances_with_errors
  for(sp1 in organisms_list) {
    for(sp2 in organisms_list) {
      denominator = variance_estimates[variance_estimates$species == sp1 & variance_estimates$tissue == tissue, "TR_var"] +
        variance_estimates[variance_estimates$species == sp2 & variance_estimates$tissue == tissue, "TR_var"] -
        error_estimate[error_estimate$species == sp1 & error_estimate$type == "TR" & error_estimate$tissue == tissue, "medianVar"] -
        error_estimate[error_estimate$species == sp2 & error_estimate$type == "TR" & error_estimate$tissue == tissue, "medianVar"]
      TR_distances_with_errors_normalized[sp1,sp2] = TR_distances_with_errors[sp1,sp2]/(denominator/2)
    }
  }
  
  TR_tree = nj( TR_distances_with_errors_normalized )
  
  bootstrap_support1 = c(bootstrap_support1,prop.clades(control_tree,TR_tree)[1])
  bootstrap_support2 = c(bootstrap_support2,prop.clades(control_tree,TR_tree)[2])
  bootstrap_support3 = c(bootstrap_support3,prop.clades(control_tree,TR_tree)[3])
  
}

prop_1 <- function(x) {
  sum(x == 1, na.rm = TRUE) / length(x)
}

boots_formative <- c(
  prop_1(bootstrap_support1),
  prop_1(bootstrap_support2),
  prop_1(bootstrap_support3)
)
plot(TR_tree_formative, main = "formative")
nodelabels(boots_formative)



## set tissue: expanded
# estimate of within-species vriance and measurement error ####
error_estimate = tibble( colname = colnames(expression_combined) ) %>%
  filter( ! str_detect( colname, "median" ) ) %>%
  filter( colname != "Gene_ID" ) %>%
  separate( colname, c( "species", "tissue", "type", "rep" ), "(\\.|_)", remove=FALSE  ) %>%
  group_by( species, tissue, type ) %>%
  summarise( meanVar = mean( genefilter::rowVars( expression_combined[,colname,drop=FALSE] ), na.rm=TRUE ) ) 
error_estimate = data.frame(error_estimate)

# count the sample number of each group
group_number = tibble( colname = colnames(expression_combined) ) %>%
  filter( ! str_detect( colname, "median" ) ) %>%
  filter( colname != "Gene_ID" ) %>%
  separate( colname, c( "species", "type", "rep" ), "(\\.|_)", remove=FALSE  ) %>%
  group_by( species, type ) %>% summarise(n = n())

## divide by group number
error_estimate$dividedVar = error_estimate$meanVar/group_number$n

# estimate within group variance
variance_estimates_TR = tibble( colname = colnames(expression_combined) ) %>%
  filter( str_detect( colname, "median" ) ) %>%
  filter( colname != "Gene_ID" ) %>%
  separate( "colname", c( "median", "species", "tissue", "type" ),  "(\\.|_)", remove=FALSE ) %>%
  spread( type, colname ) %>%
  group_by_all() %>%
  summarise( TR_var = var( expression_combined[ , TR ], na.rm=TRUE ) ) %>%
  ungroup
variance_estimates = data.frame(variance_estimates_TR)

tissue = "expanded"
group_number[which(group_number$type=="expanded"), 1]
organisms_list = c("cattle","mouse","human","pig","marmoset") #expanded

#### transcriptome layer divergence calculations, corrected ####
TR_distances_with_errors = crossing( sp1 = organisms_list, sp2 = organisms_list ) %>%
  filter( sp1 != sp2 ) %>%
  rowwise() %>%
  group_by_all() %>%
  reframe(var_divergence = 
            var(expression_combined[ , paste("median_",sp1,"_",tissue,".","TR",sep = "")] - 
                  expression_combined[ , paste("median_",sp2,"_",tissue,".","TR",sep = "")], na.rm=TRUE ) - 
            error_estimate[error_estimate$species == sp1 & error_estimate$tissue == tissue & error_estimate$type == "TR", "dividedVar"] -
            error_estimate[error_estimate$species == sp2 & error_estimate$tissue == tissue & error_estimate$type == "TR", "dividedVar"]) %>%
  #ungroup %>%
  spread( sp2, var_divergence, fill=0 ) %>%
  column_to_rownames( "sp1" ) %>%
  as.matrix

TR_distances_with_errors_normalized = TR_distances_with_errors

variance_sum = 0
error_sum = 0

for(sp in organisms_list) {
  variance_sum = variance_sum + variance_estimates[variance_estimates$species == sp & variance_estimates$tissue == tissue, "TR_var"] 
  error_sum = error_sum + error_estimate[error_estimate$species == sp & error_estimate$type == "TR" & error_estimate$tissue == tissue, "dividedVar"]
}
denumerator = (variance_sum - error_sum)/length(organisms_list)

for(sp1 in organisms_list) {
  for(sp2 in organisms_list) {
    TR_distances_with_errors_normalized[sp1,sp2] = TR_distances_with_errors[sp1,sp2]/denumerator
  }
}

TR_distances_with_errors_normalized_expanded <- TR_distances_with_errors_normalized
TR_tree_expanded = nj(TR_distances_with_errors_normalized_expanded)
ggtree(TR_tree_expanded)+geom_tiplab()

# bootstrap
TR_tree_expanded
bootstrap_support1 = c()
bootstrap_support2 = c()
bootstrap_support3 = c()

control_tree <- TR_tree_expanded

for(i in 1:100) {
  expression_bootstrapped = expression_combined[sample(nrow(expression_combined),replace=T),]
  
  #### estimate of within-species vriance and measurement error ####
  
  error_estimate = tibble( colname = colnames(expression_bootstrapped) ) %>%
    filter( grepl("_R", colnames(expression_bootstrapped) )) %>%
    separate( colname, c( "species", "tissue", "type", "rep" ), "(\\.|_)", remove=FALSE  ) %>%
    group_by( species, tissue, type ) %>%
    summarise( medianVar = median( genefilter::rowVars( expression_bootstrapped[,colname,drop=FALSE] ), na.rm=TRUE ) )
  error_estimate = data.frame(error_estimate)
  
  ####  variance estimates ####
  
  variance_estimates_RP_TR = tibble( colname = colnames(expression_bootstrapped) ) %>%
    filter( str_detect( colname, "median" ) ) %>%
    filter( colname != "Gene_ID" ) %>%
    separate( "colname", c( "median", "species", "tissue", "type" ),  "(\\.|_)", remove=FALSE ) %>%
    #  filter( tissue == tissue) %>%
    spread( type, colname ) %>%
    group_by_all() %>%
    summarise( TR_var = var( expression_bootstrapped[ , TR ], na.rm=TRUE ) ) %>%
    ungroup
  variance_estimates = data.frame(variance_estimates_RP_TR)
  
  #### expression divergence calculations ####
  
  TR_distances_with_errors = crossing( sp1 = organisms_list, sp2 = organisms_list ) %>%
    filter( sp1 != sp2 ) %>%
    rowwise() %>%
    group_by_all() %>%
    summarise(var_divergence =
                var(expression_bootstrapped[ , paste("median_",sp1,"_",tissue,".","TR",sep = "")] -
                      expression_bootstrapped[ , paste("median_",sp2,"_",tissue,".","TR",sep = "")], na.rm=TRUE ) -
                error_estimate[error_estimate$species == sp1 & error_estimate$tissue == tissue & error_estimate$type == "TR", "medianVar"] -
                error_estimate[error_estimate$species == sp2 & error_estimate$tissue == tissue & error_estimate$type == "TR", "medianVar"]) %>%
    ungroup %>%
    spread( sp2, var_divergence, fill=0 ) %>%
    column_to_rownames( "sp1" ) %>%
    as.matrix
  
  TR_distances_with_errors_normalized = TR_distances_with_errors
  for(sp1 in organisms_list) {
    for(sp2 in organisms_list) {
      denominator = variance_estimates[variance_estimates$species == sp1 & variance_estimates$tissue == tissue, "TR_var"] +
        variance_estimates[variance_estimates$species == sp2 & variance_estimates$tissue == tissue, "TR_var"] -
        error_estimate[error_estimate$species == sp1 & error_estimate$type == "TR" & error_estimate$tissue == tissue, "medianVar"] -
        error_estimate[error_estimate$species == sp2 & error_estimate$type == "TR" & error_estimate$tissue == tissue, "medianVar"]
      TR_distances_with_errors_normalized[sp1,sp2] = TR_distances_with_errors[sp1,sp2]/(denominator/2)
    }
  }
  
  TR_tree = nj( TR_distances_with_errors_normalized )
  
  bootstrap_support1 = c(bootstrap_support1,prop.clades(control_tree,TR_tree)[1])
  bootstrap_support2 = c(bootstrap_support2,prop.clades(control_tree,TR_tree)[2])
  bootstrap_support3 = c(bootstrap_support3,prop.clades(control_tree,TR_tree)[3])
  
}

prop_1 <- function(x) {
  sum(x == 1, na.rm = TRUE) / length(x)
}

boots_expanded <- c(
  prop_1(bootstrap_support1),
  prop_1(bootstrap_support2),
  prop_1(bootstrap_support3)
)
plot(TR_tree_expanded, main = "expanded")
nodelabels(boots_expanded)

p <- ggtree(TR_tree_expanded, layout = "circular") + 
  geom_tiplab(size = 2.5) +   # Adjust size as needed
  geom_text2(aes(subset = !isTip, label = boots_expanded), hjust = -0.5, size = 3) + 
  ggtitle("Expanded")
# Print the plot
print(p)


# draw the plot
library(ggplot2)
library(patchwork)


create_tree_plot <- function(tree, boots, title, tree_color, label_color) {
  tree_plot <- ggtree(tree, layout = "circular", color = tree_color) + 
    geom_tiplab(size = 5.5, color = "black", fontface = "italic")
  
  internal_nodes <- tree_plot$data %>% filter(!isTip) %>% arrange(node)
  
  tree_plot <- tree_plot + 
    geom_text2(aes(x = x, y = y, label = ifelse(isTip, NA, boots)), 
               data = internal_nodes, 
               hjust = -0.5, size = 4, color = label_color) + 
    ggtitle(title) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank()
    )
  
  return(tree_plot)
}

# Define colors for each tree
colors <- list(
  expanded = c(tree_color = "black", label_color = "black"),  # Light blue and darker blue
  naive = c(tree_color = "black", label_color = "black"),     # Light orange and darker orange
  formative = c(tree_color = "black", label_color = "black"), # Light green and darker green
  primed = c(tree_color = "black", label_color = "black")     # Light pink and darker pink
)

# Create individual tree plots with different colors
expanded_plot <- create_tree_plot(TR_tree_expanded, boots_expanded, "Expanded", colors$expanded["tree_color"], colors$expanded["label_color"])
naive_plot <- create_tree_plot(TR_tree_naive, boots_naive, "Naïve", colors$naive["tree_color"], colors$naive["label_color"])
formative_plot <- create_tree_plot(TR_tree_formative, boots_formative, "Formative", colors$formative["tree_color"], colors$formative["label_color"])
primed_plot <- create_tree_plot(TR_tree_primed, boots_primed, "Primed", colors$primed["tree_color"], colors$primed["label_color"])


# Combine the plots into a single figure
combined_plot <- (expanded_plot | naive_plot) / (formative_plot | primed_plot)

# Print the combined plot
print(combined_plot)

# Save the combined plot to a file
ggsave("combined_phylogenetic_trees.png", plot = combined_plot, width = 8, height = 10, dpi = 300)
ggsave("combined_phylogenetic_trees.svg", plot = combined_plot, width = 8, height = 8, dpi = 300)
# Print the combined plot
print(combined_plot)


#save.image("07.11_tree.rds")
# load("07.11_tree.rds")
