####################################################################
######Test GBL place                                  ##############
######Limma-voom new annotation with roughskin sculpin##############
#########Jan 28, 2026                               ##############
####################################################################
setwd("/Users/Sherry/Desktop/RNA_bioinformatics")

### VARIANCE PARTITIONING FOR HABITAT ASSIGNMENT COMPARISON ###################
# Load required libraries
library(edgeR)
library(variancePartition)
library(BiocParallel)
library(tidyverse)
library(tibble)

# Import metadata
metadata <- read.delim("PS_gill_RNA_metafile3.txt", sep = "\t") %>%
  filter(!Sample %in% c("CR4SW", "LCR1SW")) %>%
  mutate(
    Region = factor(Region, levels = c("coastR", "coastL", "inlandL")),
    Treatment = factor(Treatment, levels = c("AW", "SW")),
    Population = factor(Population)
  )

head(metadata)
tail(metadata)
# === COUNTS DATA IMPORT ===
# Import and filter counts
TPM_quantmerge <- read.delim("DE_Results_limma_voom_2025/real_annotation/salmon_tx_TPM.txt", row.names = 1)
quantmerge <- read.delim("DE_Results_limma_voom_2025/real_annotation/salmon_tx_NumReads.txt", row.names = 1) %>%
  .[rowMeans(TPM_quantmerge) > 1.0, metadata$Sample] %>%
  as.matrix() %>% `mode<-`("integer") %>% replace(is.na(.), 0)

# Get common samples after filtering
common_samples <- intersect(colnames(quantmerge), metadata$Sample)

# Subset and reorder both objects
quantmerge <- quantmerge[, common_samples, drop = FALSE]
metadata <- metadata %>% 
  filter(Sample %in% common_samples) %>%
  arrange(match(Sample, common_samples)) %>%
  column_to_rownames("Sample")

# Verify alignment
stopifnot(identical(colnames(quantmerge), rownames(metadata)))

### 2. NORMALIZATION AND FILTERING ############################################

# Create DGEList and normalize
dge_new_arrange <- DGEList(counts = quantmerge)
dge_new_arrange <- calcNormFactors(dge_new_arrange)

# Check BEFORE filtering
cat("=== BEFORE FILTERING ===\n")
cat("Total genes:", nrow(dge_new_arrange), "\n")
cat("Genes with zero counts across all samples:", sum(rowSums(dge_new_arrange$counts) == 0), "\n")
cat("Total samples:", ncol(dge_new_arrange), "\n")

# Create filter criteria - keep genes expressed in at least 10% of samples
keep_new_arrange <- rowSums(dge_new_arrange$counts > 0) >= 0.1 * ncol(dge_new_arrange)
cat("Genes passing expression filter (≥10% of samples):", sum(keep_new_arrange), "\n")
cat("Genes failing expression filter:", sum(!keep_new_arrange), "\n")

# APPLY FILTER
dge_new_arrange_filtered <- dge_new_arrange[keep_new_arrange, , keep.lib.sizes = FALSE]

# Check AFTER filtering
cat("\n=== AFTER FILTERING ===\n")
cat("Genes remaining:", nrow(dge_new_arrange_filtered), "\n")
cat("Percentage of original genes retained:", round(nrow(dge_new_arrange_filtered)/nrow(dge_new_arrange)*100, 2), "%\n")

# Final alignment verification
cat("\n=== SAMPLE ALIGNMENT CHECK ===\n")
cat("Count data samples:", paste(colnames(dge_new_arrange_filtered), collapse = ", "), "\n")
cat("Metadata samples:", paste(rownames(metadata), collapse = ", "), "\n")
stopifnot(identical(colnames(dge_new_arrange_filtered), rownames(metadata)))

### 3. VOOM TRANSFORMATION ####################################################

cat("\n=== PERFORMING VOOM TRANSFORMATION ===\n")

# Define the model formula for voom
formula_voom <- ~ Region + Treatment + (1|Population)

# Perform voom transformation with dream weights
vobj_new_arrange <- voomWithDreamWeights(
  counts = dge_new_arrange_filtered, 
  formula = formula_voom, 
  data = metadata,
  BPPARAM = SerialParam()
)

cat("Voom transformation completed successfully\n")

### 4. VARIANCE PARTITIONING COMPARISON #######################################

### 4.1 CREATE ALTERNATIVE METADATA WITH DIFFERENT ASSIGNMENTS ###

# Original metadata (assuming current assignment)
metadata_original <- metadata

# GBL is the ambiguous population - test coastR vs coastL assignment
ambiguous_pop <- "GBL"
option1_region <- "coastR"  # Alternative assignment
option2_region <- "coastL"  # Current assignment

# Create two alternative metadata versions
metadata_option1 <- metadata_original
metadata_option2 <- metadata_original

# Modify the region assignment for GBL population
metadata_option1$Region[metadata_option1$Population == ambiguous_pop] <- option1_region
metadata_option2$Region[metadata_option2$Population == ambiguous_pop] <- option2_region

# Verify the changes
cat("=== ORIGINAL ASSIGNMENT (GBL in coastL) ===\n")
print(table(metadata_original$Population, metadata_original$Region))

cat("\n=== OPTION 1: GBL assigned to coastR ===\n")
print(table(metadata_option1$Population, metadata_option1$Region))

cat("\n=== OPTION 2: GBL assigned to coastL (current) ===\n")
print(table(metadata_option2$Population, metadata_option2$Region))

### 4.2 DEFINE THE FULL MODEL FOR VARIANCE PARTITIONING ###
formula_vp <- ~ (1|Region) + (1|Treatment) + (1|Region:Treatment) + (1|Population)
cat("Using variance partitioning model:", deparse(formula_vp), "\n")

### 4.3 VARIANCE PARTITIONING FOR EACH ASSIGNMENT ###
# For computational efficiency, use a subset of genes
set.seed(123)
test_genes <- sample(1:nrow(vobj_new_arrange), 2000)  # 2000 random genes

cat("\n=== RUNNING VARIANCE PARTITIONING FOR OPTION 1 (GBL in coastR) ===\n")
var_part_option1 <- fitExtractVarPartModel(
  vobj_new_arrange[test_genes, ], 
  formula_vp, 
  metadata_option1,
  BPPARAM = SerialParam()
)

##plotting volin graph for variance contribution####
library(tidyr)
library(dplyr)

vp_long <- var_part_option1 %>%
  as.data.frame() %>%
  mutate(gene = rownames(.)) %>%
  pivot_longer(
    cols = -gene,
    names_to = "factor",
    values_to = "variance"
  )
vp_long$factor <- factor(
  vp_long$factor,
  levels = c("Residuals", "Region:Treatment", "Population","Treatment","Region")
)
library(ggplot2)

ggplot(vp_long, aes(x = factor, y = variance, fill = factor)) +
  geom_violin(color = "black", trim = TRUE) +
  stat_summary(fun = median, geom = "crossbar", width = 0.15, color = "black") +
  theme_bw() +
  coord_flip() +
  ylab("Variance explained (%)") +
  xlab("") +
  scale_fill_brewer(palette = "Set2") +
  theme(legend.position = "none")
#######################################################
cat("\n=== RUNNING VARIANCE PARTITIONING FOR OPTION 2 (GBL in coastL) ===\n")

var_part_option2 <- fitExtractVarPartModel(
  vobj_new_arrange[test_genes, ], 
  formula_vp, 
  metadata_option2,
  BPPARAM = SerialParam()
)

### 4.4 COMPARE MODEL PERFORMANCE ###

# Calculate average variance explained by each component
avg_var_option1 <- colMeans(var_part_option1)
avg_var_option2 <- colMeans(var_part_option2)

# Create comparison dataframe
comparison_df <- data.frame(
  VarianceComponent = names(avg_var_option1),
  GBL_in_coastR = avg_var_option1,
  GBL_in_coastL = avg_var_option2,
  Difference = avg_var_option2 - avg_var_option1  # coastL - coastR
)

cat("\n=== VARIANCE PARTITIONING COMPARISON ===\n")
print(comparison_df)

# Key metrics for comparison
cat("\n=== KEY COMPARISON METRICS ===\n")
cat("Total explained variance (1 - Residuals):\n")
cat("GBL in coastR:", 1 - avg_var_option1["Residuals"], "\n")
cat("GBL in coastL:", 1 - avg_var_option2["Residuals"], "\n")

cat("\nRegion variance (main effect):\n")
cat("GBL in coastR:", avg_var_option1["Region"], "\n")
cat("GBL in coastL:", avg_var_option2["Region"], "\n")

cat("\nRegion:Treatment interaction variance:\n")
cat("GBL in coastR:", avg_var_option1["Region:Treatment"], "\n")
cat("GBL in coastL:", avg_var_option2["Region:Treatment"], "\n")

cat("\nPopulation variance:\n")
cat("GBL in coastR:", avg_var_option1["Population"], "\n")
cat("GBL in coastL:", avg_var_option2["Population"], "\n")

### 4.5 STATISTICAL COMPARISON ###

# Compare if the differences are significant using paired t-tests
cat("\n=== STATISTICAL COMPARISONS (paired t-tests) ===\n")

for(component in names(avg_var_option1)) {
  if(component != "Residuals") {
    test_result <- t.test(var_part_option1[[component]], 
                          var_part_option2[[component]], 
                          paired = TRUE)
    cat(component, ": p-value =", format.pval(test_result$p.value, digits = 3), "\n")
  }
}

### 4.6 VISUALIZE THE COMPARISON ###

# Prepare data for plotting
plot_data <- comparison_df %>%
  filter(VarianceComponent != "Residuals") %>%
  pivot_longer(cols = c(GBL_in_coastR, GBL_in_coastL), 
               names_to = "Assignment", 
               values_to = "Variance")

# Bar plot comparison
ggplot(plot_data, aes(x = VarianceComponent, y = Variance, fill = Assignment)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Variance Partitioning: GBL Habitat Assignment Comparison",
       subtitle = "Testing whether GBL fits better in coastR or coastL",
       x = "Variance Component", y = "Average Proportion of Variance") +
  scale_fill_manual(values = c("GBL_in_coastR" = "#E41A1C", "GBL_in_coastL" = "#377EB8")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

### 4.7 INTERPRETATION AND DECISION CRITERIA ###

cat("\n=== INTERPRETATION GUIDE ===\n")
cat("Which assignment is better for GBL population?\n\n")

cat("1. BETTER EXPLAINED VARIANCE:\n")
if((1 - avg_var_option1["Residuals"]) > (1 - avg_var_option2["Residuals"])) {
  cat("   ✓ GBL in coastR explains MORE total variance\n")
} else {
  cat("   ✓ GBL in coastL explains MORE total variance\n")
}

cat("\n2. REGION EFFECT:\n")
if(avg_var_option1["Region"] > avg_var_option2["Region"]) {
  cat("   ✓ GBL in coastR gives STRONGER Region effects\n")
} else {
  cat("   ✓ GBL in coastL gives STRONGER Region effects\n")
}

cat("\n3. POPULATION EFFECT:\n")
if(avg_var_option1["Population"] < avg_var_option2["Population"]) {
  cat("   ✓ GBL in coastR explains population differences BETTER via habitat\n")
} else {
  cat("   ✓ GBL in coastL explains population differences BETTER via habitat\n")
}

cat("\n4. INTERACTION EFFECTS:\n")
if(avg_var_option1["Region:Treatment"] > avg_var_option2["Region:Treatment"]) {
  cat("   ✓ GBL in coastR gives STRONGER habitat-specific treatment responses\n")
} else {
  cat("   ✓ GBL in coastL gives STRONGER habitat-specific treatment responses\n")
}

### 5. SAVE RESULTS ##########################################################

saveRDS(list(
  option1 = var_part_option1,
  option2 = var_part_option2,
  comparison = comparison_df,
  metadata_option1 = metadata_option1,
  metadata_option2 = metadata_option2,
  vobj = vobj_new_arrange
), file = "GBL_habitat_assignment_comparison.rds")

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Results saved to: GBL_habitat_assignment_comparison.rds\n")

##generate and save the plots###########################################
# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(variancePartition)

# Create violin plot for variance partition results
create_variance_violin_plot <- function(var_part_results, 
                                        scenario_name = "Variance Partition",
                                        remove_residuals = TRUE) {
  
  # Convert to long format
  plot_data <- var_part_results %>%
    as.data.frame() %>%
    pivot_longer(cols = everything(), 
                 names_to = "Component", 
                 values_to = "Variance") %>%
    mutate(Variance = Variance * 100)  # Convert to percentage
  
  # Remove residuals if requested
    if (remove_residuals) {
    plot_data <- plot_data %>% filter(Component != "Residuals")
  }
  
  # Calculate mean variance for each component
  mean_variance <- plot_data %>%
    group_by(Component) %>%
    summarise(Mean = mean(Variance)) %>%
    arrange(Mean)
  
  # Set factor levels based on mean variance (for ordering)
  plot_data$Component <- factor(plot_data$Component, 
                                levels = mean_variance$Component)
  
  # Create violin plot
  p <- ggplot(plot_data, aes(x = Component, y = Variance)) +
    geom_violin(fill = "lightblue", alpha = 0.7, scale = "width", width = 0.8) +
    geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA) +
    geom_point(data = mean_variance, aes(x = Component, y = Mean), 
               color = "red", size = 2, shape = 18) +
    geom_text(data = mean_variance, 
              aes(x = Component, y = Mean, 
                  label = sprintf("%.1f%%", Mean)),
              vjust = -1, size = 3.5, fontface = "bold", color = "darkred") +
    labs(
      title = "Distribution of Variance Explained Across Genes",
      subtitle = scenario_name,
      x = "Variance Component",
      y = "Variance Explained (%)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 10)),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank()
    ) +
    scale_y_continuous(limits = c(0, NA))  # Start from 0%
  
  return(p)
}
# 1. Individual violin plots
violin_plot1 <- create_variance_violin_plot(var_part_option1, "GBL in coastR")
violin_plot2 <- create_variance_violin_plot(var_part_option2, "GBL in coastL")

print(violin_plot1)
print(violin_plot2)
# 1. Individual violin plots
violin_plot1 <- create_variance_violin_plot(var_part_option1, "GBL in coastR")
violin_plot2 <- create_variance_violin_plot(var_part_option2, "GBL in coastL")

print(violin_plot1)
print(violin_plot2)

ggsave("DE_Results_limma_voom_2025/real_annotation/variance_violin_coastR.pdf", violin_plot1, 
       width = 10, height = 6, device = "pdf")
ggsave("DE_Results_limma_voom_2025/real_annotation/variance_violin_coastL.pdf", violin_plot2, 
       width = 10, height = 6, device = "pdf")

# Create comparison violin plot for two scenarios
create_comparison_violin_plot <- function(var_part_option1, var_part_option2,
                                          name1 = "GBL in coastR", 
                                          name2 = "GBL in coastL",
                                          remove_residuals = TRUE) {
  
  # Process first scenario
  plot_data1 <- var_part_option1 %>%
    as.data.frame() %>%
    pivot_longer(cols = everything(), 
                 names_to = "Component", 
                 values_to = "Variance") %>%
    mutate(Variance = Variance * 100,
           Scenario = name1)
  
  # Process second scenario
  plot_data2 <- var_part_option2 %>%
    as.data.frame() %>%
    pivot_longer(cols = everything(), 
                 names_to = "Component", 
                 values_to = "Variance") %>%
    mutate(Variance = Variance * 100,
           Scenario = name2)
  
  # Combine data
  plot_data_combined <- bind_rows(plot_data1, plot_data2)
  
  # Remove residuals if requested
  if (remove_residuals) {
    plot_data_combined <- plot_data_combined %>% filter(Component != "Residuals")
  }
  
  # Calculate mean variance for ordering
  mean_variance <- plot_data_combined %>%
    group_by(Component, Scenario) %>%
    summarise(Mean = mean(Variance)) %>%
    ungroup()
  
  # Set factor levels based on overall mean variance
  component_order <- mean_variance %>%
    group_by(Component) %>%
    summarise(OverallMean = mean(Mean)) %>%
    arrange(OverallMean) %>%
    pull(Component)
  
  plot_data_combined$Component <- factor(plot_data_combined$Component, 
                                         levels = component_order)
  
  # Create comparison violin plot
  p <- ggplot(plot_data_combined, aes(x = Component, y = Variance, fill = Scenario)) +
    geom_violin(alpha = 0.7, scale = "width", width = 0.8, 
                position = position_dodge(0.8)) +
    geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA,
                 position = position_dodge(0.8)) +
    stat_summary(fun = mean, geom = "point", aes(group = Scenario),
                 position = position_dodge(0.8), 
                 color = "black", size = 2, shape = 18) +
    labs(
      title = "Variance Partition Comparison",
      subtitle = "Distribution Across Genes",
      x = "Variance Component",
      y = "Variance Explained (%)",
      fill = "GBL Assignment"
    ) +
    scale_fill_manual(values = c("GBL in coastR" = "#E41A1C", 
                                 "GBL in coastL" = "#377EB8")) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5, margin = margin(b = 10)),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "top",
      legend.title = element_text(face = "bold"),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank()
    ) +
    scale_y_continuous(limits = c(0, NA))
  
  return(p)
}
# 2. Combined comparison violin plot

comparison_violin <- create_comparison_violin_plot(var_part_option1, var_part_option2)
print(comparison_violin)
ggsave("DE_Results_limma_voom_2025/real_annotation/variance_violin_comparison.pdf", comparison_violin, 
       width = 12, height = 7, device = "pdf")
##################################################################################################################
##This PCA is with the model which correcting for populations as random effect, this is not perfect,ideally should be done with raw transformed data
#For differential gene expression analysis#########
# Load required packages
library(limma)
library(edgeR)
library(variancePartition)
library(dplyr)
library(ggplot2)
library(gplots)
library(VennDiagram)
library(scales)
library(BiocParallel)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
setwd("/Users/Sherry/Desktop/RNA_bioinformatics")

### 1. DATA LOADING AND PREPROCESSING ########################################

# Import metadata
metadata <- read.delim("PS_gill_RNA_metafile3.txt", sep = "\t") %>%
  filter(!Sample %in% c("CR4SW", "LCR1SW")) %>%
  mutate(
    Region = factor(Region, levels = c("coastR", "coastL", "inlandL")),
    Treatment = factor(Treatment, levels = c("AW", "SW")),
    Population = factor(Population)
  )

# Import and filter counts
TPM_quantmerge <- read.delim("DE_Results_limma_voom_2025/real_annotation/salmon_tx_TPM.txt", row.names = 1)
quantmerge <- read.delim("DE_Results_limma_voom_2025/real_annotation/salmon_tx_NumReads.txt", row.names = 1) %>%
  .[rowMeans(TPM_quantmerge) > 1.0, metadata$Sample] %>%
  as.matrix() %>% `mode<-`("integer") %>% replace(is.na(.), 0)

# Get common samples after filtering
common_samples <- intersect(colnames(quantmerge), metadata$Sample)

# Subset and reorder both objects
quantmerge <- quantmerge[, common_samples, drop = FALSE]
metadata <- metadata %>% 
  filter(Sample %in% common_samples) %>%
  arrange(match(Sample, common_samples)) %>%
  column_to_rownames("Sample")

# Verify alignment
stopifnot(identical(colnames(quantmerge), rownames(metadata)))

### 2. CREATE ALTERNATIVE METADATA FOR BOTH ASSIGNMENTS ######################
# Original metadata (GBL in coastL - current assignment)
metadata_coastL <- metadata

# Alternative metadata (GBL in coastR - test assignment)
metadata_coastR <- metadata
metadata_coastR$Region[metadata_coastR$Population == "GBL"] <- "coastR"

# Verify the assignments
cat("=== GBL in coastL (Current) ===\n")
print(table(metadata_coastL$Population, metadata_coastL$Region))

cat("\n=== GBL in coastR (Alternative) ===\n")
print(table(metadata_coastR$Population, metadata_coastR$Region))

### 3. NORMALIZATION AND FILTERING ###########################################

process_counts <- function(counts_matrix, metadata_df) {
  # Create DGEList and normalize
  dge <- DGEList(counts = counts_matrix)
  dge <- calcNormFactors(dge)
  
  # Filter genes expressed in at least 10% of samples
  keep <- rowSums(dge$counts > 0) >= 0.1 * ncol(dge)
  dge_filtered <- dge[keep, , keep.lib.sizes = FALSE]
  
  # Verify alignment
  stopifnot(identical(colnames(dge_filtered), rownames(metadata_df)))
  
  return(dge_filtered)
}

# Process counts for both assignments
dge_coastL <- process_counts(quantmerge, metadata_coastL)
dge_coastR <- process_counts(quantmerge, metadata_coastR)

cat("Genes after filtering - GBL in coastL:", nrow(dge_coastL), "\n")
cat("Genes after filtering - GBL in coastR:", nrow(dge_coastR), "\n")

### 4. VOOM TRANSFORMATION FOR BOTH ASSIGNMENTS ##############################

perform_voom <- function(dge_obj, metadata_df, assignment_name) {
  cat("\n=== VOOM TRANSFORMATION - GBL in", assignment_name, "===\n")
  
  # Define model with interaction term
  formula <- ~ Region + Treatment + Region:Treatment + (1|Population)
  
  vobj <- voomWithDreamWeights(
    counts = dge_obj,
    formula = formula,
    data = metadata_df,
    BPPARAM = SerialParam()
  )
  
  cat("Voom transformation completed for", assignment_name, "\n")
  return(vobj)
}

vobj_coastL <- perform_voom(dge_coastL, metadata_coastL, "coastL")
vobj_coastR <- perform_voom(dge_coastR, metadata_coastR, "coastR")

pca_coastL <- vobj_coastL$E
pca_coastR <- vobj_coastR$E
write.csv(pca_coastL, "DE_Results_limma_voom_2025/real_annotation/pca_coastL.csv", row.names = TRUE)
write.csv(pca_coastR, "DE_Results_limma_voom_2025/real_annotation/pca_coastR.csv", row.names = TRUE)

##background genes#####
background_genes_coastL<- rownames(dge_coastL)
write.table(background_genes_coastL, 
            "DE_Results_limma_voom_2025/real_annotation/background_gene_set_coastL.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
background_genes_coastR<- rownames(dge_coastR)
write.table(background_genes_coastR, 
            "DE_Results_limma_voom_2025/real_annotation/background_gene_set_coastR.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

####PCA analysis and graph########################################################################################

pca_result_coastR <- prcomp(t(pca_coastR), scale. = TRUE, center = TRUE)

library(ggplot2)
library(stringr)
library(dplyr)
library(ggrepel)
# Calculate PCA variance percentages
pca_var_perc_R <- round(summary(pca_result_coastR)$importance[2, 1:2] * 100, 1)

# Calculate centroids
centroids <- as.data.frame(pca_result_coastR$x) %>%
  mutate(Region = metadata_coastR$Region) %>%
  group_by(Region) %>%
  summarize(PC1 = mean(PC1), PC2 = mean(PC2))

# Create the plot with publication-quality settings
pca_plot_coastR <- ggplot(as.data.frame(pca_result_coastR$x), aes(x = PC1, y = PC2)) +
  geom_point(
    aes(color = metadata_coastR$Region, shape = metadata_coastR$Treatment), 
    size = 3, 
    position = position_jitter(width = 0.1, height = 0.1),
    alpha = 0.8
  ) +
  geom_point(
    data = centroids,
    aes(x = PC1, y = PC2, color = Region),
    size = 5, shape = 4, stroke = 1.2, show.legend = FALSE
  ) +
  stat_ellipse(
    aes(group = metadata_coastR$Region, color = metadata_coastR$Region),
    level = 0.95, 
    linetype = "dashed",
    linewidth = 0.6
  ) +
  scale_color_manual(
    values = c("coastR" = "grey20",    # Dark grey
               "coastL" = "grey50",    # Medium grey  
               "inlandL" = "grey80"),  # Light grey
    labels = function(x) str_replace(x, "Region", "")
  ) +
  scale_shape_manual(values = c("AW" = 16, "SW" = 17)) +
  labs(
    x = paste0("PC1 (", pca_var_perc_R[1], "%)"),
    y = paste0("PC2 (", pca_var_perc_R[2], "%)"),
    color = "Region",
    shape = "Treatment"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    legend.background = element_blank(),
    legend.key = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    aspect.ratio = 1
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 3, alpha = 1)),
    shape = guide_legend(override.aes = list(size = 3, alpha = 1))
  )

ggsave("DE_Results_limma_voom_2025/real_annotation/PCA_plot_coastR.pdf", 
       plot = pca_plot_coastR,
       width = 8, height = 6, device = "pdf")

##########PCA with sample names to be changed##################################################################
metadata_coastL$Sample <- rownames(metadata_coastL)
metadata_coastR$Sample <- rownames(metadata_coastR)
# --- 1. Build PCA data frame ---
pca_result_coastL <- prcomp(t(pca_coastL), scale. = TRUE, center = TRUE)
pca_coastL_df <- as.data.frame(pca_result_coastL$x)

pca_coastL_df$Sample    <- metadata_coastL$Sample
pca_coastL_df$Region    <- metadata_coastL$Region
pca_coastL_df$Treatment <- metadata_coastL$Treatment

# --- 2. Variance explained ---
pca_var_perc_L <- round(summary(pca_result_coastL)$importance[2, 1:2] * 100, 1)

# --- 3. Centroids ---
centroids <- pca_coastL_df %>%
  group_by(Region) %>%
  summarize(PC1 = mean(PC1), PC2 = mean(PC2))

# --- 4. PCA plot with labels ---
pca_plot_L <- ggplot(pca_coastL_df, aes(x = PC1, y = PC2)) +
  geom_point(
    aes(color = Region, shape = Treatment),
    size = 3,
    position = position_jitter(width = 0.1, height = 0.1),
    alpha = 0.8
  ) +
  geom_text_repel(
    aes(label = Sample),
    size = 3,
    box.padding = 0.4,
    point.padding = 0.2,
    segment.size = 0.2,
    max.overlaps = Inf
  ) +
  geom_point(
    data = centroids,
    aes(x = PC1, y = PC2, color = Region),
    size = 5, shape = 4, stroke = 1.2, show.legend = FALSE
  ) +
  stat_ellipse(
    aes(group = Region, color = Region),
    level = 0.95,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  scale_color_manual(
    values = c("coastR" = "grey20",
               "coastL" = "grey50",
               "inlandL" = "grey80")
  ) +
  scale_shape_manual(values = c("AW" = 16, "SW" = 17)) +
  labs(
    x = paste0("PC1 (", pca_var_perc_L[1], "%)"),
    y = paste0("PC2 (", pca_var_perc_L[2], "%)"),
    color = "Region",
    shape = "Treatment"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    legend.background = element_blank(),
    legend.key = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    aspect.ratio = 1
  )

# Save
ggsave(
  "DE_Results_limma_voom_2025/real_annotation/PCA_plot_coastL_with_labels.pdf",
  plot = pca_plot_L,
  width = 8, height = 6,
  device = "pdf"
)

##GBL in coastR###
metadata_coastR$Sample <- rownames(metadata_coastR)
# --- 1. Build PCA data frame ---
pca_result_coastR <- prcomp(t(pca_coastR), scale. = TRUE, center = TRUE)
pca_coastR_df <- as.data.frame(pca_result_coastR$x)

pca_coastR_df$Sample    <- metadata_coastR$Sample
pca_coastR_df$Region    <- metadata_coastR$Region
pca_coastR_df$Treatment <- metadata_coastR$Treatment

# --- 2. Variance explained ---
pca_var_perc_R <- round(summary(pca_result_coastR)$importance[2, 1:2] * 100, 1)

# --- 3. Centroids ---
centroids <- pca_coastR_df %>%
  group_by(Region) %>%
  summarize(PC1 = mean(PC1), PC2 = mean(PC2))

# --- 4. PCA plot with labels ---
pca_plot_R <- ggplot(pca_coastR_df, aes(x = PC1, y = PC2)) +
  geom_point(
    aes(color = Region, shape = Treatment),
    size = 3,
    position = position_jitter(width = 0.1, height = 0.1),
    alpha = 0.8
  ) +
  geom_text_repel(
    aes(label = Sample),
    size = 3,
    box.padding = 0.4,
    point.padding = 0.2,
    segment.size = 0.2,
    max.overlaps = Inf
  ) +
  geom_point(
    data = centroids,
    aes(x = PC1, y = PC2, color = Region),
    size = 5, shape = 4, stroke = 1.2, show.legend = FALSE
  ) +
  stat_ellipse(
    aes(group = Region, color = Region),
    level = 0.95,
    linetype = "dashed",
    linewidth = 0.6
  ) +
  scale_color_manual(
    values = c("coastR" = "grey20",
               "coastL" = "grey50",
               "inlandL" = "grey80")
  ) +
  scale_shape_manual(values = c("AW" = 16, "SW" = 17)) +
  labs(
    x = paste0("PC1 (", pca_var_perc_R[1], "%)"),
    y = paste0("PC2 (", pca_var_perc_R[2], "%)"),
    color = "Region",
    shape = "Treatment"
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    legend.background = element_blank(),
    legend.key = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    aspect.ratio = 1
  )

# Save
ggsave(
  "DE_Results_limma_voom_2025/real_annotation/PCA_plot_coastR_with_labels.pdf",
  plot = pca_plot_R,
  width = 8, height = 6,
  device = "pdf"
)
####PCA with no model for GBL assign in coastL##########
### 1. DATA LOADING AND PREPROCESSING ########################################

# Import metadata
metadata_coastL <- read.delim("PS_gill_RNA_metafile3.txt", sep = "\t") %>%
  filter(!Sample %in% c("CR4SW", "LCR1SW")) %>%
  mutate(
    Region = factor(Region, levels = c("coastR", "coastL", "inlandL")),
    Treatment = factor(Treatment, levels = c("AW", "SW")),
    Population = factor(Population)
  )

# Import and filter counts
TPM_quantmerge <- read.delim("DE_Results_limma_voom_2025/real_annotation/salmon_tx_TPM.txt", row.names = 1)
quantmerge <- read.delim("DE_Results_limma_voom_2025/real_annotation/salmon_tx_NumReads.txt", row.names = 1) %>%
  .[rowMeans(TPM_quantmerge) > 1.0, metadata_coastL$Sample] %>%
  as.matrix() %>% `mode<-`("integer") %>% replace(is.na(.), 0)

# Get common samples after filtering
common_samples <- intersect(colnames(quantmerge), metadata_coastL$Sample)
cat("Found", length(common_samples), "common samples\n")
# Subset and reorder both objects
quantmerge <- quantmerge[, common_samples, drop = FALSE]
metadata_coastL <- metadata_coastL %>% 
  filter(Sample %in% common_samples) %>%
  arrange(match(Sample, common_samples)) %>%
  column_to_rownames("Sample")

# Verify alignment
stopifnot(identical(colnames(quantmerge), rownames(metadata_coastL)))
### 4. Create DGEList ###

dge_coastL_pca <- DGEList(counts = quantmerge)
cat("DGEList samples:", ncol(dge_coastL_pca), "vs metadata:", nrow(metadata_coastL), "\n")
stopifnot(identical(colnames(dge_coastL_pca$counts), rownames(metadata_coastL)))

### 5. Normalize ###
dge_coastL <- calcNormFactors(dge_coastL_pca)

cat("=== BEFORE FILTERING ===\n")
cat("Total genes:", nrow(dge_coastL), "\n")
cat("Zero-count genes:", sum(rowSums(dge_coastL$counts) == 0), "\n")

### Expression filter ###
keep_coastL <- rowSums(dge_coastL$counts > 0) >= 0.10 * ncol(dge_coastL)
cat("Genes passing expression filter:", sum(keep_coastL), "\n")
dge_coastL_filtered <- dge_coastL[keep_coastL, , keep.lib.sizes = FALSE]

cat("=== AFTER FILTERING ===\n")
cat("Genes remaining:", nrow(dge_coastL_filtered), "\n")
cat("Percent retained:", round(nrow(dge_coastL_filtered)/nrow(dge_coastL)*100, 2), "%\n")
####For PCA analysis, no model included##############

vobj_coastL_pca <- voom(dge_coastL_filtered, design = NULL, plot = FALSE)
pca_res_CL <- prcomp(t(vobj_coastL_pca$E), center = TRUE, scale. = FALSE)
pc_scores_CL <- as.data.frame(pca_res_CL$x)
head(pc_scores_CL)
write.csv(pc_scores_CL, "DE_Results_limma_voom_2025/real_annotation/DE_coastR/pca_scores_CL.csv", row.names = TRUE)

# PCA on voom-normalized logCPM
# Calculate PCA variance percentages
pca_var_perc_CL <- round(summary(pca_res_CL)$importance[2, 1:2] * 100, 1)
pca_CL_df <- as.data.frame(pca_res_CL$x[, 1:2])
#Build PCA data frame with metadata
pca_CL_df$Sample    <- rownames(pca_CL_df)
pca_CL_df$Region <- metadata_coastL[rownames(pca_CL_df), "Region"]
pca_CL_df$Treatment <- metadata_coastL[rownames(pca_CL_df), "Treatment"]
pca_CL_df$Population <- metadata_coastL[rownames(pca_CL_df), "Population"]

#Variance explained
pca_var_perc_CL <- round(summary(pca_res_CL)$importance[2, 1:2] * 100, 1)
#Compute centroids (by Region)
centroids <- pca_CL_df %>%
  group_by(Region) %>%
  summarize(
    PC1 = mean(PC1),
    PC2 = mean(PC2),
    .groups = "drop"
  )
#PCA plot with labels, centroids, and 95% ellipses
pca_plot_CL <- ggplot(pca_CL_df, aes(x = PC1, y = PC2)) +
  
  # Points
  geom_point(
    aes(color = Region, shape = Treatment),
    size = 3,
    alpha = 0.85,
    position = position_jitter(width = 0.1, height = 0.1)
  ) +
  
  # Sample labels
  geom_text_repel(
    aes(label = Sample),
    size = 3,
    box.padding = 0.4,
    point.padding = 0.25,
    segment.size = 0.2,
    max.overlaps = Inf
  ) +
  
  # Region centroids
  geom_point(
    data = centroids,
    aes(x = PC1, y = PC2, color = Region),
    size = 5,
    shape = 4,
    stroke = 1.3,
    show.legend = FALSE
  ) +
  
  # 95% confidence ellipses
  stat_ellipse(
    aes(group = Region, color = Region),
    level = 0.95,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  
  # Labels
  labs(
    x = paste0("PC1 (", pca_var_perc_CL[1], "%)"),
    y = paste0("PC2 (", pca_var_perc_CL[2], "%)"),
    color = "Region",
    shape = "Treatment",
    title = "PCA of Gill RNA-seq"
  ) +
  
  # Styling
  scale_color_manual(
    values = c(
      "coastR"  = "grey20",
      "coastL"  = "grey50",
      "inlandL" = "grey80"
    )
  ) +
  scale_shape_manual(values = c("AW" = 16, "SW" = 17)) +
  
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    aspect.ratio = 1
  )
#Save plot
ggsave(
  "DE_Results_limma_voom_2025/real_annotation/PCA_coastL_with_centroids_ellipses_labels_no_model.pdf",
  plot = pca_plot_CL,
  width = 8,
  height = 6
)
######
######
###Warning message: In .standard_transform(ret) : No testable fixed effects were included in the model.
#Running topTable() will fail.
##This is because there is only population in coastL under SW, the Region and interaction effects cannot be statistically separated from the Population random effect###
###Need to block population, account for the population effect by assuming the correlation between populations and region are all the same
#####################################################
#####New analysis, real DE analysis #################
#####################################################
# Load required packages
library(limma)
library(edgeR)
library(variancePartition)
library(dplyr)
library(ggplot2)
library(gplots)
library(VennDiagram)
library(scales)
library(tibble)
setwd("/Users/Sherry/Desktop/RNA_bioinformatics")
### 1. DATA IMPORT AND PREPROCESSING ##########################################
### 1. Load metadata ###
metadata_coastR <- read.delim("PS_gill_RNA_metafile_new_region.txt", sep = "\t") %>%
  filter(!Sample %in% c("CR4SW", "LCR1SW")) %>%
  mutate(
    Region = factor(Region, levels = c("coastR", "coastL", "inlandL")),
    Treatment = factor(Treatment, levels = c("AW", "SW")),
    Population = factor(Population)
  )

### 2. Load counts ###
TPM_quantmerge <- read.delim("DE_Results_limma_voom_2025/real_annotation/salmon_tx_TPM.txt", row.names = 1)

quantmerge <- read.delim("DE_Results_limma_voom_2025/real_annotation/salmon_tx_NumReads.txt", row.names = 1) %>%
  .[rowMeans(TPM_quantmerge) > 1.0, metadata_coastR$Sample] %>%
  as.matrix() %>% 
  `mode<-`("integer") %>% 
  replace(is.na(.), 0)

### 3. Align samples ###
common_samples <- intersect(colnames(quantmerge), metadata_coastR$Sample)
cat("Found", length(common_samples), "common samples\n")

quantmerge <- quantmerge[, common_samples, drop = FALSE]

metadata_coastR <- metadata_coastR %>%
  filter(Sample %in% common_samples) %>%
  column_to_rownames("Sample")

stopifnot(identical(colnames(quantmerge), rownames(metadata_coastR)))

### 4. Create DGEList ###
dge_coastR <- DGEList(counts = quantmerge)
cat("DGEList samples:", ncol(dge_coastR), "vs metadata:", nrow(metadata_coastR), "\n")

stopifnot(identical(colnames(dge_coastR$counts), rownames(metadata_coastR)))

### 5. Normalize ###
dge_coastR <- calcNormFactors(dge_coastR)

cat("=== BEFORE FILTERING ===\n")
cat("Total genes:", nrow(dge_coastR), "\n")
cat("Zero-count genes:", sum(rowSums(dge_coastR$counts) == 0), "\n")

### Expression filter ###
keep_coastR <- rowSums(dge_coastR$counts > 0) >= 0.10 * ncol(dge_coastR)
cat("Genes passing expression filter:", sum(keep_coastR), "\n")

dge_coastR_filtered <- dge_coastR[keep_coastR, , keep.lib.sizes = FALSE]

cat("=== AFTER FILTERING ===\n")
cat("Genes remaining:", nrow(dge_coastR_filtered), "\n")
cat("Percent retained:", round(nrow(dge_coastR_filtered)/nrow(dge_coastR)*100, 2), "%\n")
####For PCA analysis, no model included##############

vobj_coastR_pca <- voom(dge_coastR_filtered, design = NULL, plot = FALSE)
pca_res_pca <- prcomp(t(vobj_coastR_pca$E), center = TRUE, scale. = FALSE)
pc_scores <- as.data.frame(pca_res_pca$x)
head(pc_scores)
write.csv(pc_scores, "DE_Results_limma_voom_2025/real_annotation/DE_coastR/pca_scores.csv", row.names = TRUE)
### 6. DESIGN AND DE ANALYSIS WITH POPULATION BLOCKING #########################
# 6.1 Design matrix for fixed effects: Region, Treatment, Region:Treatment
design_coastR <- model.matrix(~ Region * Treatment, data = metadata_coastR)
rownames(design_coastR) <- rownames(metadata_coastR)

# 6.2 voom transformation (logCPM, variance modeling)
vobj_coastR <- voom(dge_coastR_filtered, design = design_coastR, plot = TRUE)
expression_real_anno <- vobj_coastR$E
expr_out <- cbind(Gene = rownames(expression_real_anno),
                  expression_real_anno)
write.table(expr_out, "DE_Results_limma_voom_2025/real_annotation/DE_coastR/norm_couts_real_annot.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
# 6.3 Estimate intra-population correlation using duplicateCorrelation
# This accounts for non-independence of samples within populations
corfit <- duplicateCorrelation(vobj_coastR, design_coastR, block = metadata_coastR$Population)
cat("Consensus correlation among populations:", corfit$consensus, "\n")

# 6.4 Fit linear model with population as blocking factor
fit_coastR <- lmFit(vobj_coastR, design_coastR, block = metadata_coastR$Population,
                    correlation = corfit$consensus)

cat("Available coefficients:\n")
print(colnames(coef(fit_coastR)))

### 6.5 Empirical Bayes moderation
fit_coastR <- eBayes(fit_coastR)

# 6.6 Define coefficients for joint F-tests
region_coefs <- c("RegioncoastL", "RegioninlandL")  # Region main effects
interaction_coefs <- c("RegioncoastL:TreatmentSW", "RegioninlandL:TreatmentSW")  # Interaction effects

# Find indices of coefficients
coef_idx <- match(region_coefs, colnames(coef(fit_coastR)))

# Joint F-test
region_F <- topTable(fit_coastR, coef = coef_idx, number = Inf, sort.by = "F", adjust.method = "BH")

# Similarly for interaction
interaction_idx <- match(interaction_coefs, colnames(coef(fit_coastR)))
interaction_F <- topTable(fit_coastR, coef = interaction_idx, number = Inf, sort.by = "F", adjust.method = "BH")

# 6.8 t-test for salinity/Treatment effect
salinity_t <- topTable(fit_coastR, coef = "TreatmentSW", number = Inf, adjust.method = "BH")

# Optional: save results
write.csv(region_F, "DE_Results_limma_voom_2025/real_annotation/DE_coastR/Region_Ftest.csv", row.names = TRUE)
write.csv(interaction_F, "DE_Results_limma_voom_2025/real_annotation/DE_coastR/Region_Treatment_interaction_Ftest.csv", row.names = TRUE)
write.csv(salinity_t, "DE_Results_limma_voom_2025/real_annotation/DE_coastR/Salinity_Ttreatment_Ttest.csv", row.names = TRUE)

# Coefficients for Region main effect

region_coefs <- c("RegioncoastL", "RegioninlandL")

# Assume you already have:
# salinity_t <- topTable(fit_coastR, coef = "TreatmentSW", number = Inf, adjust.method = "BH")

# 1. Filter for significant DE genes (FDR < 0.05)
sig_salinity <- salinity_t[salinity_t$adj.P.Val < 0.05, ]

# 2. Optional: add gene names if not already rownames

sig_salinity$Gene <- rownames(sig_salinity)

# 3. Save to CSV
write.csv(sig_salinity, "DE_Results_limma_voom_2025/real_annotation/DE_coastR/Salinity_SWvsAW_sig_DE_genes.csv", row.names = FALSE)

# region_F <- topTable(fit_coastR, coef = region_coefs, number = Inf, test = "F", adjust.method = "BH")

# 1. Filter by FDR < 0.05
sig_region <- region_F[region_F$adj.P.Val < 0.05, ]

# 2. Add gene names
sig_region$Gene <- rownames(sig_region)

# 3. Save
write.csv(sig_region, "DE_Results_limma_voom_2025/real_annotation/DE_coastR/Region_Ftest_sig_genes.csv", row.names = FALSE)

sig_interaction <- interaction_F[interaction_F$adj.P.Val < 0.05, ]
sig_interaction$Gene <- rownames(sig_interaction)
write.csv(sig_interaction, "DE_Results_limma_voom_2025/real_annotation/DE_coastR/Region_Treatment_interaction_sig_genes.csv", row.names = FALSE)

### regional F test within each salinity treatment####################################
# Recreate design with valid names
design_coastR <- model.matrix(~ Region * Treatment, data = metadata_coastR)
colnames(design_coastR) <- make.names(colnames(design_coastR))

# Re-run voom (no need to plot again)
vobj_coastR <- voom(dge_coastR_filtered, design_coastR, plot = FALSE)

# Re-estimate correlation
corfit <- duplicateCorrelation(
  vobj_coastR,
  design_coastR,
  block = metadata_coastR$Population
)

# Refit model
fit_coastR <- lmFit(
  vobj_coastR,
  design_coastR,
  block = metadata_coastR$Population,
  correlation = corfit$consensus
)

fit_coastR <- eBayes(fit_coastR)
colnames(design_coastR)

##within AW####
region_AW_idx <- match(
  c("RegioncoastL", "RegioninlandL"),
  colnames(coef(fit_coastR))
)

region_AW_F <- topTable(
  fit_coastR,
  coef = region_AW_idx,
  number = Inf,
  sort.by = "F",
  adjust.method = "BH"
)

write.csv(region_AW_F, "DE_Results_limma_voom_2025/real_annotation/DE_coastR/region_AW_F.csv", row.names = TRUE)

sig_region_AW_F <- region_AW_F[region_AW_F$adj.P.Val < 0.05, ]
sig_region_AW_F$Gene <- rownames(sig_region_AW_F)
write.csv(sig_region_AW_F, "DE_Results_limma_voom_2025/real_annotation/DE_coastR/sig_region_AW_F.csv", row.names = FALSE)
###within SW####
contrast_SW_region <- makeContrasts(
  coastL_vs_coastR_SW  = RegioncoastL + RegioncoastL.TreatmentSW,
  inlandL_vs_coastR_SW = RegioninlandL + RegioninlandL.TreatmentSW,
  levels = colnames(coef(fit_coastR))
)

fit_SW <- contrasts.fit(fit_coastR, contrast_SW_region)
fit_SW <- eBayes(fit_SW)

region_SW_F <- topTable(
  fit_SW,
  coef = 1:2,        # both contrasts jointly
  number = Inf,
  sort.by = "F",
  adjust.method = "BH"
)

write.csv(region_SW_F, "DE_Results_limma_voom_2025/real_annotation/DE_coastR/region_SW_F.csv", row.names = TRUE)

sig_region_SW_F <- region_SW_F[region_SW_F$adj.P.Val < 0.05, ]
sig_region_SW_F$Gene <- rownames(sig_region_SW_F)
write.csv(sig_region_SW_F, "DE_Results_limma_voom_2025/real_annotation/DE_coastR/sig_region_SW_F.csv", row.names = FALSE)
### 7. CONTRASTS: Salinity effect within each Region #########################
library(limma)

# 7.1 Check coefficient names
colnames(design_coastR)
# 7.2 Define contrasts: SW vs AW within each Region
contrast_matrix <- makeContrasts(
  SW_vs_AW_coastR = TreatmentSW,
  SW_vs_AW_coastL = TreatmentSW + `RegioncoastL.TreatmentSW`,
  SW_vs_AW_inlandL = TreatmentSW + `RegioninlandL.TreatmentSW`,
  levels = design_coastR
)

# 7.3 Fit contrasts
fit2 <- contrasts.fit(fit_coastR, contrast_matrix)
fit2 <- eBayes(fit2)

# 7.4 Extract DE genes for each Region
sw_aw_coastR <- topTable(fit2, coef = "SW_vs_AW_coastR", number = Inf, adjust.method = "BH")
sw_aw_coastL <- topTable(fit2, coef = "SW_vs_AW_coastL", number = Inf, adjust.method = "BH")
sw_aw_inlandL <- topTable(fit2, coef = "SW_vs_AW_inlandL", number = Inf, adjust.method = "BH")

# 7.5 Optional: save results
write.csv(sw_aw_coastR, "DE_Results_limma_voom_2025/real_annotation/DE_coastR/Salinity_SWvsAW_coastR.csv", row.names = TRUE)
write.csv(sw_aw_coastL, "DE_Results_limma_voom_2025/real_annotation/DE_coastR/Salinity_SWvsAW_coastL.csv", row.names = TRUE)
write.csv(sw_aw_inlandL, "DE_Results_limma_voom_2025/real_annotation/DE_coastR/Salinity_SWvsAW_inlandL.csv", row.names = TRUE)

# SW vs AW within coastR
sig_sw_coastR <- sw_aw_coastR[sw_aw_coastR$adj.P.Val < 0.05, ]
sig_sw_coastR$Gene <- rownames(sig_sw_coastR)
write.csv(sig_sw_coastR, "DE_Results_limma_voom_2025/real_annotation/DE_coastR/Salinity_SWvsAW_coastR_sig_DE_genes.csv", row.names = FALSE)

# SW vs AW within coastL
sig_sw_coastL <- sw_aw_coastL[sw_aw_coastL$adj.P.Val < 0.05, ]
sig_sw_coastL$Gene <- rownames(sig_sw_coastL)
write.csv(sig_sw_coastL, "DE_Results_limma_voom_2025/real_annotation/DE_coastR/Salinity_SWvsAW_coastL_sig_DE_genes.csv", row.names = FALSE)

# SW vs AW within inlandL
sig_sw_inlandL <- sw_aw_inlandL[sw_aw_inlandL$adj.P.Val < 0.05, ]
sig_sw_inlandL$Gene <- rownames(sig_sw_inlandL)
write.csv(sig_sw_inlandL, "DE_Results_limma_voom_2025/real_annotation/DE_coastR/Salinity_SWvsAW_inlandL_sig_DE_genes.csv", row.names = FALSE)

########Venn diagram####################
sig_treatment <- read.csv("DE_Results_limma_voom_2025/real_annotation/DE_coastR/Salinity_SWvsAW_sig_DE_genes.csv")$Gene
sig_habitat <- read.csv("DE_Results_limma_voom_2025/real_annotation/DE_coastR/Region_Ftest_sig_genes.csv")$Gene
sig_interaction <- read.csv("DE_Results_limma_voom_2025/real_annotation/DE_coastR/Region_Treatment_interaction_sig_genes.csv")$Gene
library(eulerr)
library(polylabelr)
venn_data <- list(
  Treatment = sig_treatment,
  Habitat = sig_habitat,
  Interaction = sig_interaction
)

fit <- euler(venn_data)  # Computes proportional areas
##show counts instead of percentage###########
plot(fit,
     quantities = list(type = "counts", fontsize = 9, fontface = "bold"),
     fills = c("#1F77B4", "#FF7F0E", "#2CA02C"),  # Matplotlib default colors
     alpha = 0.6,
     labels = list(fontsize = 11, fontfamily = "sans"))
##If you get warnings about "missing intersections", force exact calculations:
fit <- euler(venn_data, input = "disjoint", shape = "ellipse")
### Save Publication-Ready Figure

pdf("DE_Results_limma_voom_2025/real_annotation/DE_coastR/Venn_diagram_coastR.pdf", width = 6, height = 6)
plot(fit,
     quantities = list(type = "counts", fontsize = 9),
     fills = c("blue", "red", "green"),
     labels = list(fontsize = 11))
dev.off()

#######using this version model, with populations counted as random effect##############################
### 9. Save PCA input ###
# PCA on voom-normalized logCPM
pca_res <- prcomp(t(vobj_coastR$E), center = TRUE, scale. = FALSE)
pca_data_coastR <- vobj_coastR$E
write.csv(pca_data_coastR, 
         "DE_Results_limma_voom_2025/real_annotation/DE_coastR/pca_data_coastR_new_modle.csv")
### 10. Background genes ###
background_genes_coastR <- rownames(dge_coastR_filtered)
write.table(background_genes_coastR,
            "DE_Results_limma_voom_2025/real_annotation/DE_coastR/background_gene_set_coastR_new_model.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
### 11. PCA ###
#library(ggplot2)
#library(stringr)
#library(dplyr)
#library(ggrepel)
vobj_coastR_pca <- voom(dge_coastR_filtered, design = NULL, plot = FALSE)
pca_res_pca <- prcomp(t(vobj_coastR_pca$E), center = TRUE, scale. = FALSE)
pc_scores <- as.data.frame(pca_res_pca$x)
head(pc_scores)
# PCA on voom-normalized logCPM
#pca_res <- prcomp(t(vobj_coastR$E), center = TRUE, scale. = FALSE)
# Calculate PCA variance percentages
pca_var_perc_pca <- round(summary(pca_res_pca)$importance[2, 1:2] * 100, 1)
pca_pca_df <- as.data.frame(pca_res_pca$x[, 1:2])
#Build PCA data frame with metadata
pca_pca_df$Sample    <- rownames(pca_pca_df)
pca_pca_df$Region <- metadata_coastR[rownames(pca_pca_df), "Region"]
pca_pca_df$Treatment <- metadata_coastR[rownames(pca_pca_df), "Treatment"]
pca_pca_df$Population <- metadata_coastR[rownames(pca_pca_df), "Population"]

#Variance explained
pca_var_perc_pca <- round(summary(pca_res_pca)$importance[2, 1:2] * 100, 1)
#Compute centroids (by Region)
centroids <- pca_pca_df %>%
  group_by(Region) %>%
  summarize(
    PC1 = mean(PC1),
    PC2 = mean(PC2),
    .groups = "drop"
  )
#PCA plot with labels, centroids, and 95% ellipses
pca_plot_pca <- ggplot(pca_pca_df, aes(x = PC1, y = PC2)) +
  
  # Points
  geom_point(
    aes(color = Region, shape = Treatment),
    size = 3,
    alpha = 0.85,
    position = position_jitter(width = 0.1, height = 0.1)
  ) +
  
  # Sample labels
  geom_text_repel(
    aes(label = Sample),
    size = 3,
    box.padding = 0.4,
    point.padding = 0.25,
    segment.size = 0.2,
    max.overlaps = Inf
  ) +
  
  # Region centroids
  geom_point(
    data = centroids,
    aes(x = PC1, y = PC2, color = Region),
    size = 5,
    shape = 4,
    stroke = 1.3,
    show.legend = FALSE
  ) +
  
  # 95% confidence ellipses
  stat_ellipse(
    aes(group = Region, color = Region),
    level = 0.95,
    linetype = "dashed",
    linewidth = 0.7
  ) +
  
  # Labels
  labs(
    x = paste0("PC1 (", pca_var_perc_pca[1], "%)"),
    y = paste0("PC2 (", pca_var_perc_pca[2], "%)"),
    color = "Region",
    shape = "Treatment",
    title = "PCA of Gill RNA-seq"
  ) +
  
  # Styling
  scale_color_manual(
    values = c(
      "coastR"  = "grey20",
      "coastL"  = "grey50",
      "inlandL" = "grey80"
    )
  ) +
  scale_shape_manual(values = c("AW" = 16, "SW" = 17)) +
  
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    panel.grid.major = element_line(color = "grey90", linewidth = 0.2),
    aspect.ratio = 1
  )
#Save plot
ggsave(
  "DE_Results_limma_voom_2025/real_annotation/DE_coastR/PCA_coastR_with_centroids_ellipses_labels_no_model.pdf",
  plot = pca_plot_pca,
  width = 8,
  height = 6
)

table(pca_df$Region)
###############################################################################################################

### 8. HEATMAPS FOR INTERACTION GENES ########################################

create_interaction_heatmaps <- function(vobj, metadata_df, DE_results, assignment_name) {
  cat("\n=== CREATING INTERACTION HEATMAPS -", assignment_name, "===\n")
  
  # Get genes with significant interaction effects
  interaction_genes <- unique(c(
    rownames(DE_results$Interaction_coastR)[DE_results$Interaction_coastR$adj.P.Val < 0.05],
    rownames(DE_results$Interaction_coastL)[DE_results$Interaction_coastL$adj.P.Val < 0.05],
    rownames(DE_results$Interaction_inlandL)[DE_results$Interaction_inlandL$adj.P.Val < 0.05]
  ))
  
  if(length(interaction_genes) == 0) {
    cat("No significant interaction genes found for", assignment_name, "\n")
    return(NULL)
  }
  
  # Limit to top 50 genes for clarity
  top_interaction_genes <- head(interaction_genes, 50)
  expr_data <- vobj$E[top_interaction_genes, ]
  
  # Create annotation
  annotation_col <- data.frame(
    Region = metadata_df$Region,
    Treatment = metadata_df$Treatment,
    Population = metadata_df$Population
  )
  rownames(annotation_col) <- rownames(metadata_df)
  
  # Colors
  region_colors <- c(coastR = "#E41A1C", coastL = "#377EB8", inlandL = "#4DAF4A")
  treatment_colors <- c(AW = "#FBB4AE", SW = "#B3CDE3")
  
  annotation_colors <- list(
    Region = region_colors,
    Treatment = treatment_colors
  )
  
  # Create heatmap
  pheatmap(
    expr_data,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    scale = "row",
    show_rownames = FALSE,
    show_colnames = TRUE,
    cluster_cols = TRUE,
    main = paste("Region × Treatment Interaction Genes - GBL in", assignment_name),
    color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(100)
  )
  
  cat("Created interaction heatmap with", length(top_interaction_genes), "genes for", assignment_name, "\n")
  return(top_interaction_genes)
}

# Create interaction heatmaps for both assignments
interaction_genes_coastL <- create_interaction_heatmaps(vobj_coastL, metadata_coastL, DE_coastL, "coastL")
interaction_genes_coastR <- create_interaction_heatmaps(vobj_coastR, metadata_coastR, DE_coastR, "coastR")

### 9. SAVE COMPREHENSIVE RESULTS ###########################################

cat("\n=== SAVING COMPREHENSIVE RESULTS ===\n")

saveRDS(list(
  DE_coastL = DE_coastL,
  DE_coastR = DE_coastR,
  metadata_coastL = metadata_coastL,
  metadata_coastR = metadata_coastR,
  vobj_coastL = vobj_coastL,
  vobj_coastR = vobj_coastR,
  interaction_genes_coastL = interaction_genes_coastL,
  interaction_genes_coastR = interaction_genes_coastR
), file = "GBL_assignment_DE_with_interactions.rds")

# Save detailed interaction results
write.csv(DE_coastL$Interaction_coastR_vs_coastL, "Interaction_coastR_vs_coastL_GBL_in_coastL.csv")
write.csv(DE_coastR$Interaction_coastR_vs_coastL, "Interaction_coastR_vs_coastL_GBL_in_coastR.csv")

cat("DE analysis with interaction testing completed!\n")
cat("Key files saved:\n")
cat("- GBL_assignment_DE_with_interactions.rds (complete results)\n")
cat("- Interaction_coastR_vs_coastL_GBL_in_coastL.csv\n")
cat("- Interaction_coastR_vs_coastL_GBL_in_coastR.csv\n")
cat("\nCheck the interaction term results to see which assignment gives more biologically meaningful interaction patterns!\n")

#####try heatmeaps for regional effect with the new GBL placement#####

##1. Load and prepare data
library(tidyverse)
library(cowplot)
library(viridis)
# Read significant habitat genes
# ---------------------------
# Load significant genes
# ---------------------------
sig_genes <- read_csv(
  "DE_Results_limma_voom_2025/new_arrange/again/sig_habitat_joint_Ftest.csv",
  show_col_types = FALSE
)

# ---------------------------
# Load normalized expression (voom logCPM)
# ---------------------------
logCPM_matrix <- vobj_coastR$E

# Verify matching genes
genes_use <- intersect(sig_genes$Gene, rownames(logCPM_matrix))
cat(length(genes_use), "out of", nrow(sig_genes), "genes matched\n")

if (length(genes_use) == 0) {
  stop("No significant genes found in expression matrix")
}

# ---------------------------
# Region colors
# ---------------------------
region_colors <- c(
  coastR  = "#1f77b4",
  coastL  = "#aec7e8",
  inlandL = "#ff7f0e"
)

# ---------------------------
# Prepare heatmap dataframe
# ---------------------------
metadata_coastR_df <- metadata_coastR %>%
  tibble::rownames_to_column("Sample")

heatmap_data <- logCPM_matrix[genes_use, ] %>%
  t() %>%
  scale() %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  pivot_longer(
    -Gene,
    names_to = "Sample",
    values_to = "Zscore"
  ) %>%
  left_join(metadata_coastR_df, by = "Sample") %>%
  mutate(
    Region = factor(Region, levels = c("coastR", "coastL", "inlandL")),
    Sample = factor(Sample, levels = unique(Sample))
  )

# ---------------------------
# Annotation bar (Region)
# ---------------------------
region_anno <- heatmap_data %>%
  distinct(Sample, Region) %>%
  ggplot(aes(x = Sample, y = 1, fill = Region)) +
  geom_tile() +
  scale_fill_manual(values = region_colors) +
  theme_void() +
  theme(
    legend.position = "top",
    legend.title = element_blank()
  )

# ---------------------------
# Gene ordering (mean Z-score)
# ---------------------------
gene_order <- heatmap_data %>%
  group_by(Gene) %>%
  summarize(mean_z = mean(Zscore), .groups = "drop") %>%
  arrange(mean_z) %>%
  pull(Gene)

# ---------------------------
# Main heatmap
# ---------------------------
main_heatmap <- ggplot(
  heatmap_data,
  aes(
    x = Sample,
    y = factor(Gene, levels = gene_order),
    fill = Zscore
  )
) +
  geom_tile() +
  scale_fill_viridis(
    option = "magma",
    limits = c(-3, 3),
    oob = scales::squish,
    name = "Z-score"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(x = NULL, y = NULL)

# ---------------------------
# Combine annotation + heatmap
# ---------------------------
final_plot <- plot_grid(
  region_anno,
  main_heatmap,
  ncol = 1,
  align = "v",
  rel_heights = c(0.08, 1)
)

# ---------------------------
# Save
# ---------------------------
ggsave(
  "DE_Results_limma_voom_2025/new_arrange/habitat_heatmap_publication_ready.pdf",
  final_plot,
  width = 16,
  height = 12,
  device = "pdf"
)
### a 2nd version####
library(pheatmap)
library(RColorBrewer)
library(dplyr)

# ---------------------------
# 1. Load significant genes
# ---------------------------
sig_genes <- readr::read_csv(
  "DE_Results_limma_voom_2025/new_annotation/DE_coastR/Region_Ftest_sig_genes.csv",
  show_col_types = FALSE
)

# ---------------------------
# 2. Load voom logCPM
# ---------------------------
logCPM <- vobj_coastR$E

genes_use <- intersect(sig_genes$Gene, rownames(logCPM))
cat(length(genes_use), "genes matched\n")
stopifnot(length(genes_use) > 0)

# ---------------------------
# 3. Subset & row Z-score
# ---------------------------
mat <- logCPM[genes_use, ]
mat_z <- t(scale(t(mat)))   # NO clipping, NO limits

# ---------------------------
# 4. Column annotation (Region only)
# ---------------------------
annotation_col$Region <- as.character(annotation_col$Region)

annotation_col <- metadata_coastR %>%
  dplyr::select(Region, Treatment) %>%
  as.data.frame()


annotation_col <- annotation_col[colnames(mat_z), , drop = FALSE]

# ---------------------------
# 5. Annotation colors (blue / green / red)
# ---------------------------
annotation_colors <- list(
  Region = c(
    coastR  = "#1f78b4",  # blue
    coastL  = "#33a02c",  # green
    inlandL = "#e31a1c"   # red
  )
)
# Clip Z-scores for visualization
mat_z[mat_z >  3] <-  3
mat_z[mat_z < -3] <- -3
# ---------------------------
# 6. Heatmap color palette
#    lighter orange → white → blue
# ---------------------------
heatmap_colors <- colorRampPalette(
  c("orange", "white", "blue")
)(101)

# ---------------------------
# 7. Draw heatmap
# ---------------------------
pheatmap(
  mat_z,
  color              = heatmap_colors,
  annotation_col     = annotation_col,
  annotation_colors  = annotation_colors,
  cluster_rows       = TRUE,
  cluster_cols       = FALSE,
  show_rownames      = FALSE,
  show_colnames      = TRUE,
  border_color       = NA,
  fontsize_col       = 9,
  main               = "Omnibus Regional Differential Expression",
  filename           = "DE_Results_limma_voom_2025/new_annotation/habitat_heatmap_Dec.pdf",
  width              = 16,
  height             = 12
)

###build a heatmap for regional  effect stratified salinity, need to change################################

library(pheatmap)
library(RColorBrewer)
library(dplyr)
##Load file
AW.gene <- read.csv("DE_Results_limma_voom_2025/real_annotation/DE_coastR/sig_region_AW_F.csv") 
# Verify the structure
head(AW.gene)
SW.gene <- read.csv("DE_Results_limma_voom_2025/real_annotation/DE_coastR/sig_region_SW_F.csv")
head(SW.gene)
# Get all unique genes from both files
all_genes <- unique(c(AW.gene$Gene, SW.gene$Gene))
logCPM <- vobj_coastR$E

genes_use <- intersect(all_genes, rownames(logCPM))
cat(length(genes_use), "genes matched\n")
stopifnot(length(genes_use) > 0)

# ---------------------------
# 3. Subset & row Z-score
# ---------------------------
mat <- logCPM[genes_use, ]
mat_z <- t(scale(t(mat)))   # NO clipping, NO limits

# ---------------------------
# 4. Column annotation (Region only)
# ---------------------------
annotation_col <- metadata_coastR %>%
  dplyr::select(Region, Treatment) %>%
  as.data.frame()

# Ensure column order matches mat_z
annotation_col <- annotation_col[colnames(mat_z), , drop = FALSE]

# ---------------------------
# 5. Annotation colors (blue / green / red)
# ---------------------------
annotation_col$Region <- factor(annotation_col$Region,
                                levels = c("Coastal River", "Coastal Lake", "Interior Lake"),
                                labels = c("coastR", "coastL", "inlandL"))
annotation_colors <- list(
  Region = c(
    coastR = "#1f78b4",  # blue
    coastL  = "#33a02c",  # green
    inlandL = "#e31a1c"   # red
  ),
  Treatment = c(
    AW = "#f7f7f7",
    SW = "#cccccc"
  )
)
# Clip Z-scores for visualization
mat_z[mat_z >  3] <-  3
mat_z[mat_z < -3] <- -3
# ---------------------------
# 6. Heatmap color palette
#    lighter orange → white → blue
# ---------------------------
heatmap_colors <- colorRampPalette(
  c("orange", "white", "blue")
)(101)

# ---------------------------
# 7. Draw heatmap
# ---------------------------
pheatmap(
  mat_z,
  color             = heatmap_colors,
  annotation_col    = annotation_col,
  annotation_colors = annotation_colors,
  cluster_rows      = TRUE,
  cluster_cols      = FALSE,
  show_rownames     = FALSE,
  show_colnames     = TRUE,
  border_color      = NA,
  fontsize_col      = 9,
  main              = "Omnibus Regional Differential Expression",
  filename          = "DE_Results_limma_voom_2025/real_annotation/DE_coastR/habitat_heatmap_AW_SW.pdf",
  width             = 16,
  height            = 12
)

####GO enrichment analysis#####################################
##Get annotation, GO id for the background genes#######
background <- read.table("DE_Results_limma_voom_2025/real_annotation/background_gene_set_coastR.txt", header = FALSE, stringsAsFactors = FALSE)
colnames(background) <- "transcript_id"
gene2go <- read.table("DE_Results_limma_voom_2025/real_annotation/whole_annotation_transcript_to_go.tsv", sep = "\t", header = TRUE, stringsAsFactor = FALSE)
head(gene2go)
gene2go_bg <- gene2go[gene2go$transcript_id %in% background$transcript_id, ]
head(gene2go_bg)
write.table(gene2go_bg, "DE_Results_limma_voom_2025/real_annotation/background_gene_to_GOid.tsv", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
gene_len <- read.table("DE_Results_limma_voom_2025/real_annotation/gene_effective_lengths.tsv", header = TRUE, stringsAsFactor = FALSE)
background_gene_length <- gene_len[gene_len$Transcript %in% background$transcript_id, ]
head(background_gene_length)
nrow(background)
nrow(background_gene_length)
write.table(background_gene_length, "DE_Results_limma_voom_2025/real_annotation/background_gene_length.tsv", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
length(unique(gene2go_bg$transcript_id))
##To get the same genes with GOid and length##
gene2go_bg <- read.table("DE_Results_limma_voom_2025/real_annotation/background_gene_to_GOid.tsv", 
                         header = FALSE, stringsAsFactors = FALSE)
genes_with_go <- unique(gene2go_bg$V1)
background_gene_length <- read.table("DE_Results_limma_voom_2025/real_annotation/background_gene_length.tsv",
                                     header = FALSE, stringsAsFactors = FALSE)
genes_with_len <- background_gene_length$V1

valid_genes <- intersect(genes_with_go, genes_with_len)
length(valid_genes)
write.table(
  valid_genes, 
  file = "DE_Results_limma_voom_2025/real_annotation/valid_genes.txt", 
  quote = FALSE,      # no quotes around IDs
  row.names = FALSE,  # no row numbers
  col.names = FALSE   # no header
)
# Keep only genes with GO + length
gene2go_bg_filtered <- gene2go_bg[gene2go_bg$V1 %in% valid_genes, ]
length(unique(gene2go_bg_filtered$V1))
write.table(gene2go_bg_filtered, "DE_Results_limma_voom_2025/real_annotation/gene2go_bg_filtered.tsv", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
background_gene_length_filter <- background_gene_length[background_gene_length$V1 %in% valid_genes, ]
nrow(background_gene_length_filter)
write.table(background_gene_length_filter, "DE_Results_limma_voom_2025/real_annotation/background_gene_length_filtered.tsv", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
######GO enrichment for reginal effect#########
library(dplyr)
library(tidyr)
library(goseq) #AS the goseq can't be used here, the following work was done on HPC R4 conda environment

# 1. Load background genes
background_genes <- read.delim("DE_Results_limma_voom_2025/real_annotation/background_gene_set_coastR.txt", header = FALSE)
colnames(background_genes) <- "Transcript"
head(background_genes$Transcript)
str(background_genes)
# 2. Load gene lengths (ensure it has 'Transcript' and 'Length' columns)
#background_gene_length_filtered <- read.delim("DE_Results_limma_voom_2025/real_annotation/background_gene_length_filtered.tsv")
#rownames(background_gene_length_filtered) <- background_gene_length_filtered$Transcript  # Set row names
background_gene_length_filtered <- read.delim(
  "DE_Results_limma_voom_2025/real_annotation/background_gene_length_filtered.tsv",
  header = FALSE,
  stringsAsFactors = FALSE,
  col.names = c("Transcript", "Length")
)
head(background_gene_length_filtered)
#background_gene_length_filtered <- background_gene_length_filtered[, c("Transcript", "Length")]
#background_gene_length_filtered <- background_gene_length_filtered[, !duplicated(colnames(background_gene_length_filtered))]
#background_gene_length_filtered$Length <- as.numeric(background_gene_length_filtered$Length)
rownames(background_gene_length_filtered) <- background_gene_length_filtered$Transcript
background_gene_length_filtered$Transcript <- NULL
# 3. Load GO annotations (if needed later)
transcript_to_go <- read.delim("DE_Results_limma_voom_2025/real_annotation/gene2go_bg_filtered.tsv")
colnames(transcript_to_go) <- c("Transcript", "GO_ID", "source")
# 4. Load significant region genes (correct sep if CSV)
sig_region_genes_real <- read.delim(
  "DE_Results_limma_voom_2025/real_annotation/DE_coastR/Region_Ftest_sig_genes.csv",
  sep = ",",
  header = TRUE,
  stringsAsFactors = FALSE
)
# Clean whitespace from the gene column
sig_region_genes_real$Gene <- trimws(sig_region_genes_real$Gene)

head(sig_region_genes_real)

# 5. Create binary vector (1=DE, 0=not DE)
region_vector_real <- as.integer(background_genes$Transcript %in% sig_region_genes_real$Gene)
names(region_vector_real) <- background_genes$Transcript

# 7. Match gene lengths to habitat_vector_re
gene_lengths_matched <- background_gene_length_filtered[names(region_vector_real), "Length"]
str(gene_lengths_matched)
class(gene_lengths_matched)
# 8. Verify lengths match
stopifnot(length(region_vector_real) == length(gene_lengths_matched))
stopifnot(all(names(region_vector_real) == names(gene_lengths_matched)))

# 9. Run nullp()
region_real_pwf <- nullp(region_vector_real, bias.data = gene_lengths_matched)
region_go_real <- goseq(
  region_real_pwf,
  gene2cat = transcript_to_go,
  method = "Wallenius",
  use_genes_without_cat = FALSE  # More stringent
)
write.table(region_go_real, 
            "DE_Results_limma_voom_2025/real_annotation/DE_coastR/region_go_real.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

###########################################################################################3
############KEGG pathway enrichment analysis################################################
library(KEGGREST)
library(dplyr)
library(tidyr)
library(pbapply)
uni_kegg <-read.csv("DE_Results_limma_voom_2025/real_annotation/DE_coastR/uni_kegg.csv")
head(uni_kegg)
uni_kegg_clean <- uni_kegg |>
  filter(KEGG != "")
nrow(uni_kegg)
nrow(uni_kegg_clean)
head(uni_kegg_clean)
write.table(uni_kegg_clean, "DE_Results_limma_voom_2025/real_annotation/DE_coastR/uni_kegg_clean.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
combine_anno <- read.delim("DE_Results_limma_voom_2025/real_annotation/combined_annotations_outer_clean.tsv")
head(combine_anno)
colnames(uni_kegg_clean)
uni_kegg_clean <- uni_kegg_clean %>%
  rename(uniprot_accession = Uniprot_accession)
background_gene_uniprot_kegg <- left_join(combine_anno, uni_kegg_clean, by="uniprot_accession")
head(background_gene_uniprot_kegg)
write.table(background_gene_uniprot_kegg, "DE_Results_limma_voom_2025/real_annotation/DE_coastR/background_gene_uniprot_kegg.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
background_gene_uniprot_kegg_clean <- background_gene_uniprot_kegg |>
  filter(KEGG !="")
nrow(background_gene_uniprot_kegg)
nrow(background_gene_uniprot_kegg_clean)
background_gene_uniprot_kegg_clean2 <- background_gene_uniprot_kegg_clean %>%
  select(protein_id, KEGG)
head(background_gene_uniprot_kegg_clean2)
write.table(background_gene_uniprot_kegg_clean2, "DE_Results_limma_voom_2025/real_annotation/DE_coastR/background_gene_uniprot_kegg_clean.tsv", sep = "\t",row.names = FALSE, quote = FALSE )
# Split multiple KEGG IDs into separate rows
proteinid_kegg_long <- background_gene_uniprot_kegg_clean2 %>%
  separate_rows(KEGG, sep = ";") %>%   # split on ";"
  filter(KEGG != "")                   # remove any leftover empty strings
head(proteinid_kegg_long)
write.table(proteinid_kegg_long, "DE_Results_limma_voom_2025/real_annotation/DE_coastR/proteinid_kegg_long.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
KEGG_ids <- unique(proteinid_kegg_long$KEGG)
head(KEGG_ids)
# Function to query KO number for a given KEGG gene ID
getKO <- function(gene_id) {
  url <- paste0("http://rest.kegg.jp/get/", gene_id)
  lines <- tryCatch(readLines(url, warn = FALSE), 
                    error = function(e) return(NA))
  
  # If lines failed, skip
  if (length(lines) == 1 && is.na(lines)) return(NA)
  
  # Find ORTHOLOGY line
  ko_line <- grep("^ORTHOLOGY", lines, value = TRUE)
  if (length(ko_line) > 0) {
    ko <- strsplit(ko_line, "\\s+")[[1]][2]
    return(ko)  # Fixed: returning ko_re instead of ko
  } else {
    return(NA)
  }
}

# Apply function to each gene ID with rate limit
results_ko <- pbsapply(KEGG_ids, function(x) {
  Sys.sleep(0.5)  # respect KEGG server limits
  getKO(x)
})
# Create KO mapping table
ko_mapping <- data.frame(
  KEGG = names(results_ko),
  KO = unname(results_ko),
  stringsAsFactors = FALSE
)
# Merge with background data
final_table_bg <- merge(background_gene_uniprot_kegg_clean2, 
                        ko_mapping, 
                        by="KEGG", 
                        all.x=TRUE)
# Save results
write.table(final_table_bg,
            "DE_Results_limma_voom_2025/real_annotation/DE_coastR/background_kegg_kos.tsv",
            sep="\t",
            quote=FALSE,
            row.names=FALSE)
##bulk mapping
#Step 1 — bulk KO mapping (NO loops)
# Split IDs into chunks (KEGG limit ≈ 100–200 per request)
chunks <- split(KEGG_ids, ceiling(seq_along(KEGG_ids)/100))

get_chunk <- function(ids){
  url <- paste0("https://rest.kegg.jp/link/ko/", paste(ids, collapse="+"))
  read.delim(url, header=FALSE, stringsAsFactors=FALSE)
}

ko_list <- lapply(chunks, get_chunk)

ko_df <- do.call(rbind, ko_list)
colnames(ko_df) <- c("KEGG", "KO")
head(ko_df)
# clean prefixes
ko_df$KEGG <- sub("genes:", "", ko_df$KEGG)
ko_df$KO   <- sub("ko:", "", ko_df$KO)
head(ko_df)
#Step 2 — merge back
background_gene_uniprot_kegg_clean2$KEGG_clean <- sub(";.*", "", background_gene_uniprot_kegg_clean2$KEGG)
protein_kegg_ko <- merge(background_gene_uniprot_kegg_clean2,
                         ko_df,
                         by.x = "KEGG_clean",
                         by.y = "KEGG",
                         all.x=TRUE)
head(protein_kegg_ko)
write.table(protein_kegg_ko,
            "DE_Results_limma_voom_2025/real_annotation/DE_coastR/background_kegg_kos.tsv",
            sep="\t",
            quote=FALSE,
            row.names=FALSE)
##KEGG enrichment with clusterProfiler
library(clusterProfiler)
background_KOs <- protein_kegg_ko %>%
  distinct(protein_id, KO)
sig_region_genes_real <- read.delim(
  "DE_Results_limma_voom_2025/real_annotation/DE_coastR/Region_Ftest_sig_genes.csv",
  sep = ",",
  header = TRUE,
  stringsAsFactors = FALSE
)
head(sig_region_genes_real)
sig_region_kos <- sig_region_genes_real %>%
  left_join(background_KOs, by = c("Gene" = "protein_id")) %>%
  filter(!is.na(KO))
head(sig_region_kos)
sig_region_work_KOs <- unique(sig_region_kos$KO)
background_work_KOs <- unique(background_KOs$KO)
kegg_region_enrich <- enrichKEGG(
  gene     = sig_region_work_KOs,
  universe = background_work_KOs,
  organism = "ko",
  keyType  = "kegg"
)
head(kegg_region_enrich)
write.table(as.data.frame(kegg_region_enrich),
            "DE_Results_limma_voom_2025/real_annotation/DE_coastR/kegg_region_enrichment.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)
##treatment##
sig_treat_genes_real <- read.delim(
  "DE_Results_limma_voom_2025/real_annotation/DE_coastR/Salinity_SWvsAW_sig_DE_genes.csv",
  sep = ",",
  header = TRUE,
  stringsAsFactors = FALSE
)
head(sig_treat_genes_real)
sig_treat_kos <- sig_treat_genes_real %>%
  left_join(background_KOs, by = c("Gene" = "protein_id")) %>%
  filter(!is.na(KO))
head(sig_treat_kos)
sig_treat_work_KOs <- unique(sig_treat_kos$KO)
kegg_treat_enrich <- enrichKEGG(
  gene     = sig_treat_work_KOs,
  universe = background_work_KOs,
  organism = "ko",
  keyType  = "kegg"
)
head(kegg_treat_enrich)
write.table(as.data.frame(kegg_treat_enrich),
            "DE_Results_limma_voom_2025/real_annotation/DE_coastR/kegg_treat_enrichment.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)
###interaction of treatment and habitat(region)###
sig_inter_genes_real <- read.delim(
  "DE_Results_limma_voom_2025/real_annotation/DE_coastR/Region_Treatment_interaction_sig_genes.csv",
  sep = ",",
  header = TRUE,
  stringsAsFactors = FALSE
)
head(sig_inter_genes_real)
sig_inter_kos <- sig_inter_genes_real %>%
  left_join(background_KOs, by = c("Gene" = "protein_id")) %>%
  filter(!is.na(KO))
head(sig_inter_kos)
sig_inter_work_KOs <- unique(sig_inter_kos$KO)
kegg_inter_enrich <- enrichKEGG(
  gene     = sig_inter_work_KOs,
  universe = background_work_KOs,
  organism = "ko",
  keyType  = "kegg"
)
head(kegg_inter_enrich)
write.table(as.data.frame(kegg_inter_enrich),
            "DE_Results_limma_voom_2025/real_annotation/DE_coastR/kegg_inter_enrichment.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)
##DEGs  across three regions under AW####
sig_AW_genes_real <- read.delim(
  "DE_Results_limma_voom_2025/real_annotation/DE_coastR/sig_region_AW_F.csv",
  sep = ",",
  header = TRUE,
  stringsAsFactors = FALSE
)
head(sig_AW_genes_real)
sig_AW_kos <- sig_AW_genes_real %>%
  left_join(background_KOs, by = c("Gene" = "protein_id")) %>%
  filter(!is.na(KO))
head(sig_AW_kos)
sig_AW_work_KOs <- unique(sig_AW_kos$KO)
kegg_AW_enrich <- enrichKEGG(
  gene     = sig_AW_work_KOs,
  universe = background_work_KOs,
  organism = "ko",
  keyType  = "kegg"
)
head(kegg_AW_enrich)
write.table(as.data.frame(kegg_AW_enrich),
            "DE_Results_limma_voom_2025/real_annotation/DE_coastR/kegg_AW_enrichment.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)
##DEGs  across three regions under SW####
sig_SW_genes_real <- read.delim(
  "DE_Results_limma_voom_2025/real_annotation/DE_coastR/sig_region_SW_F.csv",
  sep = ",",
  header = TRUE,
  stringsAsFactors = FALSE
)
head(sig_SW_genes_real)
sig_SW_kos <- sig_SW_genes_real %>%
  left_join(background_KOs, by = c("Gene" = "protein_id")) %>%
  filter(!is.na(KO))
head(sig_SW_kos)
sig_SW_work_KOs <- unique(sig_SW_kos$KO)
kegg_SW_enrich <- enrichKEGG(
  gene     = sig_SW_work_KOs,
  universe = background_work_KOs,
  organism = "ko",
  keyType  = "kegg"
)
head(kegg_SW_enrich)
write.table(as.data.frame(kegg_SW_enrich),
            "DE_Results_limma_voom_2025/real_annotation/DE_coastR/kegg_SW_enrichment.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)
###overlapped gene across three effects###
library(purrr)
overlapped_genes <- list(
  Treatment = sig_treat_genes_real,
  Habitat = sig_region_genes_real,
  Interaction = sig_inter_genes_real
) %>%
  # Extract the 'Gene' column from each data frame in the list
  map(~ .x$Gene) %>% 
  # Now intersect the resulting vectors
  reduce(intersect)
head(overlapped_genes)
overlap_df <- data.frame(Gene = overlapped_genes)
write.table(overlap_df, "DE_Results_limma_voom_2025/real_annotation/DE_coastR/overlapped_gene_across_effects.tsv", sep = "\t", quote = FALSE,row.names = FALSE,col.names = FALSE)
head(overlap_df)
overlap_kos <- overlap_df %>%
  left_join(background_KOs, by = c("Gene" = "protein_id")) %>%
  filter(!is.na(KO))
head(overlap_kos)
overlap_work_KOs <- unique(overlap_kos$KO)
kegg_overlap_enrich <- enrichKEGG(
  gene     = overlap_work_KOs,
  universe = background_work_KOs,
  organism = "ko",
  keyType  = "kegg"
)
head(kegg_overlap_enrich)
write.table(as.data.frame(kegg_overlap_enrich),
            "DE_Results_limma_voom_2025/real_annotation/DE_coastR/kegg_overlap_enrichment.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)
######################################