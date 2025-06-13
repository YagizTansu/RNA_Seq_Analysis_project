# Load required libraries
library(edgeR)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(stringr)

# Create results directory if it doesn't exist
dir.create("results", showWarnings = TRUE)

# Read metadata
if (!file.exists("data/Metadata2025.csv")) {
    stop("Metadata file not found")
}
metadata <- read.csv("data/Metadata2025.csv", sep="\t", header=TRUE)

# Read counts data
if (!file.exists("data/Counts2025.csv")) {
    stop("Counts file not found")
}
counts <- read.csv("data/Counts2025.csv", sep="\t", row.names=1)

# Create DGEList object
dge <- DGEList(counts=counts, 
               samples=metadata)

# Check the dimensions and structure
dim(dge)
head(dge$counts)
head(dge$samples)

# Filter low expressed genes (>10 reads in at least 1 sample)
keep <- rowSums(dge$counts > 10) >= 1
dge <- dge[keep, , keep.lib.sizes=FALSE]
dim(dge)

# Normalize the data
dge <- calcNormFactors(dge)
dge$samples$norm.factors

# Calculate CPM values (non-log scaled)
cpm_values <- cpm(dge, log=FALSE)

# Calculate average expression per individual
donor_averages <- list()
for (donor in unique(metadata$Donor)) {
    donor_samples <- metadata$Donor == donor
    donor_averages[[donor]] <- rowMeans(cpm_values[, donor_samples], na.rm = TRUE)
}

# Convert to data frame for easier comparison
avg_expr <- data.frame(
    S7 = donor_averages$S7,
    S12 = donor_averages$S12,
    S13 = donor_averages$S13
)
rownames(avg_expr) <- rownames(cpm_values)

# Classify genes based on expression patterns
classify_expression <- function(row) {
    max_expr <- max(row)
    others <- row[row != max_expr]
    avg_others <- mean(others)
    
    if (max_expr >= 4 * max(others)) {
        return("individual_specific")
    } else if (max_expr >= 2 * avg_others) {
        return("individual_elevated")
    } else {
        return("not_elevated")
    }
}

# Apply classification to each gene
expression_classes <- apply(avg_expr, 1, classify_expression)

# Compute statistics per individual
individual_stats <- data.frame(
    Individual = unique(metadata$Donor),
    Specific_Genes = 0,
    Elevated_Genes = 0
)

for (donor in unique(metadata$Donor)) {
    # Get genes where this individual has maximum expression
    max_expr <- apply(avg_expr, 1, function(row) which.max(row) == which(colnames(avg_expr) == donor))
    
    # Count specific and elevated genes for this individual
    individual_stats$Specific_Genes[individual_stats$Individual == donor] <- 
        sum(expression_classes == "individual_specific" & max_expr)
    
    individual_stats$Elevated_Genes[individual_stats$Individual == donor] <- 
        sum(expression_classes == "individual_elevated" & max_expr)
}

# Print individual statistics
print("Individual-specific and elevated gene counts:")
print(individual_stats)

# Save individual statistics
write.csv(individual_stats, "results/individual_gene_stats.csv", row.names = FALSE)

# Create summary table
expression_summary <- table(expression_classes)

# Create detailed results table
expression_results <- data.frame(
    Gene = rownames(avg_expr),
    S7_avg = avg_expr$S7,
    S12_avg = avg_expr$S12,
    S13_avg = avg_expr$S13,
    Class = expression_classes
)

# Save results
write.csv(expression_results, "results/individual_expression_patterns.csv")

# Print summary
print("Expression pattern summary:")
print(expression_summary)

# Create a bar plot of expression classes
class_plot <- ggplot(data.frame(Class = expression_classes)) +
    geom_bar(aes(x = Class, fill = Class)) +
    theme_bw() +
    ggtitle("Distribution of Expression Patterns") +
    xlab("Expression Class") +
    ylab("Number of Genes")

# Save the plot
ggsave("results/expression_classes_plot.png", class_plot, width=8, height=6)
# Write plot description
cat("Expression Classes Plot Description:\n",
    "This bar plot shows the distribution of genes across three expression categories:\n",
    "- Individual-specific: Genes with expression levels 4x higher in one individual\n",
    "- Individual-elevated: Genes with expression levels 2x higher in one individual\n",
    "- Not elevated: Genes with similar expression across individuals\n",
    file="results/expression_classes_plot_description.txt")

# Create PCA plot
# Get normalized log-CPM values
logcpm <- cpm(dge, log=TRUE)

# Perform PCA
pca <- prcomp(t(logcpm))
pca_data <- as.data.frame(pca$x)
pca_data$Tissue <- metadata$Tissue
pca_data$Sample <- metadata$Sample
pca_data$Donor <- metadata$Donor # Added Donor column

# Create PCA plot using ggplot2
pca_plot <- ggplot(pca_data, aes(x=PC1, y=PC2, color=Tissue, label=Donor)) + # Changed label to Donor
  geom_point(size=3) +
  geom_text(hjust=0.5, vjust=-0.5) +
  theme_bw() +
  xlab(paste0("PC1 (", round(summary(pca)$importance[2,1]*100, 1), "%)")) +
  ylab(paste0("PC2 (", round(summary(pca)$importance[2,2]*100, 1), "%)")) +
  ggtitle("PCA Plot by Tissue")

# Save the plot
ggsave("results/PCA_plot.png", pca_plot, width=10, height=8)
# Write PCA plot description
cat("PCA Plot Description:\n",
    "Principal Component Analysis (PCA) plot showing sample clustering by tissue type.\n",
    "- Points are colored by tissue type\n",
    "- Labels show donor IDs (S7, S12, S13)\n",
    "- PC1 and PC2 represent the two main axes of variation in the data\n",
    "- Samples clustering together indicate similar expression profiles\n",
    file="results/PCA_plot_description.txt")

# Create MDS plot
# Calculate MDS
mds <- plotMDS(dge, plot=FALSE)
mds_data <- data.frame(
    Dim1 = mds$x,
    Dim2 = mds$y,
    Tissue = metadata$Tissue,
    Donor = metadata$Donor
)

# Create MDS plot using ggplot2
mds_plot <- ggplot(mds_data, aes(x=Dim1, y=Dim2, color=Tissue, label=Donor)) +
    geom_point(size=3) +
    geom_text(hjust=0.5, vjust=-0.5) +
    theme_bw() +
    xlab("Leading LogFC Dim 1") +
    ylab("Leading LogFC Dim 2") +
    ggtitle("MDS Plot by Tissue")

# Save the MDS plot
ggsave("results/MDS_plot.png", mds_plot, width=10, height=8)

# Write MDS plot description
cat("MDS Plot Description:\n",
    "Multi-dimensional Scaling (MDS) plot showing sample relationships based on expression profiles.\n",
    "- Points are colored by tissue type\n",
    "- Labels show donor IDs (S7, S12, S13)\n",
    "- Distances between samples represent leading fold-change differences\n",
    "- Similar samples cluster together\n",
    file="results/MDS_plot_description.txt")

# Differential Expression Analysis: Liver vs Brain
# Create a new group factor for the selected tissues
selected_tissues <- c("liver", "brain")
selected_samples <- dge$samples$Tissue %in% selected_tissues
dge_selected <- dge[, selected_samples]
dge_selected$samples$group <- factor(dge_selected$samples$Tissue)

# Design matrix
design <- model.matrix(~group, data=dge_selected$samples)

# Estimate dispersion
dge_selected <- estimateCommonDisp(dge_selected)
dge_selected <- estimateTagwiseDisp(dge_selected)

# Perform exact test for differential expression
et <- exactTest(dge_selected, pair=c("liver", "brain"))

# Get top differentially expressed genes
top_genes <- topTags(et, n=Inf)
sig_genes <- top_genes$table[top_genes$table$FDR < 0.05,]

# Create and filter topTags object with FDR <= 0.01
top_genes_stringent <- topTags(et, n=Inf)
deg_list_stringent <- top_genes_stringent$table[top_genes_stringent$table$FDR <= 0.01,]

# Create a summary of results
de_summary <- data.frame(
  Total_DE_Genes = nrow(sig_genes),
  Upregulated = sum(sig_genes$logFC > 0),
  Downregulated = sum(sig_genes$logFC < 0)
)

# Create a more detailed summary including the stringent threshold
de_summary_stringent <- data.frame(
  FDR_Threshold = c(0.05, 0.01),
  Total_DE_Genes = c(nrow(sig_genes), nrow(deg_list_stringent)),
  Upregulated = c(sum(sig_genes$logFC > 0), sum(deg_list_stringent$logFC > 0)),
  Downregulated = c(sum(sig_genes$logFC < 0), sum(deg_list_stringent$logFC < 0))
)

# Save the stringent DEG list
write.csv(deg_list_stringent, "results/liver_vs_brain_DE_genes_FDR01.csv")

# Create volcano plot
volcano_plot <- ggplot(top_genes$table, aes(x=logFC, y=-log10(FDR))) +
  geom_point(aes(color=FDR < 0.05), size=1) +
  scale_color_manual(values=c("grey", "red")) +
  theme_bw() +
  ggtitle("Volcano Plot: Liver vs Brain") +
  xlab("Log2 Fold Change") +
  ylab("-log10 FDR") +
  theme(legend.position="none")

# Save results
write.csv(sig_genes, "results/liver_vs_brain_DE_genes.csv")
ggsave("results/volcano_plot.png", volcano_plot, width=8, height=6)
# Write volcano plot description
cat("Volcano Plot Description:\n",
    "Volcano plot showing differential expression between liver and brain samples.\n",
    "- X-axis: Log2 fold change in expression\n",
    "- Y-axis: Statistical significance (-log10 FDR)\n",
    "- Red points: Significantly differentially expressed genes (FDR < 0.05)\n",
    "- Grey points: Non-significant genes\n",
    file="results/volcano_plot_description.txt")

# Classify genes into 4 categories
gene_classes <- data.frame(top_genes$table) %>%
  mutate(class = case_when(
    FDR <= 0.01 & logFC > 0 ~ "DE_UP",
    FDR <= 0.01 & logFC < 0 ~ "DE_DOWN",
    FDR > 0.01 & logFC > 0 ~ "notDE_UP",
    FDR > 0.01 & logFC < 0 ~ "notDE_DOWN"
  ))

# Convert class to factor with specific order
gene_classes$class <- factor(gene_classes$class, 
                           levels = c("DE_UP", "DE_DOWN", "notDE_UP", "notDE_DOWN"))

# Create boxplot
class_boxplot <- ggplot(gene_classes, aes(x=class, y=logFC, fill=class)) +
  geom_boxplot() +
  theme_bw() +
  scale_fill_brewer(palette="Set2") +
  ggtitle("Distribution of logFC by Gene Class") +
  xlab("Gene Class") +
  ylab("Log2 Fold Change") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the boxplot
ggsave("results/gene_class_boxplot.png", class_boxplot, width=8, height=6)
# Write boxplot description
cat("Gene Class Boxplot Description:\n",
    "Distribution of log2 fold changes across different gene categories:\n",
    "- DE_UP: Significantly upregulated genes (FDR ≤ 0.01)\n",
    "- DE_DOWN: Significantly downregulated genes (FDR ≤ 0.01)\n",
    "- notDE_UP: Non-significant genes with positive log fold change\n",
    "- notDE_DOWN: Non-significant genes with negative log fold change\n",
    file="results/gene_class_boxplot_description.txt")

# Print summary of gene counts in each class
class_summary <- gene_classes %>%
  group_by(class) %>%
  summarise(count = n())
print("Number of genes in each class:")
print(class_summary)

# Print summary
print("Differential Expression Analysis Summary:")
print(de_summary)

# Print updated summary
print("Differential Expression Analysis Summary at different FDR thresholds:")
print(de_summary_stringent)

# Analyze overlap between DEGs and individual patterns
deg_genes <- rownames(sig_genes)
individual_specific_genes <- rownames(expression_results)[expression_results$Class == "individual_specific"]
individual_elevated_genes <- rownames(expression_results)[expression_results$Class == "individual_elevated"]

# Calculate overlaps
deg_specific_overlap <- intersect(deg_genes, individual_specific_genes)
deg_elevated_overlap <- intersect(deg_genes, individual_elevated_genes)

# Create summary
deg_individual_summary <- data.frame(
    Category = c("Individual-specific DEGs", "Individual-elevated DEGs"),
    Count = c(length(deg_specific_overlap), length(deg_elevated_overlap)),
    Percentage_of_DEGs = c(
        length(deg_specific_overlap) / length(deg_genes) * 100,
        length(deg_elevated_overlap) / length(deg_genes) * 100
    )
)

# Print results
print("\nOverlap between DEGs and Individual Expression Patterns:")
print(deg_individual_summary)

# Save results
write.csv(deg_individual_summary, "results/deg_individual_overlap_summary.csv", row.names = FALSE)

# Function to analyze gene expression patterns
analyze_gene <- function(gene_name, dge, metadata, expression_results) {
    # Check if gene exists in the data
    if (!gene_name %in% rownames(dge$counts)) {
        stop("Gene not found in dataset")
    }
    
    # Get CPM values for the gene
    cpm_values <- cpm(dge, log=FALSE)[gene_name,]
    
    # Calculate average expression per individual
    avg_by_individual <- tapply(cpm_values, metadata$Donor, mean)
    
    # Get gene classification information
    gene_info <- expression_results[expression_results$Gene == gene_name,]
    
    # Identify individuals where gene is elevated or specific
    elevated_in <- NULL
    specific_in <- NULL
    
    if (gene_info$Class == "individual_elevated") {
        elevated_in <- names(which.max(c(gene_info$S7_avg, gene_info$S12_avg, gene_info$S13_avg)))
    } else if (gene_info$Class == "individual_specific") {
        specific_in <- names(which.max(c(gene_info$S7_avg, gene_info$S12_avg, gene_info$S13_avg)))
    }
    
    # Create expression plot
    plot_data <- data.frame(
        Expression = cpm_values,
        Individual = metadata$Donor,
        Tissue = metadata$Tissue
    )
    
    expression_plot <- ggplot(plot_data, aes(x=Individual, y=Expression)) +
        geom_boxplot(aes(fill=Individual)) +
        geom_point(aes(color=Tissue)) +
        theme_bw() +
        ggtitle(paste("Expression of", gene_name)) +
        ylab("CPM") +
        theme(legend.position="right")
    
    # Print results
    cat("\nResults for", gene_name, ":\n")
    cat("\nAverage expression by individual:\n")
    print(avg_by_individual)
    
    if (!is.null(elevated_in)) {
        cat("\nGene is elevated in:", elevated_in, "\n")
    }
    
    if (!is.null(specific_in)) {
        cat("\nGene is specific to:", specific_in, "\n")
    }
    
    # Explicitly print the plot
    print(expression_plot)
    
    # Save the plot to a file
    ggsave(
        filename = paste0("results/", gene_name, "_expression.png"),
        plot = expression_plot,
        width = 8,
        height = 6
    )
    # Write gene-specific plot description
    cat(paste0("Gene Expression Plot Description for ", gene_name, ":\n",
        "Boxplot showing expression levels across different individuals and tissues:\n",
        "- Y-axis: Expression level in CPM (Counts Per Million)\n",
        "- X-axis: Individual donors (S7, S12, S13)\n",
        "- Points: Individual tissue samples\n",
        "- Colors: Different tissue types\n"),
        file=paste0("results/", gene_name, "_expression_description.txt"))
    
    # Return results as a list
    return(list(
        averages = avg_by_individual,
        elevated_in = elevated_in,
        specific_in = specific_in,
        plot = expression_plot
    ))
}

# Example usage - make sure to assign the result AND print the plot:
# result <- analyze_gene("gene_name", dge, metadata, expression_results)

result <- analyze_gene("AL157440", dge, metadata, expression_results)
print(result$plot)  # Explicitly print the plot
