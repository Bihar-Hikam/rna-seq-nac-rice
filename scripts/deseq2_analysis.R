# =========================================
# RNA-Seq DEG Analysis Pipeline (DESeq2)
# =========================================

# ===== LOAD LIBRARIES =====
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

packages <- c("DESeq2", "tidyverse", "pheatmap", "EnhancedVolcano")
for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% c("DESeq2", "EnhancedVolcano")) {
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
  }
  library(pkg, character.only = TRUE)
}

# ===== CONFIG =====
base_dir <- "~/transkriptomik/counts"
output_dir <- "~/transkriptomik/results_deseq2"
dir.create(output_dir, showWarnings = FALSE)

# ===== FUNCTION: RUN DESEQ2 =====
run_deseq2 <- function(sample_names, conditions, prefix) {
  
  # Create metadata
  sample_info <- data.frame(
    sampleName = sample_names,
    condition = conditions
  )
  rownames(sample_info) <- sample_info$sampleName
  sample_info$condition <- as.factor(sample_info$condition)
  
  # Read count files
  count_list <- lapply(sample_names, function(file) {
    df <- read.table(file.path(base_dir, file), header = FALSE, stringsAsFactors = FALSE)
    df <- df[!grepl("^__", df$V1), ]
    colnames(df) <- c("gene_id", file)
    return(df)
  })
  
  # Merge all counts
  count_data <- Reduce(function(x, y) merge(x, y, by = "gene_id"), count_list)
  
  # Set rownames
  rownames(count_data) <- count_data$gene_id
  count_data <- count_data[, -1]
  
  # Ensure order match
  sample_info <- sample_info[colnames(count_data), ]
  
  # ===== DESEQ2 =====
  dds <- DESeqDataSetFromMatrix(
    countData = count_data,
    colData = sample_info,
    design = ~ condition
  )
  
  dds <- DESeq(dds)
  res <- results(dds)
  resOrdered <- res[order(res$padj), ]
  
  # Save results
  write.csv(as.data.frame(resOrdered),
            file.path(output_dir, paste0("DEG_", prefix, ".csv")))
  
  # ===== VISUALIZATION =====
  
  # PCA
  vsd <- vst(dds)
  png(file.path(output_dir, paste0("PCA_", prefix, ".png")))
  plotPCA(vsd, intgroup = "condition")
  dev.off()
  
  # Volcano
  res_df <- as.data.frame(res)
  res_df$gene <- rownames(res_df)
  
  res_df$pvalue[res_df$pvalue == 0] <- 1e-300
  
  res_df <- res_df %>%
    mutate(threshold = case_when(
      pvalue < 0.05 & log2FoldChange > 0 ~ "Upregulated",
      pvalue < 0.05 & log2FoldChange < 0 ~ "Downregulated",
      TRUE ~ "Not Significant"
    ))
  
  p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = threshold)) +
    geom_point(alpha = 0.8, size = 2) +
    scale_color_manual(values = c("Upregulated" = "red",
                                  "Downregulated" = "blue",
                                  "Not Significant" = "grey")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    theme_minimal() +
    labs(title = paste("Volcano Plot -", prefix))
  
  ggsave(file.path(output_dir, paste0("Volcano_", prefix, ".png")), p)
  
  # Filter DEG
  deg_filtered <- res_df %>% filter(pvalue < 0.05)
  
  write.csv(deg_filtered,
            file.path(output_dir, paste0("DEG_filtered_", prefix, ".csv")),
            row.names = FALSE)
  
  # Stats
  cat("====", prefix, "====\n")
  cat("Total DEG:", nrow(deg_filtered), "\n")
  cat("Upregulated:", sum(deg_filtered$log2FoldChange > 0), "\n")
  cat("Downregulated:", sum(deg_filtered$log2FoldChange < 0), "\n\n")
}

# =========================================
# RUN ANALYSIS
# =========================================

# ===== N22 (Drought vs Control) =====
run_deseq2(
  sample_names = c(
    "N22_R1_Drought_counts.txt",
    "N22_R1_Kontrol_counts.txt",
    "N22_R2_Drought_counts.txt",
    "N22_R2_Kontrol_counts.txt",
    "N22_R3_Drought_counts.txt",
    "N22_R3_Kontrol_counts.txt"
  ),
  conditions = c("Drought", "Control", "Drought", "Control", "Drought", "Control"),
  prefix = "N22"
)

# ===== Pokkali (Salt vs Control) =====
run_deseq2(
  sample_names = c(
    "Pok_R1_Salt_counts.txt",
    "Pok_R1_Kontrol_counts.txt",
    "Pok_R2_Salt_counts.txt",
    "Pok_R2_Kontrol_counts.txt",
    "Pok_R3_Salt_counts.txt",
    "Pok_R3_Kontrol_counts.txt"
  ),
  conditions = c("Salt", "Control", "Salt", "Control", "Salt", "Control"),
  prefix = "Pokkali"
)
