# =========================================
# Venn Diagram Analysis (Gene Overlap)
# =========================================

# ===== LOAD LIBRARIES =====
if (!requireNamespace("VennDiagram", quietly = TRUE)) install.packages("VennDiagram")
if (!requireNamespace("gridExtra", quietly = TRUE)) install.packages("gridExtra")

library(VennDiagram)
library(grid)

# ===== CONFIG =====
base_dir <- "~/transkriptomik/venn"
output_dir <- "~/transkriptomik/results_venn"
dir.create(output_dir, showWarnings = FALSE)

# ===== LOAD DATA =====
n22_genes <- readLines(file.path(base_dir, "N22_Fix.txt"))
pokkali_genes <- readLines(file.path(base_dir, "Pokkali_Fix.txt"))
nac_genes <- readLines(file.path(base_dir, "NAC.txt"))

# ===== CALCULATE INTERSECTIONS =====
n22_pok <- intersect(n22_genes, pokkali_genes)
n22_nac <- intersect(n22_genes, nac_genes)
pok_nac <- intersect(pokkali_genes, nac_genes)
all_three <- Reduce(intersect, list(n22_genes, pokkali_genes, nac_genes))

# ===== SAVE INTERSECTION RESULTS =====
writeLines(n22_pok, file.path(output_dir, "N22_Pokkali_overlap.txt"))
writeLines(n22_nac, file.path(output_dir, "N22_NAC_overlap.txt"))
writeLines(pok_nac, file.path(output_dir, "Pokkali_NAC_overlap.txt"))
writeLines(all_three, file.path(output_dir, "All_overlap.txt"))

# ===== CREATE VENN DIAGRAM =====
venn.plot <- draw.triple.venn(
  area1 = length(n22_genes),
  area2 = length(pokkali_genes),
  area3 = length(nac_genes),
  n12 = length(n22_pok),
  n13 = length(n22_nac),
  n23 = length(pok_nac),
  n123 = length(all_three),
  category = c("N22", "Pokkali", "NAC"),
  fill = c("red", "green", "blue"),
  cex = 1.2,
  cat.cex = 1.2,
  cat.pos = c(0, 0, 180),
  cat.dist = c(0.05, 0.05, 0.05)
)

# ===== SAVE PLOT =====
png(file.path(output_dir, "Venn_Diagram.png"), width = 800, height = 800)
grid.draw(venn.plot)
dev.off()

# ===== SUMMARY =====
cat("==== VENN SUMMARY ====\n")
cat("N22 genes:", length(n22_genes), "\n")
cat("Pokkali genes:", length(pokkali_genes), "\n")
cat("NAC genes:", length(nac_genes), "\n\n")

cat("N22 ∩ Pokkali:", length(n22_pok), "\n")
cat("N22 ∩ NAC:", length(n22_nac), "\n")
cat("Pokkali ∩ NAC:", length(pok_nac), "\n")
cat("All three:", length(all_three), "\n")
