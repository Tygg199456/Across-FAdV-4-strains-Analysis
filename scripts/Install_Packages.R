# ============================================================================
# Install Required R Packages for FAdV-4 Analysis
# ============================================================================

cat("\n==============================================================================\n")
cat("Installing Required R Packages for FAdV-4 Analysis\n")
cat("==============================================================================\n\n")

# Install CRAN packages / 安装CRAN包
cat("[Step 1/3] Installing CRAN packages\n")
cat("==============================================\n")

cran_packages <- c(
  "dplyr",
  "tidyr",
  "ggplot2",
  "ggrepel",
  "pheatmap",
  "gridExtra",
  "scales",
  "RColorBrewer"
)

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("Installing %s...\n", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org/")
    cat(sprintf("OK %s installed\n", pkg))
  } else {
    cat(sprintf("OK %s already installed\n", pkg))
  }
}

cat("\n")

# Install Bioconductor packages / 安装Bioconductor包
cat("[Step 2/3] Installing Bioconductor packages\n")
cat("==============================================\n")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  cat("Installing BiocManager...\n")
  install.packages("BiocManager", repos = "https://cloud.r-project.org/")
  cat("OK BiocManager installed\n\n")
}

bioc_packages <- c(
  "limma",
  "edgeR",
  "clusterProfiler",
  "org.Gg.eg.db",
  "enrichplot",
  "AnnotationDbi"
)

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("Installing %s...\n", pkg))
    BiocManager::install(pkg, update = FALSE)
    cat(sprintf("OK %s installed\n", pkg))
  } else {
    cat(sprintf("OK %s already installed\n", pkg))
  }
}

cat("\n")

# Verify installation / 验证安装
cat("[Step 3/3] Verifying Package Installation\n")
cat("==============================================\n")

all_packages <- c(cran_packages, bioc_packages)
failed_packages <- c()

for (pkg in all_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    failed_packages <- c(failed_packages, pkg)
    cat(sprintf("FAILED %s\n", pkg))
  } else {
    cat(sprintf("OK %s\n", pkg))
  }
}

cat("\n")
if (length(failed_packages) == 0) {
  cat("==============================================================================\n")
  cat("All packages installed successfully!\n")
  cat("==============================================================================\n\n")
} else {
  cat("==============================================================================\n")
  cat("WARNING: Failed to install:\n")
  cat("==============================================\n")
  for (pkg in failed_packages) {
    cat(sprintf("  - %s\n", pkg))
  }
  cat("\nPlease install these packages manually\n")
  cat("==============================================================================\n")
  quit(save = "no", status = 1)
}
