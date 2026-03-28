#!/usr/bin/env Rscript
# ============================================================================
# Unified Analysis Configuration for FAdV-4 Transcriptome Meta-Analysis
# ============================================================================
# Version: v2.0
# Author: Expert R Development Team
# Date: 2026-02-03
#
# Description: Master configuration file for the entire analysis pipeline
# Usage: source("config.R")  # Load at the beginning of each Part script
# ============================================================================

# ============================================================================
# Section 1: Analysis Parameters / 分析参数设置
# ============================================================================

# Significance thresholds (based on literature: Zhang et al. 2018, Li et al. 2025)
ANALYSIS_PARAMS <- list(
  # FDR threshold - unified 0.05 for all analyses / FDR阈值
  fdr_threshold = 0.05,

  # Log2 fold change threshold / log2折叠变化阈值
  # Updated from 1.0 to 1.5 based on power analysis (Part1, ~70-75% power)
  logfc_threshold = 1.5,

  # Effect size consistency threshold / 效应量一致性阈值
  # effect_ratio > 0.30 means effect sizes differ by < 70% across datasets
  effect_ratio_threshold = 0.30,

  # Minimum effect size requirement / 最小效应量要求
  min_effect_size = 0.5
)

# ============================================================================
# Section 2: Data Filtering Parameters / 数据过滤参数
# ============================================================================

FILTER_PARAMS <- list(
  # CPM threshold - gene expression filtering / CPM阈值
  cpm_threshold = 1.0,

  # Minimum sample count / 最小样本数
  min_samples = 3,

  # Low expression gene filtering threshold / 低表达基因过滤阈值
  low_count_threshold = 10,

  # Sample filtering minimum reads / 样本过滤最小reads数
  min_library_size = 1000000  # 1 million reads
)

# ============================================================================
# Section 3: Enrichment Analysis Parameters / 富集分析参数
# ============================================================================

ENRICHMENT_PARAMS <- list(
  # GO enrichment parameters / GO富集分析参数
  go = list(
    pvalue_cutoff = 0.05,
    qvalue_cutoff = 0.05,
    min_gssize = 10,
    max_gssize = 500,      # Modified to 500 to capture large GO terms
    p_adjust_method = "BH"
  ),

  # KEGG pathway parameters / KEGG通路分析参数
  kegg = list(
    pvalue_cutoff = 0.05,
    qvalue_cutoff = 0.05,
    min_gssize = 10,
    max_gssize = 500,      # Modified to 500 to capture large pathways
    p_adjust_method = "BH",
    organism = "gga"  # Gallus gallus (chicken)
  ),

  # GSEA parameters / GSEA参数
  gsea = list(
    n_perm = 1000,           # Permutation次数
    min_gssize = 10,
    max_gssize = 500,
    pvalue_cutoff = 0.05
  )
)

# ============================================================================
# Section 4: Leave-One-Out Validation Parameters / 留一法验证参数
# ============================================================================

VALIDATION_PARAMS <- list(
  # Leave-one-out thresholds / 留一法验证阈值
  # NOTE: Current implementation uses same thresholds as main analysis (FDR=0.05)
  # These parameters are reserved for future alternative validation strategies
  # 当前实现使用与主分析相同的阈值（FDR=0.05）
  # 这些参数保留用于未来的替代验证策略
  loo_fdr_threshold = 0.10,   # Reserved: Not currently used
  loo_logfc_threshold = 1.5,  # Reserved: Not currently used (updated for consistency)
  loo_effect_ratio = 0.40     # Reserved: Not currently used
)

# ============================================================================
# Section 5: Visualization Parameters / 可视化参数
# ============================================================================

VIZ_PARAMS <- list(
  # Figure dimensions / 图表尺寸
  figure_width = 10,
  figure_height = 8,
  figure_dpi = 300,

  # Color schemes / 颜色方案
  colors = list(
    control = "#2E86AB",      # Blue
    favd4 = "#A23B72",        # Magenta
    upregulated = "#D64045",  # Red
    downregulated = "#1B998B" # Green
  ),

  # Dot plot parameters / 点图参数
  dotplot = list(
    show_category = 20,       # Number of items to display
    font_size = 12
  )
)

# ============================================================================
# Section 6: Path Configuration / 路径配置
# ============================================================================

PATHS <- list(
  # Main directories / 主要目录
  # UPDATED: Unified results directory structure
  intermediate = "results/intermediate",
  tables = "results/tables",
  figures = "results/figures",
  enrichment = "results/tables/enrichment",
  ppi = "results/tables/ppi",
  qc = "results/figures/qc",

  # Input data paths (relative to project root) / 输入数据路径
  source_data = "source_data/processed_data",

  # Key filenames / 关键文件名
  files = list(
    # Intermediate directory
    expr_stats = "Normalized_Expression_For_Statistics.rds",
    expr_vis = "Normalized_Expression_For_Visualization.rds",

    # Tables directory
    core_degs = "Integrated_Core_DEGs.csv",
    exploratory_degs = "Integrated_Exploratory_DEGs.csv",
    rankagg_results = "Integrated_RankAggregation_Results.csv",
    deg_106839 = "GSE106839_DEGs.csv",
    deg_299945 = "GSE299945_DEGs.csv",
    validation_overlap = "Validation_Overlap_CoreDEGs.csv",
    hcs_results = "HCS_Scoring_Results.csv",
    hcs_core = "HCS_Core_DEGs.csv",

    # Metadata / 元数据
    sample_metadata = "source_data/sample_metadata.csv"
  )
)

# ============================================================================
# Section 7: QC and Validation Parameters / 质量控制和验证参数
# ============================================================================

QC_PARAMS <- list(
  # Gene ID conversion success rate requirement / 基因ID转换成功率要求
  min_id_conversion_rate = 0.75,  # At least 75% conversion rate

  # Batch effect detection / 批次效应检测
  batch_detection_method = "PCA",  # PCA or UMAP

  # Outlier sample detection / 离群样本检测
  outlier_threshold = 3,  # 3 standard deviations

  # Sample correlation check / 样本相关性检查
  min_correlation = 0.7   # Minimum sample correlation coefficient
)

# ============================================================================
# Section 8: Computational Parameters / 计算参数
# ============================================================================

COMPUTE_PARAMS <- list(
  # Parallel computing / 并行计算
  n_cores = parallel::detectCores() - 1,  # Keep one core free

  # Memory management / 内存管理
  max_memory_gb = 16,  # Maximum memory usage (GB)

  # Random seed / 随机种子
  random_seed = 12345,

  # Numerical stability / 数值稳定性
  min_pvalue = 1e-300,  # Minimum p-value (avoid log(0))
  epsilon = 1e-10       # Small constant to avoid division by zero
)

# ============================================================================
# Section 9: Logging and Reporting Parameters / 日志和报告参数
# ============================================================================

LOG_PARAMS <- list(
  # Log level: "DEBUG", "INFO", "WARNING", "ERROR" / 日志级别
  log_level = "INFO",

  # Save detailed logs / 是否保存详细日志
  verbose = TRUE,

  # Save intermediate results / 是否保存中间结果
  save_intermediate = TRUE,

  # Report format / 报告格式
  report_format = "markdown"  # or "html", "pdf"
)

# ============================================================================
# Section 10: Helper Functions / 辅助函数
# ============================================================================

#' Validate configuration completeness / 验证配置完整性
#' @return Logical indicating if configuration is valid
validate_config <- function() {
  required_sections <- c("ANALYSIS_PARAMS", "FILTER_PARAMS", "ENRICHMENT_PARAMS",
                        "VALIDATION_PARAMS", "PATHS")

  missing_sections <- c()
  for (section in required_sections) {
    if (!exists(section)) {
      missing_sections <- c(missing_sections, section)
    }
  }

  if (length(missing_sections) > 0) {
    warning(sprintf("Missing config sections: %s",
                   paste(missing_sections, collapse = ", ")))
    return(FALSE)
  }

  message("OK Configuration validated")
  return(TRUE)
}

#' Print configuration summary / 打印配置摘要
print_config_summary <- function() {
  cat("\n==============================================\n")
  cat("Analysis Configuration Summary\n")
  cat("==============================================\n\n")

  cat("Significance Thresholds:\n")
  cat(sprintf("  - FDR: %.2f\n", ANALYSIS_PARAMS$fdr_threshold))
  cat(sprintf("  - logFC: %.1f\n", ANALYSIS_PARAMS$logfc_threshold))
  cat(sprintf("  - Effect Ratio: %.2f\n", ANALYSIS_PARAMS$effect_ratio_threshold))

  cat("\nData Filtering:\n")
  cat(sprintf("  - CPM > %.1f\n", FILTER_PARAMS$cpm_threshold))
  cat(sprintf("  - Min samples: %d\n", FILTER_PARAMS$min_samples))

  cat("\nEnrichment Analysis:\n")
  cat(sprintf("  - GO maxGS: %d\n", ENRICHMENT_PARAMS$go$max_gssize))
  cat(sprintf("  - KEGG maxGS: %d\n", ENRICHMENT_PARAMS$kegg$max_gssize))

  cat("\nValidation:\n")
  cat(sprintf("  - LOO FDR: %.2f\n", VALIDATION_PARAMS$loo_fdr_threshold))

  cat("\n==============================================\n\n")
}

#' Get parameters convenience function / 获取参数的便捷函数
#' @param section Config section name / 配置部分名称
#' @param param Parameter name (optional) / 参数名称
#' @return Parameter value or list / 参数值或参数列表
get_param <- function(section, param = NULL) {
  if (!exists(section)) {
    stop(sprintf("Config section not found: %s", section))
  }

  section_obj <- get(section)

  if (is.null(param)) {
    return(section_obj)
  }

  if (!param %in% names(section_obj)) {
    stop(sprintf("Parameter not found: %s$%s", section, param))
  }

  return(section_obj[[param]])
}

#' Build file paths convenience function / 构建文件路径的便捷函数
#' @param file_type File type (from PATHS$files) / 文件类型
#' @param base_dir Base directory (optional) / 基础目录
#' @return Complete file path / 完整文件路径
get_path <- function(file_type, base_dir = NULL) {
  if (!file_type %in% names(PATHS$files)) {
    stop(sprintf("Unknown file type: %s", file_type))
  }

  filename <- PATHS$files[[file_type]]

  # Determine directory based on file type / 根据文件类型确定目录
  if (grepl("^expr_", file_type)) {
    dir <- PATHS$intermediate
  } else if (grepl("^deg_", file_type) || file_type %in%
             c("core_degs", "exploratory_degs", "rankagg_results",
               "validation_overlap", "hcs_results", "hcs_core")) {
    dir <- PATHS$tables
  } else if (file_type == "sample_metadata") {
    return(filename)  # Metadata has full path
  } else {
    dir <- PATHS$tables  # Default
  }

  if (is.null(base_dir)) {
    base_dir <- getwd()
  }

  file.path(base_dir, dir, filename)
}

# ============================================================================
# Initialization / 初始化
# ============================================================================

# Auto-validate configuration / 自动验证配置
if (interactive()) {
  validate_config()
  print_config_summary()
}

# Export main variables for other scripts / 导出主要变量
CONFIG_VERSION <- "2.0"
CONFIG_DATE <- "2026-02-03"

cat(sprintf("OK Configuration loaded (v%s, %s)\n", CONFIG_VERSION, CONFIG_DATE))
