# STRING Network - Ranked by STRING Degree (TOP150 Only)

library(igraph)
library(ggraph)
library(ggplot2)
library(dplyr)
library(scales)

# Set working directory
setwd("/Users/tgw/Desktop/FADV_new/FAdV4_LimmaAnalysis")

# Output directories
figures_dir <- "results/figures"
tables_dir <- "results/tables"
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

cat("========================================\n")
cat("STRING Network - TOP150 Only, TOP100 Labeled\n")
cat("========================================\n\n")

# ============================================================================
# 1. Load Data
# ============================================================================
cat("[1/6] Loading data...\n")

# Read STRING interactions
string_data <- read.table("data/processed/string_interactions.tsv",
                          header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE,
                          comment.char = "", quote = "")

# Remove # prefix and fix X. prefix from column names
colnames(string_data) <- gsub("^#", "", colnames(string_data))
colnames(string_data) <- gsub("^X\\.", "", colnames(string_data))

# Read Consistent DEG data (489 high-confidence DEGs) - consistent with manuscript
# 标准：在两个数据集中都显著（FDR<0.05）+ 方向一致 + |avg_logFC|>1.5
deg_data <- read.csv("results/tables/Consistent_DEGs.csv",
                     header = TRUE, stringsAsFactors = FALSE,
                     check.names = FALSE)

cat("  - Total STRING interactions:", nrow(string_data), "\n")
cat("  - High-confidence DEG genes (489):", nrow(deg_data), "\n")

# ============================================================================
# 2. Calculate Degree and Rank by STRING Degree
# ============================================================================
cat("\n[2/6] Calculating degree and ranking...\n")

# Calculate degree for all genes in STRING
all_genes <- unique(c(string_data$node1, string_data$node2))
degree_table <- table(c(string_data$node1, string_data$node2))
degree_df <- data.frame(
  name = names(degree_table),
  degree = as.numeric(degree_table),
  stringsAsFactors = FALSE
) %>%
  arrange(desc(degree))

# Assign degree rank
degree_df <- degree_df %>%
  mutate(
    degree_rank = row_number(),
    tier = case_when(
      degree_rank <= 10 ~ "Top1-10",
      degree_rank <= 30 ~ "Top11-30",
      degree_rank <= 100 ~ "Top31-100",
      TRUE ~ "Top100+"
    ),
    # Shape mapping: all solid shapes
    shape = case_when(
      tier == "Top1-10" ~ 25,      # Triangle (solid)
      tier == "Top11-30" ~ 22,     # Square (solid)
      tier == "Top31-100" ~ 16,    # Circle (solid)
      tier == "Top100+" ~ 16,      # Circle (solid)
      TRUE ~ 16
    )
  )

cat("  - Total genes in STRING:", nrow(degree_df), "\n")
cat("  - degree >= 3:", sum(degree_df$degree >= 3), "\n")
cat("  - Top 1-10:", sum(degree_df$tier == "Top1-10"), "\n")
cat("  - Top 11-30:", sum(degree_df$tier == "Top11-30"), "\n")
cat("  - Top 31-100:", sum(degree_df$tier == "Top31-100"), "\n")
cat("  - Top 100+:", sum(degree_df$tier == "Top100+"), "\n")

# Filter to TOP150 by degree
degree_df_filtered <- degree_df %>%
  head(150)

cat("  - After filtering (TOP150):", nrow(degree_df_filtered), "nodes\n")

# ============================================================================
# 3. Build Network
# ============================================================================
cat("\n[3/6] Building network...\n")

# Filter edges to only include TOP150 genes
string_filtered <- string_data %>%
  filter(node1 %in% degree_df_filtered$name &
           node2 %in% degree_df_filtered$name)

cat("  - Filtered edges:", nrow(string_filtered), "\n")

# Merge with DEG data for direction and logFC
degree_df_filtered <- degree_df_filtered %>%
  left_join(deg_data %>% select(gene, avg_logFC) %>% rename(name = gene), by = "name") %>%
  mutate(
    avg_logFC = ifelse(is.na(avg_logFC), 0, avg_logFC),
    direction = ifelse(avg_logFC > 0, "Up", "Down")
  )

# Build network
protein_list <- degree_df_filtered %>%
  select(name, degree, degree_rank, tier, shape, avg_logFC, direction)

network <- graph_from_data_frame(string_filtered,
                                 vertices = protein_list,
                                 directed = FALSE)

V(network)$degree <- degree(network)
cat("  - Network nodes:", vcount(network), "\n")
cat("  - Network edges:", ecount(network), "\n")

# ============================================================================
# 4. Color Gradient Functions (Red/Blue)
# ============================================================================
cat("\n[4/6] Calculating red/blue gradients...\n")

# Up-regulated: Red gradient (light -> dark)
get_red_gradient <- function(logfc) {
  max_val <- max(logfc[logfc > 0])
  if(max_val == 0 | is.infinite(max_val)) return("#6b7280")
  normalized <- logfc / max_val
  colors <- colorRampPalette(c("#fca5a5", "#ef4444", "#dc2626", "#991b1b"))(100)
  result <- ifelse(logfc > 0,
                   colors[round(normalized * 99) + 1],
                   "#6b7280")
  return(result)
}

# Down-regulated: Blue gradient (light -> dark)
get_blue_gradient <- function(logfc) {
  min_val <- min(logfc[logfc < 0])
  if(is.infinite(min_val)) return("#6b7280")
  normalized <- abs(logfc / min_val)
  colors <- colorRampPalette(c("#93c5fd", "#3b82f6", "#2563eb", "#1e40af"))(100)
  result <- ifelse(logfc < 0,
                   colors[round(normalized * 99) + 1],
                   "#6b7280")
  return(result)
}

# Apply gradients
up_idx <- V(network)$direction == "Up" & V(network)$avg_logFC > 0
down_idx <- V(network)$direction == "Down" & V(network)$avg_logFC < 0

node_colors <- character(vcount(network))

if(sum(up_idx) > 0) {
  node_colors[up_idx] <- get_red_gradient(V(network)$avg_logFC[up_idx])
}

if(sum(down_idx) > 0) {
  node_colors[down_idx] <- get_blue_gradient(V(network)$avg_logFC[down_idx])
}

# Handle zero logFC
unknown_idx <- !(up_idx | down_idx)
node_colors[unknown_idx] <- "#6b7280"  # Gray for unknown/zero

V(network)$node_color <- node_colors

cat("  - Up-regulated:", sum(up_idx), "genes (red gradient)\n")
cat("  - Down-regulated:", sum(down_idx), "genes (blue gradient)\n")
cat("  - Unknown/Zero:", sum(unknown_idx), "genes (gray)\n")

# ============================================================================
# 5. Community-Based Layout
# ============================================================================
cat("\n[5/6] Creating community-based layout...\n")

# Detect communities for layout
community <- cluster_louvain(network)
V(network)$community <- as.numeric(membership(community))

n_comm <- length(unique(V(network)$community))

# Create node data frame
node_df <- data.frame(
  name = V(network)$name,
  degree = V(network)$degree,
  degree_rank = V(network)$degree_rank,
  community = V(network)$community,
  tier = V(network)$tier,
  stringsAsFactors = FALSE
) %>%
  arrange(community, desc(degree))

# Create custom circular layout
vertex_order <- match(node_df$name, V(network)$name)
n_nodes <- vcount(network)
angles <- seq(0, 2*pi, length.out = n_nodes + 1)[1:n_nodes]

layout_matrix <- matrix(0, nrow = n_nodes, ncol = 2)
layout_matrix[vertex_order, 1] <- cos(angles)
layout_matrix[vertex_order, 2] <- sin(angles)

# Rotate 180 degrees counterclockwise (multiply by -1)
V(network)$layout_x <- layout_matrix[, 1] * -1
V(network)$layout_y <- layout_matrix[, 2] * -1

cat("  - Circular layout with community ordering\n")
cat("  - Communities detected:", n_comm, "\n\n")

# ============================================================================
# 6. Premium Visualization
# ============================================================================
cat("[6/6] Generating final visualization...\n")

# Create label data frame (only Top 100)
label_df <- data.frame(
  name = V(network)$name,
  degree_rank = V(network)$degree_rank,
  tier = V(network)$tier,
  layout_x = V(network)$layout_x,
  layout_y = V(network)$layout_y,
  stringsAsFactors = FALSE
) %>%
  filter(degree_rank <= 100)

cat("  - Top 100 genes will be labeled\n")
cat("  - Remaining", vcount(network) - nrow(label_df), "genes unlabeled\n\n")

# Create the plot
p_final <- ggraph(network, layout = "manual",
                  x = V(network)$layout_x,
                  y = V(network)$layout_y) +
  
  # EDGES: Straight lines to avoid crossing the center
  geom_edge_link(aes(width = combined_score, alpha = combined_score, color = combined_score),
                 show.legend = FALSE) +
  scale_edge_width_continuous(range = c(0.3, 1.5), guide = "none") +
  scale_edge_alpha_continuous(range = c(0.2, 0.5), guide = "none") +
  scale_edge_color_gradient2(
    low = "#4e79a7",
    mid = "gray95",
    high = "#e15759",
    midpoint = 0.7
  ) +
  
  # NODES LAYER 1: Top 1-10 by degree (Triangle, LARGEST)
  geom_node_point(data = function(d) d[d$tier == "Top1-10", ],
                  aes(fill = node_color, color = node_color),
                  size = 6.5,
                  shape = 25,
                  stroke = 0.3,
                  alpha = 1.0,
                  show.legend = FALSE) +
  
  # NODES LAYER 2: Top 11-30 by degree (Square)
  geom_node_point(data = function(d) d[d$tier == "Top11-30", ],
                  aes(fill = node_color, color = node_color),
                  size = 4,
                  shape = 22,
                  stroke = 0.25,
                  alpha = 0.95,
                  show.legend = FALSE) +
  
  # NODES LAYER 3: Top 31-100 by degree (Circle)
  geom_node_point(data = function(d) d[d$tier == "Top31-100", ],
                  aes(fill = node_color, color = node_color),
                  size = 4,
                  shape = 19,
                  stroke = 0.15,
                  alpha = 1.0,
                  show.legend = FALSE) +
  
  # NODES LAYER 4: Top 100+ by degree (Circle, smallest)
  geom_node_point(data = function(d) d[d$tier == "Top100+", ],
                  aes(fill = node_color, color = node_color),
                  size = 3,
                  shape = 19,
                  stroke = 0.15,
                  alpha = 0.8,
                  show.legend = FALSE) +
  
  # Use identity scales to apply actual node colors
  scale_fill_identity() +
  scale_color_identity() +
  
  # LABELS: Only Top 100 genes
  geom_text(data = label_df,
            aes(x = layout_x * 1.06,
                y = layout_y * 1.06,
                label = name),
            size = 5.2,
            fontface = ifelse(label_df$degree_rank <= 30, "bold", "plain"),
            check_overlap = TRUE,
            show.legend = FALSE) +

  # Professional theme
  theme_void(base_size = 28) +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    plot.margin = margin(40, 40, 40, 60)
  ) +
  coord_fixed() +
  
  # Custom legend positioned at bottom-right
  # TOP 10 (Triangle)
  annotate("point", x = 1.25, y = -1.1, shape = 24, size = 5,
           fill = NA, color = "#6b7280", stroke = 0.5) +
  annotate("text", x = 1.12, y = -1.1, label = "TOP 10",
           hjust = 1, vjust = 0.5, size = 5, fontface = "bold", color = "#374151") +

  # TOP 11-30 (Square)
  annotate("point", x = 1.25, y = -1.03, shape = 22, size = 4.5,
           fill = NA, color = "#6b7280", stroke = 0.5) +
  annotate("text", x = 1.12, y = -1.03, label = "TOP 11-30",
           hjust = 1, vjust = 0.5, size = 5, color = "#374151") +

  # Regular (degree > 3) (Circle)
  annotate("point", x = 1.25, y = -0.96, shape = 1, size = 4,
           fill = NA, color = "#6b7280", stroke = 0.5) +
  annotate("text", x = 1.12, y = -0.96, label = "Regular (degree > 3)",
           hjust = 1, vjust = 0.5, size = 5, color = "#374151") +

  # Up-regulated (Red)
  annotate("point", x = 1.25, y = -0.89, shape = 19, size = 4,
           fill = NA, color = "#ef4444", stroke = 0.5) +
  annotate("text", x = 1.12, y = -0.89, label = "Up-regulated",
           hjust = 1, vjust = 0.5, size = 5, color = "#dc2626", fontface = "bold") +

  # Down-regulated (Blue)
  annotate("point", x = 1.25, y = -0.82, shape = 19, size = 4,
           fill = NA, color = "#2563eb", stroke = 0.5) +
  annotate("text", x = 1.12, y = -0.82, label = "Down-regulated",
           hjust = 1, vjust = 0.5, size = 5, color = "#2563eb", fontface = "bold")

# ============================================================================
# 7. Display and Save
# ============================================================================
cat("\n========================================\n")
cat("Displaying plot in R...\n")
cat("========================================\n\n")

# Display the plot in R
print(p_final)

# Save PDF
cat("\nSaving final plots...\n")
ggsave(
  filename = file.path(figures_dir, "STRING_Protein_Interaction.pdf"),
  plot = p_final,
  width = 18,
  height = 18,
  device = pdf
)

cat("  ✓ PDF saved to:", file.path(figures_dir, "STRING_Protein_Interaction.pdf"), "\n")

# Save hub gene information for heatmap generation
# Get logFC and direction from DEG data
hub_genes_output <- node_df %>%
  left_join(deg_data %>% select(gene, avg_logFC), by = c("name" = "gene")) %>%
  mutate(direction = ifelse(avg_logFC > 0, "Up", "Down")) %>%
  select(name, degree, degree_rank, tier, direction, avg_logFC) %>%
  arrange(degree_rank)

write.csv(hub_genes_output,
          file.path(tables_dir, "STRING_Protein_Interaction.csv"),
          row.names = FALSE)

cat("  ✓ Hub genes CSV saved to:", file.path(tables_dir, "STRING_Protein_Interaction.csv"), "\n")

cat("\n========================================\n")
cat("Analysis complete!\n")
cat("========================================\n")