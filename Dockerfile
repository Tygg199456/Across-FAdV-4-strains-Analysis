# Dockerfile for FAdV-4 Limma Analysis Pipeline
# Build: docker build -t fadv4-analysis:latest .
# Run: docker run -v $(pwd)/data:/app/data -v $(pwd)/results:/app/results fadv4-analysis:latest

FROM rocker/r-ver:4.2.0

# Set working directory
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libcairo2-dev \
    libpango1.0-dev \
    libjpeg-dev \
    libtiff-dev \
    && rm -rf /var/lib/apt/lists/*

# Install BiocManager
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org', quiet=TRUE)"

# Install Bioconductor packages (CORRECTED: org.Gg.eg.db for Gallus gallus)
RUN R -e "BiocManager::install(c('limma', 'edgeR', 'clusterProfiler', \
                                   'org.Gg.eg.db', 'annotate', 'DOSE', \
                                   'ReactomePA', 'enrichplot', 'pathview', \
                                   'ggplot2', 'pheatmap', 'ComplexHeatmap'), \
                                   ask=FALSE, update=TRUE, quiet=TRUE)"

# Install CRAN packages
RUN R -e "install.packages(c('tidyverse', 'ggpubr', 'reshape2', 'dplyr', \
                              'tidyr', 'readr', 'data.table', 'igraph', \
                              'ggraph', 'tidygraph', 'corrplot', 'VennDiagram', \
                              'ggrepel', 'scales', 'RColorBrewer', 'viridis', \
                              'gridExtra', 'patchwork'), \
                              repos='https://cloud.r-project.org', quiet=TRUE)"

# Copy project files
COPY scripts/ /app/scripts/
COPY data/ /app/data/
COPY results/ /app/results/

# Set permissions
RUN chmod -R 755 /app/scripts/*.R

# Create output directories if they don't exist
RUN mkdir -p /app/results/figures /app/results/tables /app/results/intermediate

# Default command
CMD ["R", "--version"]
