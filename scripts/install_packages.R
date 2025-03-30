# Install and load BiocManager if not already installed
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install required packages
BiocManager::install(c(
  "edgeR",
  "limma",
  "DESeq2"
))

# Install CRAN packages
install.packages(c(
  "ggplot2",
  "pheatmap",
  "RColorBrewer",
  "tidyverse",
  "gridExtra",
  "stringr"
))

# Check if all packages can be loaded
required_packages <- c(
  "edgeR",
  "limma",
  "DESeq2",
  "ggplot2",
  "pheatmap",
  "RColorBrewer",
  "tidyverse",
  "gridExtra",
  "stringr"
)

# Try loading each package and report status
for (pkg in required_packages) {
  status <- tryCatch({
    library(pkg, character.only = TRUE)
    "Successfully loaded"
  }, error = function(e) {
    "Failed to load"
  })
  cat(sprintf("%s: %s\n", pkg, status))
}
