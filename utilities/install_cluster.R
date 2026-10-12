# 1. Install BiocManager if you don't have it
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# 2. Install the missing translation packages
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db"))
