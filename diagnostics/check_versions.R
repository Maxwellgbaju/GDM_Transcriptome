cat("R version:", as.character(getRversion()), "\n")
for (pkg in c("DESeq2", "fgsea", "msigdbr", "limma", "dplyr")) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat(pkg, "version:", as.character(packageVersion(pkg)), "\n")
  } else {
    cat(pkg, "not installed\n")
  }
}
