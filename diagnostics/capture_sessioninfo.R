cat("=== Primary differential-expression environment ===\n")
if (requireNamespace("DESeq2", quietly = TRUE)) {
  cat("DESeq2 version:", as.character(packageVersion("DESeq2")), "\n")
}
cat("R version:", as.character(getRversion()), "\n\n")
print(sessionInfo())
