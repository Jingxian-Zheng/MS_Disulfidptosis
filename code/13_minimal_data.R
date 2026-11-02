if (.Platform$OS.type == "windows") suppressWarnings(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressPackageStartupMessages({
    library(Seurat)
    library(Matrix)
})
out <- Sys.getenv("REVISION_ROOT", unset = normalizePath(getwd(), winslash = "/", mustWork = TRUE))
a <- readRDS(Sys.getenv("ATLAS_RDS"))
m <- a@meta.data
ct <- a[["originalexp"]]@counts
nm <- ct %*% Diagonal(x = 10000/pmax(colSums(ct), 1))
nm@x <- log1p(nm@x)
g <- interaction(m$celltype, m$condition, drop = TRUE)
stats <- lapply(split(seq_len(ncol(ct)), g), function(j) list(n_nuclei = length(j), celltype = as.character(m$celltype[j[1]]),
    condition = as.character(m$condition[j[1]]), raw_mean = Matrix::rowMeans(ct[, j, drop = FALSE]), normalized_mean = Matrix::rowMeans(nm[,
        j, drop = FALSE]), n_detected = Matrix::rowSums(ct[, j, drop = FALSE] > 0)))
saveRDS(stats, file.path(out, "single_nucleus/lineage_condition_minimal_statistics.rds"))
cat("Minimal statistics for", length(stats), "lineage/condition groups saved\n")
writeLines(capture.output(sessionInfo()), file.path(out, "audit/minimal_data_sessionInfo.txt"))
