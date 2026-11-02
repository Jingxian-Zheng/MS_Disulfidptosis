if (.Platform$OS.type == "windows") suppressWarnings(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressPackageStartupMessages({
    library(Seurat)
    library(Matrix)
    library(scTenifoldNet)
    library(cli)
    library(reshape2)
    library(igraph)
    library(jsonlite)
})
out <- Sys.getenv("REVISION_ROOT", unset = normalizePath(getwd(), winslash = "/", mustWork = TRUE))
sink(file.path(out, "single_nucleus/virtual_KO_log.txt"), split = TRUE)
# Source the archived official functions directly; enrichment web-service wrappers are not called.
for (f in list.files(file.path(Sys.getenv("MS_CODE_ROOT"), "third_party/scTenifoldKnk/R"), full.names = TRUE,
    pattern = "R$")) sys.source(f, envir = .GlobalEnv)
a <- readRDS(file.path(out, "single_nucleus/stromal_counts_and_normalized.rds"))
i <- which(a$meta$condition == "MS")
ct <- a$counts[, i, drop = FALSE]
meta <- a$meta[i, ]
obj <- CreateSeuratObject(ct)
obj <- NormalizeData(obj, verbose = FALSE)
obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 1000, verbose = FALSE)
panel <- scan(file.path(Sys.getenv("LEGACY_STUDY_ROOT", unset = file.path(Sys.getenv("REVISION_ROOT", unset = getwd()),
    "legacy_inputs")), "Result1/disulfidptosis_markers.txt"), what = "", quiet = TRUE)
core <- c("SLC7A11", "SLC3A2", "RPN1", "NCKAP1")
genes <- unique(c(VariableFeatures(obj), panel, core, "SAT1", "AC012459.1"))
eligible <- rownames(ct)[Matrix::rowMeans(ct > 0) >= 0.05]
genes <- intersect(genes, eligible)
stopifnot("DSE" %in% genes)
writeLines(genes, file.path(out, "single_nucleus/virtual_KO_gene_universe.txt"))
configs <- list(list(name = "primary", seed = 20261026, keep = seq_len(ncol(ct))), list(name = "alternate_seed",
    seed = 20261027, keep = seq_len(ncol(ct))), list(name = "without_MS377", seed = 20261026, keep = which(meta$donor !=
    "MS377")))
params <- list(n_genes = length(genes), HVG = 1000, minimum_detection_fraction = 0.05, nNet = 10, nCells = 150,
    nComp = 3, td_K = 3, ma_nDim = 2, ma_method = "manifold", dr_direction = FALSE, nCores = 1, base_seed = 20261026)
write_json(params, file.path(out, "single_nucleus/virtual_KO_parameters.json"), pretty = TRUE, auto_unbox = TRUE)
for (cfg in configs) {
    cat("KO START", cfg$name, "cells", length(cfg$keep), "\n")
    flush.console()
    res <- scTenifoldKnk(countMatrix = ct[genes, cfg$keep, drop = FALSE], gKO = "DSE", qc = FALSE, nc_nNet = 10,
        nc_nCells = 150, nc_nComp = 3, td_K = 3, td_maxIter = 1000, ma_nDim = 2, ma_method = "manifold", dr_direction = FALSE,
        nCores = 1, seed = cfg$seed)
    saveRDS(res, file.path(out, "single_nucleus", paste0("DSE_virtual_KO_", cfg$name, ".rds")))
    tt <- as.data.frame(res$diffRegulation)
    tt$gene_id <- if ("gene" %in% names(tt))
        tt$gene else rownames(tt)
    write.csv(tt, file.path(out, "tables", paste0("S10_DSE_virtual_KO_", cfg$name, ".csv")), row.names = FALSE)
    print(head(tt))
    cat("KO COMPLETE", cfg$name, "\n")
    flush.console()
}
cat("VIRTUAL_KO_COMPLETE\n")
sink()
