if (.Platform$OS.type == "windows") suppressWarnings(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressPackageStartupMessages({
    library(Seurat)
    library(Matrix)
    library(edgeR)
    library(limma)
    library(jsonlite)
})
root <- Sys.getenv("LEGACY_STUDY_ROOT", unset = file.path(Sys.getenv("REVISION_ROOT", unset = getwd()), "legacy_inputs"))
out <- Sys.getenv("REVISION_ROOT", unset = normalizePath(getwd(), winslash = "/", mustWork = TRUE))
sink(file.path(out, "single_nucleus/reanalysis_log.txt"), split = TRUE)
sce <- readRDS(Sys.getenv("ATLAS_RDS"))
m <- sce@meta.data
m$donor <- sub(" .*", "", m$patient_id)
ct <- sce[["originalexp"]]@counts
stopifnot(all(abs(ct@x - round(ct@x)) < 1e-08))
cat("Counts", dim(ct), "\n")
cat("STROMAL CONDITION\n")
print(table(m$sample_id[m$celltype == "SC"], m$condition[m$celltype == "SC"]))
write.csv(as.data.frame(table(donor = m$donor, celltype = m$celltype)), file.path(out, "tables/S02_donor_cell_counts.csv"),
    row.names = FALSE)
write.csv(as.data.frame(table(sample = m$sample_id, celltype = m$celltype)), file.path(out, "tables/S02_sample_cell_counts.csv"),
    row.names = FALSE)
msc <- m[m$celltype == "SC", ]
msc$subtype <- as.character(msc$subtype)
msc$subtype[is.na(msc$subtype) | msc$subtype == "nan"] <- "Unassigned_control"
write.csv(as.data.frame(table(donor = msc$donor, subtype = msc$subtype)), file.path(out, "tables/S02_stromal_subtype_counts.csv"),
    row.names = FALSE)
write.csv(unique(m[, c("patient_id", "donor", "sample_id", "condition", "lesion_type", "age", "sex", "batch_sn")]),
    file.path(out, "tables/S02_sample_donor_metadata.csv"), row.names = FALSE)
# Aggregate original integer UMI counts by biological donor and cell type.
grp <- paste(m$donor, m$celltype, sep = "|")
lev <- unique(grp)
z <- sparseMatrix(i = seq_along(grp), j = match(grp, lev), x = 1, dims = c(length(grp), length(lev)))
pb <- ct %*% z
colnames(pb) <- lev
pmeta <- m[match(lev, grp), c("donor", "celltype", "condition", "age", "sex")]
rownames(pmeta) <- lev
pmeta$n_nuclei <- as.numeric(table(factor(grp, levels = lev)))
saveRDS(list(counts = pb, meta = pmeta), file.path(out, "single_nucleus/donor_pseudobulk_counts.rds"))
# Correct library normalization before any log-expression or marker calculation.
lib <- Matrix::colSums(ct)
norm <- ct %*% Diagonal(x = 10000/pmax(lib, 1))
norm@x <- log1p(norm@x)
sel <- scan(file.path(root, "Result1/disulfidptosis_markers.txt"), what = "", quiet = TRUE)
core <- c("SLC7A11", "SLC3A2", "RPN1", "NCKAP1")
canonical <- c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "LAMA2", "CEMIP", "ABCA9", "PDGFRA", "RGS5", "PDGFRB",
    "ACTA2", "TAGLN", "IL6", "CXCL12", "PTPRC", "VWF", "CLDN5", "MOG", "AQP4")
features <- unique(c(sel, core, canonical))
features <- intersect(features, rownames(norm))
small <- as.matrix(norm[features, , drop = FALSE])
saveRDS(list(expression = small, metadata = m, genes = features), file.path(out, "single_nucleus/selected_normalized_expression.rds"))
# Full donor x cell type module values; one donor carries one biological replicate.
mod <- data.frame(nucleus = colnames(ct), donor = m$donor, celltype = m$celltype, condition = m$condition,
    panel14 = colMeans(small[sel, , drop = FALSE]), genetic_core4 = colMeans(small[intersect(core, features),
        , drop = FALSE]))
agg <- aggregate(mod[, c("panel14", "genetic_core4")], mod[, c("donor", "celltype", "condition")], mean)
write.csv(agg, file.path(out, "tables/S03_donor_celltype_module_scores.csv"), row.names = FALSE)
# Leave-out-one-donor sensitivity and exact donor-level Wilcoxon for the sparse stromal group.
scd <- agg[agg$celltype == "SC", ]
print(scd)
mt <- lapply(c("panel14", "genetic_core4"), function(v) {
    a <- scd[scd$condition == "MS", v]
    b <- scd[scd$condition == "Control", v]
    data.frame(module = v, n_MS = length(a), n_control = length(b), mean_MS = mean(a), mean_control = mean(b),
        mean_difference = mean(a) - mean(b), wilcoxon_P = wilcox.test(a, b, exact = FALSE)$p.value)
})
mt <- do.call(rbind, mt)
mt$q_BH <- p.adjust(mt$wilcoxon_P, "BH")
write.csv(mt, file.path(out, "tables/S03_sparse_stromal_donor_module_tests.csv"), row.names = FALSE)
# Differential expression with prespecified minimum of ten nuclei per donor-celltype.
all_de <- list()
elig <- list()
for (ctype in unique(pmeta$celltype)) {
    ix <- which(pmeta$celltype == ctype & pmeta$n_nuclei >= 10)
    mm <- droplevels(pmeta[ix, ])
    mm$condition <- factor(mm$condition, levels = c("Control", "MS"))
    ng <- table(mm$condition)
    valid <- all(ng >= 3)
    elig[[ctype]] <- data.frame(celltype = ctype, minimum_nuclei = 10, n_control = unname(ng["Control"]),
        n_MS = unname(ng["MS"]), estimable = valid)
    if (!valid)
        next
    y <- DGEList(as.matrix(pb[, ix, drop = FALSE]))
    design <- model.matrix(~condition + scale(age) + factor(sex), mm)
    if (qr(design)$rank < ncol(design)) {
        design <- model.matrix(~condition, mm)
    }
    keep <- filterByExpr(y, design = design)
    y <- y[keep, , keep.lib.sizes = FALSE]
    y <- calcNormFactors(y)
    y <- estimateDisp(y, design, robust = TRUE)
    fit <- glmQLFit(y, design, robust = TRUE)
    test <- glmQLFTest(fit, coef = which(colnames(design) == "conditionMS"))
    tt <- topTags(test, n = Inf)$table
    tt$gene <- rownames(tt)
    tt$celltype <- ctype
    tt$model <- paste(colnames(design), collapse = ";")
    all_de[[ctype]] <- tt
    write.csv(tt, file.path(out, "single_nucleus", paste0("DE_", ctype, ".csv")), row.names = FALSE)
    cat("DE", ctype, nrow(tt), "FDR05", sum(tt$FDR < 0.05), "\n")
}
write.csv(do.call(rbind, elig), file.path(out, "tables/S04_pseudobulk_eligibility.csv"), row.names = FALSE)
write.csv(do.call(rbind, all_de), file.path(out, "tables/S04_donor_aware_DE_all_eligible_celltypes.csv"),
    row.names = FALSE)
# Normalized, descriptive subtype marker expression and total subtype coverage.
ix <- which(m$celltype == "SC")
scnorm <- norm[, ix, drop = FALSE]
saveRDS(list(counts = ct[, ix, drop = FALSE], normalized = scnorm, meta = msc), file.path(out, "single_nucleus/stromal_counts_and_normalized.rds"))
dots <- list()
k <- 0
for (st in unique(msc$subtype)) {
    j <- which(msc$subtype == st)
    for (g in features) {
        k <- k + 1
        v <- as.numeric(scnorm[g, j])
        dots[[k]] <- data.frame(subtype = st, gene = g, n_nuclei = length(j), n_donors = length(unique(msc$donor[j])),
            mean_log_normalized = mean(v), percent_expressing = mean(v > 0) * 100)
    }
}
write.csv(do.call(rbind, dots), file.path(out, "tables/S05_subtype_marker_expression.csv"), row.names = FALSE)
cat("SINGLE_NUCLEUS_COMPLETE\n")
sink()
