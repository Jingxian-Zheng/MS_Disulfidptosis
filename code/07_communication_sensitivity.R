if (.Platform$OS.type == "windows") suppressWarnings(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressPackageStartupMessages({
    library(Seurat)
    library(Matrix)
    library(CellChat)
    library(jsonlite)
})
out <- Sys.getenv("REVISION_ROOT", unset = normalizePath(getwd(), winslash = "/", mustWork = TRUE))
sink(file.path(out, "single_nucleus/communication_and_sensitivity_log.txt"), split = TRUE)
sce <- readRDS(Sys.getenv("ATLAS_RDS"))
m <- sce@meta.data
m$donor <- sub(" .*", "", m$patient_id)
ct <- sce[["originalexp"]]@counts
norm <- ct %*% Diagonal(x = 10000/pmax(Matrix::colSums(ct), 1))
norm@x <- log1p(norm@x)
dimnames(norm) <- dimnames(ct)
# Correct missing subtype annotation without manufacturing control subclusters.
i <- which(m$celltype == "SC")
mm <- m[i, ]
mm$subtype <- as.character(mm$subtype)
mm$subtype[is.na(mm$subtype) | mm$subtype == "nan"] <- "Unassigned_control"
write.csv(as.data.frame(table(donor = mm$donor, subtype = mm$subtype)), file.path(out, "tables/S02_stromal_subtype_counts.csv"),
    row.names = FALSE)
a <- readRDS(file.path(out, "single_nucleus/stromal_counts_and_normalized.rds"))
a$meta <- mm
saveRDS(a, file.path(out, "single_nucleus/stromal_counts_and_normalized.rds"))
rm(a)
features <- unique(c(scan(file.path(Sys.getenv("LEGACY_STUDY_ROOT", unset = file.path(Sys.getenv("REVISION_ROOT",
    unset = getwd()), "legacy_inputs")), "Result1/disulfidptosis_markers.txt"), what = "", quiet = TRUE),
    c("COL1A1", "COL1A2", "DCN", "LUM", "LAMA2", "CEMIP", "ABCA9", "ACTA2", "TAGLN", "IL6", "CXCL12", "PTPRC",
        "VWF", "PDGFRA", "PDGFRB", "RGS5")))
features <- intersect(features, rownames(norm))
dots <- list()
for (st in unique(mm$subtype)) {
    j <- which(mm$subtype == st)
    for (g in features) {
        v <- as.numeric(norm[g, i[j]])
        dots[[length(dots) + 1]] <- data.frame(subtype = st, gene = g, n_nuclei = length(j), n_donors = length(unique(mm$donor[j])),
            mean_log_normalized = mean(v), percent_expressing = mean(v > 0) * 100)
    }
}
write.csv(do.call(rbind, dots), file.path(out, "tables/S05_subtype_marker_expression.csv"), row.names = FALSE)
# Reproduce the original top-expression calculation and test stabilized alternative definitions.
cs <- list()
for (tp in unique(m$celltype)) {
    ic <- which(m$celltype == tp & m$condition == "Control")
    im <- which(m$celltype == tp & m$condition == "MS")
    if (!length(ic) || !length(im))
        next
    rawc <- Matrix::rowMeans(ct[, ic, drop = FALSE])
    rawm <- Matrix::rowMeans(ct[, im, drop = FALSE])
    lc <- Matrix::rowMeans(norm[, ic, drop = FALSE])
    lm <- Matrix::rowMeans(norm[, im, drop = FALSE])
    for (k in c(50, 100, 200)) {
        top <- unique(c(names(sort(rawc, decreasing = TRUE))[seq_len(k)], names(sort(rawm, decreasing = TRUE))[seq_len(k)]))
        freq <- Matrix::rowMeans(ct[top, c(ic, im), drop = FALSE] > 0)
        original <- log2(mean((rawm[top] + 1e-10)/(rawc[top] + 1e-10))) * mean(freq)
        cs[[length(cs) + 1]] <- data.frame(celltype = tp, topK = k, definition = "original_log2_mean_ratio_times_mean_pct",
            pseudocount = 1e-10, n_genes = length(top), score = original, n_control_cells = length(ic), n_MS_cells = length(im))
        for (pc in c(0.01, 0.1, 1)) {
            sc <- mean(log2((lm[top] + pc)/(lc[top] + pc)) * freq)
            cs[[length(cs) + 1]] <- data.frame(celltype = tp, topK = k, definition = "mean_stabilized_gene_log_ratio_times_pct",
                pseudocount = pc, n_genes = length(top), score = sc, n_control_cells = length(ic), n_MS_cells = length(im))
        }
    }
}
cs <- do.call(rbind, cs)
cs$rank <- ave(-cs$score, interaction(cs$topK, cs$definition, cs$pseudocount), FUN = rank)
write.csv(cs, file.path(out, "tables/S05_contribution_score_sensitivity.csv"), row.names = FALSE)
# Identical non-B nuclei are used in both networks.
set.seed(20261026)
sel <- unlist(lapply(split(seq_len(nrow(m)), m$celltype), function(v) if (length(v) > 2500) sample(v, 2500) else v))
write.csv(data.frame(nucleus = rownames(m)[sel], celltype = m$celltype[sel], donor = m$donor[sel]), file.path(out,
    "single_nucleus/communication_nuclei.csv"), row.names = FALSE)
DB <- subsetDB(CellChatDB.human, search = "Secreted Signaling")
run_net <- function(includeB) {
    jj <- sel
    if (!includeB)
        jj <- jj[m$celltype[jj] != "BC"]
    meta <- data.frame(celltype = droplevels(m$celltype[jj]), samples = m$donor[jj], row.names = rownames(m)[jj])
    obj <- createCellChat(object = norm[, jj, drop = FALSE], meta = meta, group.by = "celltype")
    obj@DB <- DB
    obj <- subsetData(obj)
    obj <- identifyOverExpressedGenes(obj)
    obj <- identifyOverExpressedInteractions(obj)
    obj <- computeCommunProb(obj, raw.use = TRUE, population.size = FALSE, nboot = 100, seed.use = 20261026)
    obj <- filterCommunication(obj, min.cells = 10)
    obj <- computeCommunProbPathway(obj)
    obj <- aggregateNet(obj)
    label <- if (includeB)
        "including_B" else "excluding_B"
    saveRDS(list(net = obj@net, netP = obj@netP, LR = obj@LR, groups = levels(obj@idents)), file.path(out,
        "single_nucleus", paste0("CellChat_", label, ".rds")))
    pr <- obj@net$prob
    pv <- obj@net$pval
    dt <- expand.grid(receiver = dimnames(pr)[[2]], interaction = dimnames(pr)[[3]], stringsAsFactors = FALSE)
    dt$probability <- as.numeric(pr["SC", , ])
    dt$P_cell_bootstrap <- as.numeric(pv["SC", , ])
    dt$q_BH <- p.adjust(dt$P_cell_bootstrap, "BH")
    dt$network <- label
    dt
}
a <- run_net(TRUE)
b <- run_net(FALSE)
write.csv(rbind(a, b), file.path(out, "tables/S05_B_cell_communication_sensitivity.csv"), row.names = FALSE)
v <- merge(a, b, by = c("receiver", "interaction"), suffixes = c("_withB", "_withoutB"))
v <- v[v$probability_withB > 0 | v$probability_withoutB > 0, ]
sumry <- data.frame(n_common_nonzero_pairs = nrow(v), Spearman_r = cor(v$probability_withB, v$probability_withoutB,
    method = "spearman"), mean_abs_probability_difference = mean(abs(v$probability_withB - v$probability_withoutB)),
    n_B_receiving_pairs = sum(a$receiver == "BC" & a$probability > 0))
write.csv(sumry, file.path(out, "tables/S05_B_cell_network_summary.csv"), row.names = FALSE)
print(sumry)
cat("COMMUNICATION_COMPLETE\n")
sink()
