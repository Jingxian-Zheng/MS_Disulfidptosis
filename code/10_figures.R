set.seed(20261026)  # Deterministic figure jitter only.
if (.Platform$OS.type == "windows") suppressWarnings(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressPackageStartupMessages({
    library(ggplot2)
    library(patchwork)
    library(svglite)
    library(Matrix)
    library(pROC)
})
out <- Sys.getenv("REVISION_ROOT", unset = normalizePath(getwd(), winslash = "/", mustWork = TRUE))
setwd(out)
dir.create("figures/source", showWarnings = FALSE)
dir.create("figures/illustrator", showWarnings = FALSE)
sink("figures/plot_log.txt", split = TRUE)
read <- function(s) read.csv(file.path("tables", s), check.names = FALSE)
cols <- c(Control = "#4B6D8B", MS = "#C05A3E", High = "#C05A3E", Low = "#4B6D8B")
theme_set(theme_classic(base_size = 9, base_family = "Arial") + theme(plot.title = element_text(size = 10,
    face = "bold"), axis.text = element_text(size = 8, color = "black"), legend.text = element_text(size = 8),
    legend.title = element_text(size = 8), strip.text = element_text(size = 8), plot.margin = margin(6, 7,
        5, 5), legend.position = "bottom"))
savefig <- function(p, n, h = 170) {
    ggsave(paste0("figures/source/", n, ".svg"), p, width = 183, height = h, units = "mm", device = svglite,
        bg = "white")
    ggsave(paste0("figures/source/", n, ".pdf"), p, width = 183, height = h, units = "mm", device = cairo_pdf,
        bg = "white")
    ggsave(paste0("figures/source/", n, ".png"), p, width = 183, height = h, units = "mm", dpi = 180, bg = "white")
    cat(n, "exported\n")
    flush.console()
}
# Figure 1: donor structure and descriptive gene expression.
d <- read("S02_donor_cell_counts.csv")
md <- read("S02_sample_donor_metadata.csv")
don <- unique(md[, c("donor", "condition")])
d <- merge(d, don, by = "donor")
sc <- d[d$celltype == "SC", ]
sc$donor <- factor(sc$donor, levels = don$donor)
p1 <- ggplot(sc, aes(donor, Freq, fill = condition)) + geom_col(width = 0.7) + geom_text(aes(label = Freq),
    vjust = -0.4, size = 2.8) + scale_fill_manual(values = cols) + scale_y_continuous(expand = expansion(mult = c(0,
    0.18))) + labs(x = NULL, y = "Stromal nuclei", title = "A  Stromal nuclei by donor") + theme(axis.text.x = element_text(angle = 55,
    hjust = 1))
x <- readRDS("single_nucleus/selected_normalized_expression.rds")
m <- x$metadata
panel <- scan(file.path(Sys.getenv("LEGACY_STUDY_ROOT", unset = file.path(Sys.getenv("REVISION_ROOT", unset = getwd()),
    "legacy_inputs")), "Result1/disulfidptosis_markers.txt"), what = "", quiet = TRUE)
dot <- do.call(rbind, lapply(unique(m$celltype), function(t) do.call(rbind, lapply(panel, function(g) {
    v <- x$expression[g, m$celltype == t]
    data.frame(type = t, gene = g, mean = mean(v), pct = 100 * mean(v > 0))
}))))
dot$gene <- factor(dot$gene, levels = rev(panel))
dot$type <- factor(dot$type, levels = c("SC", "EC", "BC", "MG", "AS", "OPC", "OL", "NEU", "TC"))
write.csv(dot, "figures/source/Fig1_expression_data.csv", row.names = FALSE)
p2 <- ggplot(dot, aes(type, gene, size = pct, color = mean)) + geom_point() + scale_color_viridis_c(option = "C",
    name = "Mean log expression") + scale_size_area(max_size = 4, name = "% nuclei") + labs(x = NULL, y = NULL,
    title = "B  Discovery panel expression") + guides(color = guide_colorbar(title.position = "top", barwidth = grid::unit(20,
    "mm")), size = guide_legend(title.position = "top", nrow = 1)) + theme(legend.box = "vertical")
a <- read("S03_donor_celltype_module_scores.csv")
ss <- a[a$celltype == "SC", ]
ml <- rbind(data.frame(ss[, c("donor", "condition")], module = "Panel 14", score = ss$panel14), data.frame(ss[,
    c("donor", "condition")], module = "Genetic core 4", score = ss$genetic_core4))
ml$module <- factor(ml$module, levels = c("Panel 14", "Genetic core 4"))
p3 <- ggplot(ml, aes(condition, score, color = condition)) + geom_boxplot(outlier.shape = NA, width = 0.5) +
    geom_jitter(width = 0.1, size = 1.8) + facet_wrap(~module, scales = "free_y") + scale_color_manual(values = cols) +
    labs(x = NULL, y = "Donor mean log expression", title = "C  Donor stromal scores") + theme(legend.position = "none")
p4 <- ggplot(d, aes(celltype, donor, fill = log10(Freq + 1))) + geom_tile(color = "white", size = 0.2) + scale_fill_viridis_c(name = "log10(n + 1)") +
    labs(x = NULL, y = NULL, title = "D  Nuclei by lineage and donor")
savefig((p1 | p2)/(p3 | p4) + plot_layout(heights = c(1, 1.05)), "Fig1", 200)
# Figure 2: locked predictions, calibration, and confidence intervals.
o <- read.csv("machine_learning/original_expression_sensitivity/GSE17048_oof_DRSS.csv")
rk <- read.csv("machine_learning/GSE17048_oof_DRSS.csv")
ext <- read.csv("machine_learning/GSE203241_locked_predictions.csv")
rocdf <- function(y, s, n) {
    r <- roc(y, s, levels = c(0, 1), direction = "<", quiet = TRUE)
    data.frame(FPR = 1 - r$specificities, TPR = r$sensitivities, analysis = n)
}
rd <- rbind(rocdf(o$y, o$DRSS_oof, "Original expression OOF"), rocdf(rk$y, rk$DRSS_oof, "Rank harmonization OOF"),
    rocdf(as.integer(ext$group == "MS"), ext$DRSS, "New PBMC cohort"))
p1 <- ggplot(rd, aes(FPR, TPR, color = analysis)) + geom_abline(slope = 1, intercept = 0, color = "grey70",
    linetype = 2) + geom_step(size = 0.75) + scale_color_manual(values = c("#79706E", "#C05A3E", "#4B6D8B")) +
    coord_equal() + labs(x = "False positive rate", y = "True positive rate", color = NULL, title = "A  Discrimination") +
    guides(color = guide_legend(ncol = 1))
cb <- o
cb$bin <- pmin(5, ceiling(rank(cb$DRSS_oof, ties.method = "first")/nrow(cb) * 5))
cal <- aggregate(cbind(DRSS_oof, y) ~ bin, cb, mean)
cal$n <- as.numeric(table(cb$bin))
cal$se <- sqrt(cal$y * (1 - cal$y)/cal$n)
write.csv(cal, "figures/source/Fig2_calibration_bins.csv", row.names = FALSE)
p2 <- ggplot(cal, aes(DRSS_oof, y)) + geom_abline(slope = 1, intercept = 0, color = "grey70", linetype = 2) +
    geom_errorbar(aes(ymin = pmax(0, y - 1.96 * se), ymax = pmin(1, y + 1.96 * se)), width = 0.02) + geom_line(color = cols["MS"]) +
    geom_point(size = 2) + coord_equal(xlim = c(0, 1), ylim = c(0, 1)) + labs(x = "Mean OOF prediction", y = "Observed MS proportion",
    title = "B  OOF calibration")
perf <- read("S06_model_performance.csv")
pf <- rbind(data.frame(label = "Original OOF (144)", read("S06_original_expression_nested_metrics.csv")[,
    c("AUC", "AUC_low", "AUC_high")]), data.frame(label = c("Rank OOF (144)", "New PBMC (59)", "Brain donors (8)",
    "Adult PBMC (26)", "Pediatric PBMC (33)"), perf[, c("AUC", "AUC_low", "AUC_high")]))
pf$label <- factor(pf$label, levels = rev(pf$label))
p3 <- ggplot(pf, aes(AUC, label)) + geom_vline(xintercept = 0.5, linetype = 2, color = "grey60") + geom_errorbarh(aes(xmin = AUC_low,
    xmax = AUC_high), height = 0.15) + geom_point(size = 2, color = cols["MS"]) + scale_x_continuous(limits = c(0,
    1)) + labs(y = NULL, x = "AUC and bootstrap 95% CI", title = "C  AUC across settings")
p4 <- ggplot(o, aes(group, DRSS_oof, color = group)) + geom_boxplot(outlier.shape = NA, width = 0.5) + geom_jitter(width = 0.15,
    size = 1, alpha = 0.65) + scale_color_manual(values = cols) + labs(x = NULL, y = "OOF model probability",
    title = "D  OOF predictions") + theme(legend.position = "none")
savefig((p1 | p2)/(p3 | p4), "Fig2", 185)
# Figure 3: complete MS-only analyses.
im <- read("S12_case_only_immune_complete.csv")
im$family <- factor(im$family, levels = c("checkpoint", "ssGSEA", "CIBERSORT"))
p1 <- ggplot(im, aes(Cliffs_delta, -log10(group_q_BH), color = family)) + geom_hline(yintercept = -log10(0.05),
    linetype = 2, color = "grey50") + geom_point(size = 1.8, alpha = 0.8) + labs(x = "Cliff's delta (High minus Low)",
    y = "-log10(BH q)", color = NULL, title = "A  All case-only immune tests") + scale_color_manual(values = c("#5B8B80",
    "#C05A3E", "#4B6D8B"))
mk <- read("S13_blood_marker_correlations.csv")
mk$feature <- factor(mk$feature, levels = rev(c("AIF1", "GFAP", "CD3D", "MBP")))
p2 <- ggplot(mk, aes(Spearman_rho, feature)) + geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
    geom_errorbarh(aes(xmin = rho_low, xmax = rho_high), height = 0.15) + geom_point(size = 2, color = cols["MS"]) +
    labs(x = "Spearman rho and bootstrap 95% CI", y = NULL, title = "B  Blood marker correlations")
bx <- as.matrix(read.delim("machine_learning/GSE17048_gene_log2.tsv", row.names = 1, check.names = FALSE))
meta <- read("S11_DRSS_scores_and_group_composition.csv")
ms <- meta[meta$group == "MS", ]
sct <- rbind(data.frame(score = ms$DRSS_oof, expression = bx["CD3D", ms$sample], gene = "CD3D"), data.frame(score = ms$DRSS_oof,
    expression = bx["MBP", ms$sample], gene = "MBP"))
p3 <- ggplot(sct, aes(score, expression)) + geom_point(size = 1.25, color = cols["MS"], alpha = 0.65) + facet_wrap(~gene,
    scales = "free_y") + labs(x = "OOF probability in MS cases", y = "Blood log2 expression", title = "C  Case-only blood expression")
savefig((p1 | p2)/p3 + plot_layout(heights = c(1, 1)), "Fig3", 150)
# Figure 4: untrimmed 198-drug analysis and eight preselected examples.
dr <- read("S14_all_198_drug_tests_untrimmed.csv")
eight <- read("S14_eight_previously_highlighted_drugs.csv")
pred <- read.csv("drug/MS_cases_untrimmed_predictions.csv", row.names = 1, check.names = FALSE)
p1 <- ggplot(dr, aes(Cliffs_delta, -log10(q_BH))) + geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey60") +
    geom_point(size = 1.4, color = cols["MS"], alpha = 0.7) + labs(x = "Cliff's delta (High minus Low)", y = "-log10(BH q across 198 drugs)",
    title = "A  All 198 drug predictions")
long <- do.call(rbind, lapply(eight$drug, function(v) data.frame(sample = rownames(pred), value = pred[, v],
    drug = sub("_[0-9]+$", "", v))))
long <- merge(long, ms[, c("sample", "MS_case_group")], by = "sample")
long$MS_case_group <- factor(long$MS_case_group, levels = c("Low", "High"))
eight$label <- paste0("q = ", formatC(eight$q_BH, digits = 2, format = "f"))
eight$drugname <- sub("_[0-9]+$", "", eight$drug)
p2 <- ggplot(long, aes(MS_case_group, value, color = MS_case_group)) + geom_boxplot(outlier.shape = NA, width = 0.5) +
    geom_jitter(width = 0.13, size = 0.6, alpha = 0.5) + facet_wrap(~drug, scales = "free_y", ncol = 4) +
    scale_color_manual(values = cols) + labs(x = NULL, y = "Predicted IC50 model output", color = "Case group",
    title = "B  Eight previously highlighted compounds") + geom_text(data = data.frame(drug = eight$drugname,
    label = eight$label), aes(x = 1.5, y = Inf, label = label), inherit.aes = FALSE, vjust = 1.1, size = 2.8,
    color = "black")
savefig(p1/p2 + plot_layout(heights = c(0.8, 1.7)), "Fig4", 190)
# Figure 5: all-section estimates and virtual perturbation stability.
sp <- readRDS("spatial/all_spot_scores.rds")
rep <- sp[sp$sample == "MS197D", ]
maps <- rbind(data.frame(rep[, c("x", "y")], value = rep$panel13_no_DSE, variable = "Panel 13"), data.frame(rep[,
    c("x", "y")], value = rep$DSE, variable = "DSE expression"), data.frame(rep[, c("x", "y")], value = rep$CEBPB_target_activity,
    variable = "CEBPB target activity"))
maps$value <- ave(maps$value, maps$variable, FUN = function(v) as.numeric(scale(v)))
p1 <- ggplot(maps, aes(x, -y, color = value)) + geom_point(size = 0.4) + facet_wrap(~variable) + scale_color_viridis_c(name = "Within-map z") +
    coord_equal() + theme_void(base_family = "Arial", base_size = 9) + theme(strip.text = element_text(size = 8),
    legend.position = "bottom", plot.title = element_text(size = 10, face = "bold")) + labs(title = "A  Representative section MS197D")
st <- read("S08_all_section_spatial_associations.csv")
st <- st[st$condition == "MS" & st$module == "panel13_no_DSE" & st$target %in% c("CEBPB_target_activity",
    "RELA_target_activity"), ]
st$TF <- sub("_target_activity", "", st$target)
st$FDR <- ifelse(st$spatial_GLS_q_BH < 0.05, "q < 0.05", "q >= 0.05")
p2 <- ggplot(st, aes(spatial_GLS_beta, sample, color = FDR)) + geom_vline(xintercept = 0, linetype = 2, color = "grey65") +
    geom_errorbarh(aes(xmin = GLS_low, xmax = GLS_high), height = 0.15) + geom_point(size = 1.6) + facet_wrap(~TF) +
    scale_color_manual(values = c("#C05A3E", "#6B747D")) + labs(x = "Spatial GLS coefficient and 95% CI",
    y = NULL, color = NULL, title = "B  Spatial GLS in all MS sections")
r <- read("S09_donor_region_scores.csv")
r <- r[r$condition == "MS" & r$area %in% c("LC", "LR", "PPWM"), ]
r <- r[r$donor %in% names(which(table(r$donor) == 3)), ]
r$area <- factor(r$area, levels = c("LC", "LR", "PPWM"))
p3 <- ggplot(r, aes(area, panel14, group = donor, color = donor)) + geom_line(size = 0.5) + geom_point(size = 1.6) +
    labs(x = NULL, y = "Donor mean panel 14", color = "Donor", title = "C  Paired regions") + guides(color = guide_legend(ncol = 2)) +
    theme(legend.text = element_text(size = 8))
ko <- read("S10_DSE_virtual_KO_primary.csv")
ka <- read("S10_DSE_virtual_KO_alternate_seed.csv")
kl <- read("S10_DSE_virtual_KO_without_MS377.csv")
genes <- c("TIMP3", "FMO2", "CEMIP", "APOD", "KCNT2", "SAT1", "AC012459.1")
ks <- do.call(rbind, lapply(list(Primary = ko, Alternate_seed = ka, Without_MS377 = kl), function(z) z[z$gene %in%
    genes, c("gene", "p.adj")]))
ks$run <- rep(c("Primary", "Alternate seed", "Without MS377"), each = length(genes))
write.csv(ks, "figures/source/Fig5_KO_source.csv", row.names = FALSE)
ks$gene <- factor(ks$gene, levels = rev(genes))
ks$logq <- pmin(12, -log10(pmax(ks$p.adj, 1e-12)))
p4 <- ggplot(ks, aes(run, gene, color = logq, size = p.adj < 0.05)) + geom_point() + scale_size_manual(values = c(1.3,
    3), guide = "none") + scale_color_viridis_c(name = "-log10(q)\ncapped at 12") + labs(x = NULL, y = NULL,
    title = "D  Virtual perturbation") + guides(color = guide_colorbar(title.position = "top", barwidth = grid::unit(25,
    "mm"))) + theme(axis.text.x = element_text(angle = 25, hjust = 1))
savefig(p1/p2/(p3 | p4) + plot_layout(heights = c(0.75, 1.3, 1.2)), "Fig5", 218)
# Fig 6: original vector evidence map; all arrows denote analytical flow.
box <- data.frame(x = c(1, 4, 7, 1, 4, 7), y = c(3, 3, 3, 1, 1, 1), label = c("Single-nucleus atlas\n103,794 nuclei / 12 donors\n277 stromal nuclei",
    "All spatial sections\n12 MS sections / 7 donors\n6 control sections / 6 donors", "Blood analyses\n99 MS / 45 controls\nNew PBMC cohort: 59",
    "Observed expression\nStromal-associated panel 14\nSparse control sampling", "Spatial association\nHeterogeneous TF relationships\nNo uniform rim enrichment",
    "Model-derived estimates\nOOF AUC 0.527\nImmune/drug group tests\nBH q >= 0.05"))
p <- ggplot() + geom_rect(data = box, aes(xmin = x - 1.35, xmax = x + 1.35, ymin = y - 0.65, ymax = y + 0.65,
    fill = factor(y)), color = "#5B6570", size = 0.4) + geom_text(data = box, aes(x, y, label = label), family = "Arial",
    size = 3.1, lineheight = 1.3) + scale_fill_manual(values = c("#F1F3F5", "#E5EDF2"), guide = "none") +
    geom_segment(aes(x = c(1, 4, 7), xend = c(1, 4, 7), y = 2.3, yend = 1.75), arrow = arrow(length = grid::unit(2,
        "mm")), color = "#596574") + annotate("text", x = 4, y = 4.25, label = "Data and computational analysis",
    family = "Arial", size = 4, fontface = "bold") + annotate("text", x = 4, y = -0.35, label = "Virtual DSE perturbation: TIMP3 remained significant in three network runs\nCandidate network association for experimental testing",
    family = "Arial", size = 3.1, lineheight = 1.4) + annotate("segment", x = 1, xend = 4, y = 0.25, yend = 0.05,
    color = "#596574", arrow = arrow(length = grid::unit(2, "mm"))) + coord_cartesian(xlim = c(-0.6, 8.6),
    ylim = c(-1, 4.7), clip = "off") + theme_void() + theme(plot.margin = margin(10, 10, 10, 10))
savefig(p, "Fig6", 130)
# Supporting figures: complete section maps, markers, contribution sensitivity, core comparator.
sp$sample <- factor(sp$sample, levels = sort(unique(sp$sample)))
ps <- wrap_plots(lapply(levels(sp$sample), function(nm) ggplot(sp[sp$sample == nm, ], aes(x, -y, color = panel13_no_DSE)) +
    geom_point(size = 0.17) + scale_color_viridis_c(name = "Panel 13", limits = range(sp$panel13_no_DSE)) +
    coord_equal() + theme_void(base_family = "Arial", base_size = 9) + labs(title = nm) + theme(plot.title = element_text(size = 8),
    legend.position = "bottom")), ncol = 6, guides = "collect")
savefig(ps, "S1_Fig", 160)
ann <- read("S05_subtype_annotation_evidence.csv")
g <- unique(unlist(strsplit(ann$Markers, ";")))
a <- readRDS("single_nucleus/stromal_counts_and_normalized.rds")
dt <- do.call(rbind, lapply(unique(a$meta$subtype), function(st) {
    j <- which(a$meta$subtype == st)
    do.call(rbind, lapply(g, function(gg) {
        v <- if (gg %in% rownames(a$normalized))
            as.numeric(a$normalized[gg, j]) else rep(NA_real_, length(j))
        data.frame(subtype = st, gene = gg, n_nuclei = length(j), n_donors = length(unique(a$meta$donor[j])),
            mean_log_normalized = mean(v), percent_expressing = 100 * mean(v > 0), available = gg %in% rownames(a$normalized))
    }))
}))
write.csv(dt, "tables/S05_source_subtype_marker_validation.csv", row.names = FALSE)
dt <- dt[dt$available, ]
dt$gene <- factor(dt$gene, levels = g)
ps <- ggplot(dt, aes(gene, subtype, color = mean_log_normalized, size = percent_expressing)) + geom_point() +
    scale_color_viridis_c(name = "Mean log expression") + scale_size_area(max_size = 4, name = "% nuclei") +
    labs(x = NULL, y = NULL) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + guides(color = guide_colorbar(title.position = "top",
    barwidth = grid::unit(25, "mm")), size = guide_legend(title.position = "top", nrow = 1)) + theme(legend.box = "vertical")
savefig(ps, "S2_Fig", 145)
c <- read("S05_contribution_score_sensitivity.csv")
c$setting <- ifelse(c$definition == "original_log2_mean_ratio_times_mean_pct", "Original raw-count ratio",
    paste0("Stabilized; pseudocount ", c$pseudocount))
ps <- ggplot(c, aes(topK, rank, color = celltype)) + geom_line() + geom_point(size = 1.5) + facet_wrap(~setting,
    ncol = 2) + scale_y_reverse(breaks = 1:8) + scale_x_continuous(breaks = c(50, 100, 200)) + labs(x = "Top expressed genes per group",
    y = "Contribution rank", color = "Lineage")
savefig(ps, "S3_Fig", 145)
core <- st[FALSE, ]
allst <- read("S08_all_section_spatial_associations.csv")
core <- allst[allst$module == "genetic_core4", ]
core$target <- gsub("_target_activity", "", core$target)
ps <- ggplot(core, aes(spatial_GLS_beta, sample, color = condition)) + geom_vline(xintercept = 0, linetype = 2,
    color = "grey65") + geom_errorbarh(aes(xmin = GLS_low, xmax = GLS_high), height = 0.12, na.rm = TRUE) +
    geom_point(size = 1.3, na.rm = TRUE) + facet_wrap(~target, ncol = 2) + scale_color_manual(values = cols) +
    labs(x = "Spatial GLS coefficient and 95% CI", y = NULL, color = NULL)
savefig(ps, "S4_Fig", 190)
writeLines(capture.output(sessionInfo()), "audit/figure_sessionInfo.txt")
cat("ALL_FIGURES_COMPLETE\n")
sink()
