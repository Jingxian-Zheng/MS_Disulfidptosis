if (.Platform$OS.type == "windows") suppressWarnings(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressPackageStartupMessages({
    library(GSVA)
    library(IOBR)
    library(Matrix)
    library(jsonlite)
})
out <- Sys.getenv("REVISION_ROOT", unset = normalizePath(getwd(), winslash = "/", mustWork = TRUE))
src <- Sys.getenv("LEGACY_STUDY_ROOT", unset = file.path(Sys.getenv("REVISION_ROOT", unset = getwd()), "legacy_inputs"))
sink(file.path(out, "immune/case_only_log.txt"), split = TRUE)
source(file.path(Sys.getenv("MS_CODE_ROOT"), "code/stat_helpers.R"))
meta <- read.csv(file.path(out, "machine_learning/GSE17048_metadata.csv"))
oof <- read.csv(file.path(out, "machine_learning/original_expression_sensitivity/GSE17048_oof_DRSS.csv"))
meta <- merge(meta, oof[, c("sample", "outer_fold", "DRSS_oof")], by = "sample", sort = FALSE)
meta <- meta[match(oof$sample, meta$sample), ]
x <- as.matrix(read.delim(file.path(out, "machine_learning/GSE17048_gene_log2.tsv"), row.names = 1, check.names = FALSE))
common <- intersect(rownames(x), rownames(read.delim(file.path(src, "Data/GSE17048_mRNA_exp.tsv"), row.names = 1)))
x <- x[common, meta$sample, drop = FALSE]
meta$whole_cohort_group <- ifelse(meta$DRSS_oof > median(meta$DRSS_oof), "High", "Low")
ms <- meta$group == "MS"
meta$MS_case_group <- NA_character_
cutoff <- median(meta$DRSS_oof[ms])
meta$MS_case_group[ms] <- ifelse(meta$DRSS_oof[ms] > cutoff, "High", "Low")
write.csv(meta, file.path(out, "tables/S11_DRSS_scores_and_group_composition.csv"), row.names = FALSE)
comp <- rbind(data.frame(scope = "whole_cohort", as.data.frame(table(DRSS_group = meta$whole_cohort_group,
    disease = meta$group))), data.frame(scope = "MS_cases_only", as.data.frame(table(DRSS_group = meta$MS_case_group[ms],
    disease = meta$group[ms]))))
write.csv(comp, file.path(out, "tables/S11_group_composition_counts.csv"), row.names = FALSE)
old <- read.csv(file.path(out, "audit/original_DRSS.csv"))
old <- old[match(meta$sample, old$sample), ]
oldg <- ifelse(old$old_DRSS > median(old$old_DRSS), "High", "Low")
write.csv(as.data.frame(table(group = oldg, disease = meta$group)), file.path(out, "tables/S11_original_group_composition.csv"),
    row.names = FALSE)
mets <- metric(as.integer(meta$group == "MS"), meta$DRSS_oof, "Original_expression_nested_OOF")
write.csv(mets, file.path(out, "tables/S06_original_expression_nested_metrics.csv"), row.names = FALSE)
print(mets)
orient <- read.csv(file.path(out, "tables/S06_original_ROC_direction_audit.csv"))
immune <- read.delim(file.path(src, "Data/免疫相关基因集.txt"))
gl <- split(as.character(immune[[2]]), immune[[1]])
ss <- if (file.exists(file.path(out, "immune/recomputed_ssGSEA.rds"))) readRDS(file.path(out, "immune/recomputed_ssGSEA.rds")) else gsva(x,
    gl, method = "ssgsea", ssgsea.norm = TRUE, kcdf = "Gaussian", verbose = FALSE)
saveRDS(ss, file.path(out, "immune/recomputed_ssGSEA.rds"))
write.csv(t(ss), file.path(out, "immune/recomputed_ssGSEA.csv"))
cat("ssGSEA", dim(ss), "\n")
flush.console()
cib <- if (file.exists(file.path(out, "immune/recomputed_CIBERSORT.rds"))) readRDS(file.path(out, "immune/recomputed_CIBERSORT.rds")) else deconvo_tme(eset = 2^x,
    method = "cibersort", arrays = TRUE, tumor = FALSE, perm = 200, seed = 20261026)
saveRDS(cib, file.path(out, "immune/recomputed_CIBERSORT.rds"))
write.csv(cib, file.path(out, "immune/recomputed_CIBERSORT.csv"), row.names = FALSE)
print(names(cib))
cat("CIBERSORT_COMPLETE\n")
flush.console()
# All case-only tests use the same 99 cases and the same recorded score.
mid <- meta$sample[ms]
score <- meta$DRSS_oof[ms]
grp <- meta$MS_case_group[ms]
fold <- factor(meta$outer_fold[ms])
nH <- sum(grp == "High")
nL <- sum(grp == "Low")
compare <- function(mat, family) {
    rows <- lapply(rownames(mat), function(v) {
        vals <- mat[v, mid]
        h <- vals[grp == "High"]
        l <- vals[grp == "Low"]
        ww <- wilcox.test(h, l, exact = FALSE)
        delta <- 2 * as.numeric(ww$statistic)/(length(h) * length(l)) - 1
        fr <- rank(vals)
        model <- lm(fr ~ rank(score) + fold)
        fp <- summary(model)$coefficients[2, 4]
        cc <- cor.test(score, vals, method = "spearman", exact = FALSE)
        data.frame(feature = v, family = family, n_MS = length(vals), n_High = length(h), n_Low = length(l),
            mean_High = mean(h), mean_Low = mean(l), median_High = median(h), median_Low = median(l), median_difference = median(h) -
                median(l), Cliffs_delta = delta, group_P = ww$p.value, Spearman_rho = unname(cc$estimate),
            correlation_P = cc$p.value, fold_adjusted_P = fp)
    })
    d <- do.call(rbind, rows)
    d$test_status <- ifelse(is.finite(d$Spearman_rho), "estimated", "constant_feature")
    d$group_P[!is.finite(d$group_P)] <- 1
    d$correlation_P[!is.finite(d$correlation_P)] <- 1
    d$fold_adjusted_P[!is.finite(d$fold_adjusted_P)] <- 1
    d$Cliffs_delta[!is.finite(d$Cliffs_delta)] <- 0
    d$group_q_BH <- p.adjust(d$group_P, "BH")
    d$correlation_q_BH <- p.adjust(d$correlation_P, "BH")
    d$fold_adjusted_q_BH <- p.adjust(d$fold_adjusted_P, "BH")
    d
}
cp <- readLines(file.path(out, "audit/original_R4_target_checkpoints.txt"))
cov <- data.frame(gene = cp, available = cp %in% rownames(x))
write.csv(cov, file.path(out, "tables/S12_checkpoint_coverage.csv"), row.names = FALSE)
a <- compare(x[intersect(cp, rownames(x)), , drop = FALSE], "checkpoint")
b <- compare(ss, "ssGSEA")
cib <- as.data.frame(cib)
idcol <- which(vapply(cib, function(v) is.character(v) || is.factor(v), logical(1)))[1]
if (is.na(idcol)) idcol <- 1
rownames(cib) <- as.character(cib[[idcol]])
cols <- setdiff(names(cib)[vapply(cib, is.numeric, logical(1))], grep("P.value|Pvalue|Correlation|RMSE", names(cib),
    value = TRUE, ignore.case = TRUE))
cm <- t(as.matrix(cib[, cols, drop = FALSE]))
c <- compare(cm, "CIBERSORT")
write.csv(rbind(a, b, c), file.path(out, "tables/S12_case_only_immune_complete.csv"), row.names = FALSE)
marker <- c("AIF1", "GFAP", "CD3D", "MBP")
mr <- compare(x[intersect(marker, rownames(x)), , drop = FALSE], "blood_marker")
write.csv(mr, file.path(out, "tables/S13_blood_marker_correlations.csv"), row.names = FALSE)
# Spearman confidence intervals for the four explicitly discussed markers.
for (j in seq_len(nrow(mr))) {
    set.seed(20261026 + j)
    vals <- x[mr$feature[j], mid]
    boot <- replicate(2000, {
        i <- sample(seq_along(score), replace = TRUE)
        cor(score[i], vals[i], method = "spearman")
    })
    mr$rho_low[j] <- quantile(boot, 0.025, na.rm = TRUE)
    mr$rho_high[j] <- quantile(boot, 0.975, na.rm = TRUE)
}
write.csv(mr, file.path(out, "tables/S13_blood_marker_correlations.csv"), row.names = FALSE)
# Independent marker panels omit all classifier features.
core <- c("SLC7A11", "SLC3A2", "RPN1", "NCKAP1")
stromal <- c("COL1A1", "COL1A2", "DCN", "LUM", "LAMA2", "CEMIP", "ABCA9")
std <- t(scale(t(x)))
module <- rbind(genetic_core4 = colMeans(std[intersect(core, rownames(std)), , drop = FALSE]), independent_stromal = colMeans(std[intersect(stromal,
    rownames(std)), , drop = FALSE]))
write.csv(compare(module, "specificity"), file.path(out, "tables/S13_blood_score_specificity.csv"), row.names = FALSE)
# Reanalyse the complete archived GDSC2-derived prediction matrix, retaining its model-output scale.
dr <- read.csv(file.path(src, "Result5/calcPhenotype_Output/DrugPredictions.csv"), row.names = 1, check.names = FALSE)
dr <- as.matrix(dr[mid, , drop = FALSE])
write.csv(dr, file.path(out, "drug/MS_cases_untrimmed_predictions.csv"))
drugrows <- list()
trimrows <- list()
for (v in colnames(dr)) {
    vals <- dr[, v]
    h <- vals[grp == "High"]
    l <- vals[grp == "Low"]
    ww <- wilcox.test(h, l, exact = FALSE)
    cc <- cor.test(score, vals, method = "spearman", exact = FALSE)
    q1 <- quantile(vals, 0.25)
    q3 <- quantile(vals, 0.75)
    iq <- q3 - q1
    keep <- vals >= q1 - 1.5 * iq & vals <= q3 + 2.5 * iq
    wh <- wilcox.test(vals[keep & grp == "High"], vals[keep & grp == "Low"], exact = FALSE)
    drugrows[[v]] <- data.frame(drug = v, n_MS = length(vals), n_High = nH, n_Low = nL, median_High = median(h),
        median_Low = median(l), median_difference = median(h) - median(l), Cliffs_delta = 2 * as.numeric(ww$statistic)/(nH *
            nL) - 1, P = ww$p.value, Spearman_rho = unname(cc$estimate), correlation_P = cc$p.value)
    trimrows[[v]] <- data.frame(drug = v, lower_limit = q1 - 1.5 * iq, upper_limit = q3 + 2.5 * iq, n_excluded = sum(!keep),
        n_High_excluded = sum(!keep & grp == "High"), n_Low_excluded = sum(!keep & grp == "Low"), trimmed_P = wh$p.value,
        trimmed_median_difference = median(vals[keep & grp == "High"]) - median(vals[keep & grp == "Low"]))
}
d <- do.call(rbind, drugrows)
d$q_BH <- p.adjust(d$P, "BH")
d$correlation_q_BH <- p.adjust(d$correlation_P, "BH")
tr <- do.call(rbind, trimrows)
tr$trimmed_q_BH <- p.adjust(tr$trimmed_P, "BH")
write.csv(d, file.path(out, "tables/S14_all_198_drug_tests_untrimmed.csv"), row.names = FALSE)
write.csv(tr, file.path(out, "tables/S14_IQR_exclusion_sensitivity.csv"), row.names = FALSE)
target <- c("Rapamycin", "Cyclophosphamide", "Dasatinib", "Linsitinib", "Erlotinib", "Doramapimod", "BMS-345541",
    "Venetoclax")
clean <- function(s) tolower(gsub("[^A-Za-z0-9]", "", s))
sel <- which(vapply(d$drug, function(v) any(startsWith(clean(v), clean(target))), logical(1)))
write.csv(merge(d[sel, ], tr, by = "drug"), file.path(out, "tables/S14_eight_previously_highlighted_drugs.csv"),
    row.names = FALSE)
panel <- scan(file.path(src, "Result1/disulfidptosis_markers.txt"), what = "", quiet = TRUE)
corrows <- list()
for (g in panel) for (dd in d$drug[sel]) {
    cc <- cor.test(x[g, mid], dr[, dd], method = "spearman", exact = FALSE)
    corrows[[length(corrows) + 1]] <- data.frame(gene = g, drug = dd, rho = unname(cc$estimate), P = cc$p.value)
}
cc <- do.call(rbind, corrows)
cc$q_BH <- p.adjust(cc$P, "BH")
write.csv(cc, file.path(out, "tables/S14_gene_drug_correlations.csv"), row.names = FALSE)
summary <- list(n_MS = sum(ms), case_median_threshold = cutoff, n_High = nH, n_Low = nL, n_checkpoints = nrow(a),
    n_ssgsea = nrow(b), n_CIBERSORT = nrow(c), checkpoint_FDR05 = sum(a$group_q_BH < 0.05), ssgsea_FDR05 = sum(b$group_q_BH <
        0.05), CIBERSORT_FDR05 = sum(c$group_q_BH < 0.05), n_drugs = nrow(d), drug_FDR05 = sum(d$q_BH < 0.05),
    marker_results = mr, original_ROC = orient)
write_json(summary, file.path(out, "audit/case_only_results_summary.json"), pretty = TRUE, auto_unbox = TRUE)
print(summary)
cat("CASE_ONLY_COMPLETE\n")
sink()
