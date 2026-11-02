if (.Platform$OS.type == 'windows') suppressWarnings(Sys.setlocale('LC_CTYPE', 'English_United States.utf8'))
# Recompute key revision results from distributed processed data; no refitting of ML.
out <- Sys.getenv("REVISION_ROOT")
if (!nzchar(out)) stop("Use run.py to select a prepared workspace")
read_table <- function(name) read.csv(file.path(out, "tables", name), check.names = FALSE)
checks <- list()
check <- function(name, observed, expected, tolerance = 1e-10) {
    passed <- isTRUE(all.equal(observed, expected, tolerance = tolerance, check.attributes = FALSE))
    checks[[length(checks) + 1L]] <<- data.frame(check = name, passed = passed, observed = paste(observed,
        collapse = ";"), expected = paste(expected, collapse = ";"))
}
auc <- function(y, p) {
    n1 <- sum(y == 1)
    n0 <- sum(y == 0)
    (sum(rank(p, ties.method = "average")[y == 1]) - n1 * (n1 + 1)/2)/(n1 * n0)
}
d <- read_table("S03_donor_celltype_module_scores.csv")
d <- d[d$celltype == "SC", ]
tests <- lapply(c("panel14", "genetic_core4"), function(module) {
    wilcox.test(d[d$condition == "MS", module], d[d$condition == "Control", module], exact = FALSE)$p.value
})
check("Stromal donor module BH", p.adjust(unlist(tests), "BH"), read_table("S03_sparse_stromal_donor_module_tests.csv")$q_BH)
oof <- read.csv(file.path(out, "machine_learning/original_expression_sensitivity/GSE17048_oof_DRSS.csv"))
check("Primary OOF AUC", auc(oof$y, oof$DRSS_oof), read_table("S06_original_expression_nested_metrics.csv")$AUC)
check("Primary OOF Brier", mean((oof$DRSS_oof - oof$y)^2), read_table("S06_original_expression_nested_metrics.csv")$Brier)
external <- read.csv(file.path(out, "machine_learning/GSE203241_locked_predictions.csv"))
probcol <- intersect(c("DRSS", "prediction", "probability", "prob", "DRSS_locked"), names(external))
if (!length(probcol)) probcol <- grep("DRSS|prob|pred", names(external), value = TRUE)
if (length(probcol)) check("New PBMC constant predictions", length(unique(external[[probcol[1]]])), 1L)
ms <- oof[oof$y == 1, ]
high <- ms$DRSS_oof > median(ms$DRSS_oof)
check("MS-only group sizes", c(sum(high), sum(!high)), c(49L, 50L))
x <- read.csv(file.path(out, "drug/MS_cases_untrimmed_predictions.csv"), row.names = 1, check.names = FALSE)
x <- x[ms$sample, , drop = FALSE]
drug_p <- vapply(x, function(v) wilcox.test(v[high], v[!high], exact = FALSE)$p.value, 0)
drug_q <- p.adjust(drug_p, "BH")
expected <- read_table("S14_all_198_drug_tests_untrimmed.csv")
check("All 198 drug BH values", unname(drug_q[expected$drug]), expected$q_BH)
check("Corrected drug group discoveries", sum(drug_q < 0.05), 0L)
sp <- read_table("S08_all_section_spatial_associations.csv")
check("Planned spatial family size", nrow(sp), 144L)
check("Spatial GLS BH with failed fits retained", p.adjust(sp$spatial_GLS_P, "BH", n = nrow(sp)), sp$spatial_GLS_q_BH)
check("Failed spatial fits retained", sum(is.na(sp$spatial_GLS_P)), 7L)
for (file in list.files(file.path(out, "tables"), pattern = "^S12_.*[.]csv$", full.names = TRUE)) {
    a <- read.csv(file)
    if (all(c("group_P", "group_q_BH") %in% names(a))) {
        families <- if ("family" %in% names(a)) split(a, a$family) else list(all = a)
        for (family in names(families)) {
            block <- families[[family]]
            check(paste(basename(file), family, "group BH"), p.adjust(block$group_P, "BH"), block$group_q_BH)
            check(paste(basename(file), family, "group discoveries"), sum(block$group_q_BH < .05), 0L)
        }
    }
}
genes <- lapply(c("primary", "alternate_seed", "without_MS377"), function(label) {
    a <- read_table(paste0("S10_DSE_virtual_KO_", label, ".csv"))
    a$gene[a$gene != "DSE" & a$p.adj < 0.05]
})
check("Nonself KO candidate shared across three runs", sort(Reduce(intersect, genes)), "TIMP3")
report <- do.call(rbind, checks)
write.csv(report, file.path(out, "reports/processed_validation.csv"), row.names = FALSE)
print(report[, c("check", "passed")], row.names = FALSE)
if (!all(report$passed)) stop("Processed-data validation failed; inspect reports/processed_validation.csv")
cat("All processed-data checks passed\n")
