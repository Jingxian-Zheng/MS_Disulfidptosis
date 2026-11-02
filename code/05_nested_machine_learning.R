if (.Platform$OS.type == "windows") suppressWarnings(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressPackageStartupMessages({
    library(glmnet)
    library(randomForest)
    library(e1071)
    library(gbm)
    library(kknn)
    library(MASS)
    library(rpart)
    library(jsonlite)
})
set.seed(20261026)
options(warn = 1)
root <- Sys.getenv("REVISION_ROOT", unset = normalizePath(getwd(), winslash = "/", mustWork = TRUE))
od <- file.path(root, "machine_learning")
sink(file.path(od, "nested_cv_log.txt"), split = TRUE)
readx <- function(ds) {
    as.matrix(read.csv(file.path(od, paste0(ds, "_features.csv")), row.names = 1, check.names = FALSE))
}
meta <- read.csv(file.path(od, "GSE17048_metadata.csv"))
x <- readx("GSE17048")
meta <- meta[match(rownames(x), meta$sample), ]
y <- as.integer(meta$group == "MS")
aucfun <- function(y, p) {
    if (any(!is.finite(p)) || length(unique(y)) < 2)
        return(NA_real_)
    n1 <- sum(y == 1)
    n0 <- sum(y == 0)
    (sum(rank(p, ties.method = "average")[y == 1]) - n1 * (n1 + 1)/2)/(n1 * n0)
}
folds <- function(y, k, seed) {
    set.seed(seed)
    f <- integer(length(y))
    for (g in unique(y)) {
        i <- which(y == g)
        f[i] <- sample(rep(seq_len(k), length.out = length(i)))
    }
    f
}
prepfit <- function(x) {
    med <- apply(x, 2, function(v) if (all(is.na(v)))
        0 else median(v, na.rm = TRUE))
    for (j in seq_len(ncol(x))) x[!is.finite(x[, j]), j] <- med[j]
    mu <- colMeans(x)
    ss <- apply(x, 2, sd)
    keep <- is.finite(ss) & ss > 1e-08
    list(med = med, mu = mu, sd = ss, keep = keep)
}
prepapply <- function(x, pr) {
    for (j in seq_len(ncol(x))) x[!is.finite(x[, j]), j] <- pr$med[j]
    z <- sweep(sweep(x[, pr$keep, drop = FALSE], 2, pr$mu[pr$keep], "-"), 2, pr$sd[pr$keep], "/")
    z
}
configs <- list()
add <- function(a, p) {
    configs[[length(configs) + 1]] <<- list(algorithm = a, parameter = p)
}
for (a in c("LASSO", "Ridge", "ElasticNet")) for (l in c(0.001, 0.01, 0.1, 1)) add(a, list(lambda = l, alpha = switch(a,
    LASSO = 1, Ridge = 0, ElasticNet = 0.5)))
for (m in c(3, 7)) add("RandomForest", list(mtry = m, ntree = 500))
for (C in c(0.1, 1, 10)) add("SVM", list(cost = C, gamma = 1/14))
for (nt in c(100, 300)) for (dp in c(1, 2)) add("GBM", list(n.trees = nt, depth = dp, shrinkage = 0.01, n.minobsinnode = 5))
for (k in c(3, 7, 11)) add("kNN", list(k = k))
add("NaiveBayes", list(laplace = 1))
add("LDA", list())
add("StepwiseLogistic", list())
for (cp in c(0.01, 0.05)) add("DecisionTree", list(cp = cp))
write_json(configs, file.path(od, "candidate_hyperparameter_grid.json"), pretty = TRUE, auto_unbox = TRUE)
fitone <- function(xx, yy, cfg, seed) {
    set.seed(seed)
    pr <- prepfit(xx)
    z <- prepapply(xx, pr)
    yf <- factor(ifelse(yy == 1, "MS", "Control"), levels = c("Control", "MS"))
    p <- cfg$parameter
    a <- cfg$algorithm
    dd <- data.frame(z, group = yf, check.names = FALSE)
    fit <- switch(a, LASSO = , Ridge = , ElasticNet = glmnet(z, yy, family = "binomial", alpha = p$alpha,
        lambda = c(p$lambda * 1.01, p$lambda), standardize = FALSE), RandomForest = randomForest(z, yf, ntree = p$ntree,
        mtry = min(p$mtry, ncol(z))), SVM = svm(z, yf, cost = p$cost, gamma = p$gamma, probability = TRUE,
        scale = FALSE), GBM = gbm(yy ~ ., data = data.frame(z, yy = yy), distribution = "bernoulli", n.trees = p$n.trees,
        interaction.depth = p$depth, shrinkage = p$shrinkage, n.minobsinnode = p$n.minobsinnode, bag.fraction = 0.5,
        verbose = FALSE), kNN = list(train = dd), NaiveBayes = naiveBayes(z, yf, laplace = p$laplace), LDA = lda(z,
        yf), StepwiseLogistic = suppressWarnings(step(glm(group ~ ., data = dd, family = binomial()), trace = FALSE)),
        DecisionTree = rpart(group ~ ., data = dd, method = "class", control = rpart.control(cp = p$cp, minsplit = 10,
            xval = 0)))
    list(model = fit, prep = pr, config = cfg)
}
predone <- function(obj, xx) {
    z <- prepapply(xx, obj$prep)
    p <- obj$config$parameter
    a <- obj$config$algorithm
    fit <- obj$model
    dd <- data.frame(z, check.names = FALSE)
    pred <- switch(a, LASSO = , Ridge = , ElasticNet = as.numeric(predict(fit, newx = z, s = p$lambda, type = "response")),
        RandomForest = predict(fit, z, type = "prob")[, "MS"], SVM = attr(predict(fit, z, probability = TRUE),
            "probabilities")[, "MS"], GBM = as.numeric(predict(fit, newdata = dd, n.trees = p$n.trees, type = "response")),
        kNN = kknn(group ~ ., train = fit$train, test = dd, k = p$k, distance = 2, kernel = "rectangular")$prob[,
            "MS"], NaiveBayes = predict(fit, z, type = "raw")[, "MS"], LDA = predict(fit, z)$posterior[, "MS"],
        StepwiseLogistic = predict(fit, newdata = dd, type = "response"), DecisionTree = predict(fit, newdata = dd,
            type = "prob")[, "MS"])
    pmin(pmax(as.numeric(pred), 1e-06), 1 - 1e-06)
}
selectcfg <- function(xx, yy, seed, label) {
    f <- folds(yy, 4, seed)
    scores <- matrix(NA_real_, length(configs), 4)
    for (j in seq_along(configs)) for (k in 1:4) {
        tr <- which(f != k)
        te <- which(f == k)
        res <- tryCatch({
            fit <- suppressWarnings(fitone(xx[tr, , drop = FALSE], yy[tr], configs[[j]], seed + 100 * j +
                k))
            aucfun(yy[te], predone(fit, xx[te, , drop = FALSE]))
        }, error = function(e) NA_real_)
        scores[j, k] <- res
    }
    means <- rowMeans(scores, na.rm = FALSE)
    idx <- which.max(replace(means, !is.finite(means), -Inf))
    tab <- data.frame(config = seq_along(configs), algorithm = sapply(configs, "[[", "algorithm"), parameters = sapply(configs,
        function(v) toJSON(v$parameter, auto_unbox = TRUE)), mean_inner_AUC = means, scores, check.names = FALSE)
    write.csv(tab, file.path(od, paste0("model_selection_", label, ".csv")), row.names = FALSE)
    write.csv(data.frame(sample = rownames(xx), group = yy, inner_fold = f), file.path(od, paste0("inner_folds_",
        label, ".csv")), row.names = FALSE)
    cat(label, "selected", configs[[idx]]$algorithm, toJSON(configs[[idx]]$parameter, auto_unbox = TRUE),
        "innerAUC", means[idx], "\n")
    flush.console()
    list(index = idx, table = tab)
}
outer <- folds(y, 5, 20261026)
oof <- rep(NA_real_, length(y))
outerselect <- list()
for (k in 1:5) {
    tr <- which(outer != k)
    te <- which(outer == k)
    choice <- selectcfg(x[tr, , drop = FALSE], y[tr], 20261026 + 1000 * k, paste0("outer", k))
    fit <- fitone(x[tr, , drop = FALSE], y[tr], configs[[choice$index]], 20261026 + k)
    oof[te] <- predone(fit, x[te, , drop = FALSE])
    outerselect[[k]] <- list(fold = k, config = configs[[choice$index]], outer_AUC = aucfun(y[te], oof[te]))
    saveRDS(fit, file.path(od, paste0("outer_fold_", k, "_model.rds")))
    cat("OUTER", k, "AUC", aucfun(y[te], oof[te]), "\n")
    flush.console()
}
write.csv(data.frame(sample = meta$sample, group = meta$group, y = y, outer_fold = outer, DRSS_oof = oof),
    file.path(od, "GSE17048_oof_DRSS.csv"), row.names = FALSE)
write_json(outerselect, file.path(od, "outer_model_selection.json"), pretty = TRUE, auto_unbox = TRUE)
choice <- selectcfg(x, y, 20262026, "full_training")
locked <- fitone(x, y, configs[[choice$index]], 20261026)
locked$feature_names <- colnames(x)
locked$fixed_threshold <- 0.5
locked$rank_universe <- readLines(file.path(od, "frozen_rank_gene_universe.txt"))
locked$locked_time <- as.character(Sys.time())
saveRDS(locked, file.path(od, "DRSS_locked_model.rds"))
# Only after the model is saved are external feature matrices and labels evaluated.
newmeta <- read.csv(file.path(od, "GSE203241_metadata.csv"))
newx <- readx("GSE203241")
newmeta <- newmeta[match(rownames(newx), newmeta$sample), ]
newmeta$DRSS <- predone(locked, newx)
write.csv(newmeta, file.path(od, "GSE203241_locked_predictions.csv"), row.names = FALSE)
bmeta <- read.csv(file.path(od, "GSE32915_metadata.csv"))
bx <- readx("GSE32915")
bmeta <- bmeta[match(rownames(bx), bmeta$sample), ]
bmeta$DRSS <- predone(locked, bx)
write.csv(bmeta, file.path(od, "GSE32915_locked_predictions.csv"), row.names = FALSE)
bdonor <- aggregate(DRSS ~ donor + group, bmeta, mean)
write.csv(bdonor, file.path(od, "GSE32915_donor_predictions.csv"), row.names = FALSE)
wilson <- function(a, n) {
    z <- qnorm(0.975)
    c((a/n + z * z/(2 * n) - z * sqrt(a/n * (1 - a/n)/n + z * z/(4 * n * n)))/(1 + z * z/n), (a/n + z * z/(2 *
        n) + z * sqrt(a/n * (1 - a/n)/n + z * z/(4 * n * n)))/(1 + z * z/n))
}
metric <- function(y, p, name) {
    set.seed(20261026)
    au <- aucfun(y, p)
    b <- replicate(2000, {
        i <- c(sample(which(y == 0), sum(y == 0), replace = TRUE), sample(which(y == 1), sum(y == 1), replace = TRUE))
        aucfun(y[i], p[i])
    })
    ci <- quantile(b, c(0.025, 0.975), na.rm = TRUE)
    sens <- mean(p[y == 1] >= 0.5)
    spec <- mean(p[y == 0] < 0.5)
    sci <- wilson(sum(p[y == 1] >= 0.5), sum(y == 1))
    tci <- wilson(sum(p[y == 0] < 0.5), sum(y == 0))
    lp <- qlogis(p)
    cal <- suppressWarnings(glm(y ~ lp, family = binomial()))
    cali <- suppressWarnings(glm(y ~ 1 + offset(lp), family = binomial()))
    data.frame(cohort = name, n = length(y), n_MS = sum(y == 1), n_control = sum(y == 0), AUC = au, AUC_low = ci[1],
        AUC_high = ci[2], Brier = mean((p - y)^2), sensitivity = sens, sens_low = sci[1], sens_high = sci[2],
        specificity = spec, spec_low = tci[1], spec_high = tci[2], calibration_intercept = coef(cali)[1],
        calibration_slope = coef(cal)[2], threshold = 0.5)
}
metrics <- rbind(metric(y, oof, "GSE17048_nested_OOF"), metric(as.integer(newmeta$group == "MS"), newmeta$DRSS,
    "GSE203241_external_all"), metric(as.integer(bdonor$group == "MS"), bdonor$DRSS, "GSE32915_previously_examined_donor_level"))
for (st in c("Adult", "Pediatric")) {
    i <- which(newmeta$age_stratum == st)
    metrics <- rbind(metrics, metric(as.integer(newmeta$group[i] == "MS"), newmeta$DRSS[i], paste0("GSE203241_external_",
        st)))
}
write.csv(metrics, file.path(root, "tables/S06_model_performance.csv"), row.names = FALSE)
print(metrics)
write_json(list(algorithm = locked$config$algorithm, parameters = locked$config$parameter, threshold = 0.5,
    seed = 20261026, outer_folds = 5, inner_folds = 4, external_missing_gene = "OSBPL5; imputed with training median",
    frozen_rank_universe = length(locked$rank_universe)), file.path(od, "final_model_parameters.json"), pretty = TRUE,
    auto_unbox = TRUE)
cat("ML_COMPLETE\n")
sink()
