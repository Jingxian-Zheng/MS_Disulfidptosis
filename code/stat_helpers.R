aucfun <- function(y, p) {
    if (any(!is.finite(p)) || length(unique(y)) < 2)
        return(NA_real_)
    n1 <- sum(y == 1)
    n0 <- sum(y == 0)
    (sum(rank(p, ties.method = "average")[y == 1]) - n1 * (n1 + 1)/2)/(n1 * n0)
}
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
