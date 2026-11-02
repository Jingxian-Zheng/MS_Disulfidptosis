if (.Platform$OS.type == "windows") suppressWarnings(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
suppressPackageStartupMessages({
    library(hdf5r)
    library(Matrix)
    library(nlme)
    library(jsonlite)
    library(RANN)
})
out <- Sys.getenv("REVISION_ROOT", unset = normalizePath(getwd(), winslash = "/", mustWork = TRUE))
src <- file.path(out, "spatial/source_h5ad")
sink(file.path(out, "spatial/spatial_log.txt"), split = TRUE)
reg <- read.csv(file.path(Sys.getenv("MS_CODE_ROOT"), "resources/independent_TF_target_sets.csv"))
panel <- scan(file.path(Sys.getenv("LEGACY_STUDY_ROOT", unset = file.path(Sys.getenv("REVISION_ROOT", unset = getwd()),
    "legacy_inputs")), "Result1/disulfidptosis_markers.txt"), what = "", quiet = TRUE)
core <- c("SLC7A11", "SLC3A2", "RPN1", "NCKAP1")
reg <- reg[!reg$target %in% c(panel, core), ]
write.csv(reg, file.path(out, "tables/S07_independent_TF_target_sets.csv"), row.names = FALSE)
read_elem <- function(g) {
    if (inherits(g, "H5D"))
        return(g[])
    cc <- g[["codes"]][]
    ca <- g[["categories"]][]
    v <- ca[pmax(cc + 1, 1)]
    v[cc < 0] <- NA
    v
}
getcol <- function(f, n, default) {
    if (f$exists(n))
        read_elem(f[[n]]) else default
}
moran <- function(v, xy) {
    w <- nn2(xy, k = 7)$nn.idx[, -1, drop = FALSE]
    n <- length(v)
    W <- sparseMatrix(i = rep(seq_len(n), each = 6), j = as.vector(t(w)), x = 1, dims = c(n, n))
    W <- 1 * ((W + t(W)) > 0)
    diag(W) <- 0
    z <- v - mean(v)
    calc <- function(zz) as.numeric(n/sum(W) * crossprod(zz, W %*% zz)/sum(zz^2))
    obs <- calc(z)
    set.seed(20261026)
    pp <- replicate(199, calc(sample(z)))
    c(I = obs, P = (1 + sum(pp >= obs))/200)
}
all_spots <- list()
all_tests <- list()
all_regions <- list()
all_summary <- list()
moranrows <- list()
kk <- 0
files <- sort(list.files(src, pattern = "^visium_.*h5ad$", full.names = TRUE))
for (fi in seq_along(files)) {
    p <- files[fi]
    f <- H5File$new(p, "r")
    nm <- sub(".h5ad$", "", sub("visium_", "", basename(p)))
    genes <- f[["var/_index"]][]
    barcodes <- f[["obs/_index"]][]
    sh <- f[["X"]]$attr_open("shape")$read()
    ct <- new("dgCMatrix", i = as.integer(f[["X/indices"]][]), p = as.integer(f[["X/indptr"]][]), x = as.numeric(f[["X/data"]][]),
        Dim = as.integer(rev(sh)))
    rownames(ct) <- genes
    colnames(ct) <- barcodes
    lib <- Matrix::colSums(ct)
    norm <- ct %*% Diagonal(x = 10000/pmax(lib, 1))
    norm@x <- log1p(norm@x)
    xy <- f[["obsm/spatial"]]$read()
    if (nrow(xy) == 2)
        xy <- t(xy)
    if (nrow(xy) != ncol(ct))
        stop("Coordinate shape mismatch")
    n <- ncol(ct)
    areas <- getcol(f, "obs/areas", rep(NA_character_, n))
    niches <- getcol(f, "obs/niches", rep(NA_character_, n))
    patient <- getcol(f, "obs/patient_id", rep(nm, n))
    condition <- getcol(f, "obs/condition", rep(if (grepl("^MS", nm)) "MS" else "Control", n))
    dat <- data.frame(spot = barcodes, sample = nm, donor = sub(" .*", "", patient), condition = condition,
        areas = areas, niches = niches, x = xy[, 1], y = xy[, 2], library_size = lib, panel14 = Matrix::colSums(norm[intersect(panel,
            genes), , drop = FALSE])/length(panel), panel13_no_DSE = Matrix::colSums(norm[intersect(setdiff(panel,
            "DSE"), genes), , drop = FALSE])/(length(panel) - 1), genetic_core4 = Matrix::colMeans(norm[intersect(core,
            genes), , drop = FALSE]), DSE = as.numeric(norm["DSE", ]), CEBPB_expression = as.numeric(norm["CEBPB",
            ]), RELA_expression = as.numeric(norm["RELA", ]))
    dat$PROGENy_NFkB <- as.numeric(getcol(f, "obsm/progeny_scores/NFkB", rep(NA_real_, n)))
    dat$SC_proportion <- as.numeric(getcol(f, "obsm/ctype_props/SC", rep(NA_real_, n)))
    for (tf in c("CEBPB", "RELA")) {
        rr <- reg[reg$tf == tf & reg$target %in% genes, ]
        ex <- as.matrix(norm[rr$target, , drop = FALSE])
        mu <- rowMeans(ex)
        sdv <- apply(ex, 1, sd)
        good <- is.finite(sdv) & sdv > 1e-08
        zz <- sweep(sweep(ex[good, , drop = FALSE], 1, mu[good], "-"), 1, sdv[good], "/")
        dat[[paste0(tf, "_target_activity")]] <- as.numeric(crossprod(rr$mor[good], zz)/sum(abs(rr$mor[good])))
    }
    spacing <- median(nn2(xy, k = 2)$nn.dists[, 2])
    blocksize <- 5 * spacing
    dat$block <- paste(floor(dat$x/blocksize), floor(dat$y/blocksize), sep = "_")
    numvars <- c("x", "y", "library_size", "panel13_no_DSE", "genetic_core4", "DSE", "CEBPB_target_activity",
        "RELA_target_activity", "PROGENy_NFkB")
    block <- aggregate(dat[, numvars], list(block = dat$block), mean)
    block$n_spots <- as.numeric(table(factor(dat$block, levels = block$block)))
    block$x <- block$x/blocksize
    block$y <- block$y/blocksize
    for (mod in c("panel13_no_DSE", "genetic_core4")) {
        mi <- moran(dat[[mod]], xy)
        moranrows[[length(moranrows) + 1]] <- data.frame(sample = nm, module = mod, n_spots = n, Moran_I = mi[1],
            permutation_P = mi[2])
        for (target in c("DSE", "CEBPB_target_activity", "RELA_target_activity", "PROGENy_NFkB")) {
            kk <- kk + 1
            raw <- cor.test(dat[[mod]], dat[[target]], method = "pearson")
            b <- block
            b$Y <- as.numeric(scale(b[[mod]]))
            b$X <- as.numeric(scale(b[[target]]))
            b$loglib <- log10(pmax(b$library_size, 1))
            status <- "ok"
            g <- tryCatch(gls(Y ~ X + loglib, data = b, correlation = corExp(form = ~x + y, nugget = TRUE),
                method = "REML", control = glsControl(msMaxIter = 100, returnObject = TRUE)), error = function(e) {
                status <<- conditionMessage(e)
                NULL
            })
            beta <- se <- pval <- lo <- hi <- NA_real_
            if (!is.null(g)) {
                tt <- summary(g)$tTable
                beta <- tt["X", "Value"]
                se <- tt["X", "Std.Error"]
                pval <- tt["X", "p-value"]
                lo <- beta - qt(0.975, nrow(b) - 3) * se
                hi <- beta + qt(0.975, nrow(b) - 3) * se
            }
            all_tests[[kk]] <- data.frame(sample = nm, donor = unique(dat$donor)[1], condition = unique(condition)[1],
                module = mod, target = target, n_spots = n, n_blocks = nrow(block), blocksize_pixels = blocksize,
                raw_Pearson_r = unname(raw$estimate), raw_Pearson_P = raw$p.value, spatial_GLS_beta = beta,
                GLS_low = lo, GLS_high = hi, spatial_GLS_P = pval, fit_status = status)
        }
    }
    for (ar in unique(areas)) if (!is.na(ar)) {
        ii <- which(areas == ar)
        all_regions[[length(all_regions) + 1]] <- data.frame(sample = nm, donor = unique(dat$donor)[1], condition = unique(condition)[1],
            area = ar, n_spots = length(ii), panel14 = mean(dat$panel14[ii]), panel13_no_DSE = mean(dat$panel13_no_DSE[ii]),
            genetic_core4 = mean(dat$genetic_core4[ii]))
    }
    all_summary[[fi]] <- data.frame(sample = nm, donor = unique(dat$donor)[1], condition = unique(condition)[1],
        n_spots = n, n_genes = length(genes), panel_genes_detected = sum(panel %in% genes), core_genes_detected = sum(core %in%
            genes), areas = paste(sort(unique(areas)), collapse = ";"))
    all_spots[[fi]] <- dat
    f$close_all()
    cat("SPATIAL", nm, "n", n, "blocks", nrow(block), "\n")
    flush.console()
    rm(ct, norm)
    invisible(gc())
}
tests <- do.call(rbind, all_tests)
tests$raw_Pearson_q_BH <- p.adjust(tests$raw_Pearson_P, "BH", n = nrow(tests))
tests$spatial_GLS_q_BH <- p.adjust(tests$spatial_GLS_P, "BH", n = nrow(tests))
write.csv(tests, file.path(out, "tables/S08_all_section_spatial_associations.csv"), row.names = FALSE)
regions <- do.call(rbind, all_regions)
write.csv(regions, file.path(out, "tables/S09_all_section_region_scores.csv"), row.names = FALSE)
write.csv(do.call(rbind, all_summary), file.path(out, "tables/S07_spatial_section_manifest.csv"), row.names = FALSE)
mm <- do.call(rbind, moranrows)
mm$q_BH <- p.adjust(mm$permutation_P, "BH")
write.csv(mm, file.path(out, "tables/S08_spatial_autocorrelation.csv"), row.names = FALSE)
spots <- do.call(rbind, all_spots)
saveRDS(spots, file.path(out, "spatial/all_spot_scores.rds"))
writeLines("Fixed panel denominator in current full pipeline.", file.path(out, "spatial/fixed_denominator_applied.txt"))
write.csv(spots, gzfile(file.path(out, "spatial/all_spot_scores.csv.gz"), "wt"), row.names = FALSE)
dr <- aggregate(regions[, c("panel14", "panel13_no_DSE", "genetic_core4")], regions[, c("donor", "condition",
    "area")], mean)
write.csv(dr, file.path(out, "tables/S09_donor_region_scores.csv"), row.names = FALSE)
contrasts <- list()
for (mod in c("panel14", "panel13_no_DSE", "genetic_core4")) for (ref in c("LC", "PPWM")) {
    a <- dr[dr$condition == "MS" & dr$area == "LR", c("donor", mod)]
    b <- dr[dr$condition == "MS" & dr$area == ref, c("donor", mod)]
    v <- merge(a, b, by = "donor", suffixes = c("_LR", "_ref"))
    dif <- v[[2]] - v[[3]]
    ww <- wilcox.test(dif, mu = 0, exact = FALSE)
    contrasts[[length(contrasts) + 1]] <- data.frame(module = mod, contrast = paste0("LR_minus_", ref), n_donors = nrow(v),
        mean_difference = mean(dif), median_difference = median(dif), P = ww$p.value)
}
cc <- do.call(rbind, contrasts)
cc$q_BH <- p.adjust(cc$P, "BH")
write.csv(cc, file.path(out, "tables/S09_donor_paired_region_tests.csv"), row.names = FALSE)
print(cc)
cat("SPATIAL_COMPLETE\n")
sink()
