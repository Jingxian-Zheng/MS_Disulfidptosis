if (.Platform$OS.type == 'windows') suppressWarnings(Sys.setlocale('LC_CTYPE', 'English_United States.utf8'))
# Install primary direct dependencies; use R 4.3.x for Bioconductor 3.18.
script_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- dirname(dirname(normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/")))
versions <- read.csv(file.path(root, "environment/r-packages-observed.csv"), stringsAsFactors = FALSE)
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", repos = "https://cloud.r-project.org")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
if (getRversion() < "4.3" || getRversion() >= "4.4") stop("Use R 4.3.x for the recorded Bioconductor release")
bioc <- c("edgeR", "limma", "GSVA")
for (i in seq_len(nrow(versions))) {
    name <- versions$package[i]; version <- versions$source_version[i]
    if (requireNamespace(name, quietly = TRUE) && packageVersion(name) == package_version(version)) next
    if (name == "scTenifoldNet") next
    if (name %in% bioc) {
        BiocManager::install(name, version = "3.18", update = FALSE, ask = FALSE)
    } else if (name %in% c("IOBR", "CellChat")) {
        sha <- versions$remote_sha[i]; owner <- versions$remote_owner[i]; repo <- versions$remote_repo[i]
        if (!nzchar(sha) || !nzchar(owner) || !nzchar(repo)) {
            stop(paste("Install the recorded", name, version, "source manually; no verified source commit was recorded"))
        }
        remotes::install_github(paste0(owner, "/", repo, "@", sha), upgrade = "never")
    } else {
        remotes::install_version(name, version = version, repos = "https://cloud.r-project.org", upgrade = "never")
    }
}
for (pkg in c("RSpectra", "ps", "RhpcBLASctl")) {
    if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, repos = "https://cloud.r-project.org")
}
if (!requireNamespace("scTenifoldNet", quietly = TRUE) || as.character(packageVersion("scTenifoldNet")) != "1.4.3") {
    install.packages(file.path(root, "third_party/scTenifoldNet"), repos = NULL, type = "source")
}
cat("Inspect versions against environment/r-packages-observed.csv before analysis.\n")
