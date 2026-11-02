if (.Platform$OS.type == 'windows') suppressWarnings(Sys.setlocale('LC_CTYPE', 'English_United States.utf8'))
# Parse all project and vendored R source files without loading analysis packages.
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[1], winslash = "/") else getwd()
files <- list.files(root, pattern = "[.][Rr]$", recursive = TRUE, full.names = TRUE)
files <- files[!grepl("/(work|[.]git)/", files)]
errors <- character()
for (file in files) {
    tryCatch(parse(file, encoding = "UTF-8"), error = function(e) {
        errors <<- c(errors, paste(file, conditionMessage(e)))
    })
}
if (length(errors)) stop(paste(errors, collapse = "\n"))
cat("R source files parsed:", length(files), "\n")
