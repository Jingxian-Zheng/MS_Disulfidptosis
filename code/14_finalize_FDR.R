if (.Platform$OS.type == "windows") suppressWarnings(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
p <- Sys.getenv("REVISION_ROOT", unset = normalizePath(getwd(), winslash = "/", mustWork = TRUE))
d <- read.csv(file.path(p, "tables/S08_all_section_spatial_associations.csv"))
qq <- p.adjust(d$spatial_GLS_P, "BH", n = nrow(d))
cat("Planned tests", nrow(d), "missing GLS", sum(is.na(d$spatial_GLS_P)), "max q difference", max(abs(qq -
    d$spatial_GLS_q_BH), na.rm = TRUE), "
")
print(table(d$fit_status))
d$spatial_GLS_q_BH <- qq
d$raw_Pearson_q_BH <- p.adjust(d$raw_Pearson_P, "BH", n = nrow(d))
write.csv(d, file.path(p, "tables/S08_all_section_spatial_associations.csv"), row.names = FALSE)
a <- d[d$condition == "MS" & d$module == "panel13_no_DSE", ]
print(aggregate(as.integer(a$spatial_GLS_q_BH < 0.05), list(target = a$target), sum))
writeLines(capture.output(sessionInfo()), file.path(p, "audit/final_statistics_sessionInfo.txt"))
