# Item 10: "about 0.9 more" by a run, not arithmetic: the summary with the two
# adjudication rows p043 and p169 dropped, on this review's rerun of all pairs.
out <- "dev/review-refuse-out/adj"
dir.create(file.path(out, "equiv-noadj-out"), recursive = TRUE, showWarnings = FALSE)
a <- read.delim("dev/equiv-data/adjudication.tsv", stringsAsFactors = FALSE)
a <- a[!a$id %in% c("p043", "p169"), ]
write.table(a, file.path(out, "adjudication.tsv"), sep = "\t", row.names = FALSE)
x <- readLines("dev/equiv-summary.R")
x <- sub("dev/equiv-data/adjudication.tsv\"",
         paste0(out, "/adjudication.tsv\""), x, fixed = TRUE)
tf <- "dev/review-refuse-summary-noadj.R"
writeLines(x, tf)
fs <- list.files("dev/review-refuse-equiv-out", pattern = "^p[0-9]+[.]tsv$", full.names = TRUE)
file.copy(fs, file.path(out, "equiv-noadj-out"), overwrite = TRUE)
system2("C:/Program Files/R/R-4.6.1/bin/Rscript.exe",
        c(tf, file.path(out, "equiv-noadj-out")))
s <- readLines(file.path(out, "equiv-noadj-out", "summary.md"))
cat(grep("headline", s, value = TRUE), sep = "\n")
print(s[9:13])
