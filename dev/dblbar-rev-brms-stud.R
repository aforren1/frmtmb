# Reviewer, lane dblbar: brms 2.23.0's Stan code for a numeric || term
# with dist = "student". Output: dev/dblbar-rev-log/brms-stud.txt
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-dblbar"
source(file.path(wt, "dev/dblbar-data.R"))
suppressPackageStartupMessages(library(brms))
dd <- dblbar_data()
sink(file.path(wt, "dev/dblbar-rev-log/brms-stud.txt"), split = TRUE)
for (txt in c("bf(y ~ x + (1 + x || gr(g, dist = 'student')))",
              paste0("bf(y ~ x + (1 | gr(g, dist = 'student')) + ",
                     "(0 + x | gr(g, dist = 'student')))"))) {
  cat("====", txt, "\n")
  sc <- strsplit(stancode(eval(str2lang(txt)), data = dd), "\n")[[1]]
  cat(trimws(grep("udf|dfm|r_[0-9]+_[0-9]+ = ", sc, value = TRUE)),
      sep = "\n")
}
sink()
