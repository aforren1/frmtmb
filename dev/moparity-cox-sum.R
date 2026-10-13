# dev/moparity-cox-study.R's two arms, side by side.
a <- readRDS("dev/moparity-log/cox-study-base.rds")
b <- readRDS("dev/moparity-log/cox-study-lane.rds")
for (nm in c("base", "lane")) {
  r <- if (nm == "base") a else b
  cat(sprintf("%s: seeds %d (distinct %d) | gap > 1e-6 %d, > 1e-4 %d, max %.3g | codes %s | warned %d | weight < 1e-6 at the fit %d, at the polished max %d | evaluations total %d median %.0f\n",
              nm, nrow(r), length(unique(r$seed)), sum(r$gap > 1e-6),
              sum(r$gap > 1e-4), max(r$gap),
              paste(names(table(r$code)), table(r$code), sep = "=",
                    collapse = " "),
              sum(nzchar(r$warn)), sum(r$minw < 1e-6),
              sum(r$polish_minw < 1e-6), sum(r$evals), median(r$evals)))
}
d <- b$ll - a$ll
cat(sprintf("paired lane - base logLik: > 1e-6 %d, < -1e-6 %d, max gain %.3g, max loss %.3g\n",
            sum(d > 1e-6), sum(d < -1e-6), max(d), min(d)))
cat("base seeds with a warning:", paste(a$seed[nzchar(a$warn)], collapse = " "),
    "\nbase seeds short by more than 1e-6:", paste(a$seed[a$gap > 1e-6], collapse = " "), "\n")
