# Summarize dev/review-impsamp-accuracy.R output.
out <- "C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-out/accuracy2"
f <- list.files(out, "^[^-]*-[0-9]+[.]tsv$", full.names = TRUE)
errs <- character(0)
rows <- list()
for (ff in f) {
  l <- readLines(ff)
  if (grepl("\tERROR\t", l[[1]])) {
    errs <- c(errs, l[[1]])
    next
  }
  rows[[ff]] <- read.delim(ff, stringsAsFactors = FALSE, quote = "",
                           colClasses = c(warn = "character"))
}
r <- do.call(rbind, rows)
r$dll <- r$ll_imp - r$ll_exact
r$loss <- r$ll_exact - r$ll_at_fit_exact
cat("runs:", nrow(r), "errors:", length(errs), "\n")
options(width = 200)
print(format(r[order(r$design, r$seed),
         c("design", "seed", "draws", "rounds", "grad", "ess_min", "mcse", "dll",
           "loss", "int_minus_agh", "max_z_est", "max_z_lap", "newton",
           "se_ratio_min", "se_ratio_max", "n_warn", "tfit")], digits = 3),
      row.names = FALSE)
cat("\nERRORS\n")
cat(errs, sep = "\n")
cat("\nWARNINGS\n")
for (i in which(r$n_warn > 0)) {
  cat(r$design[i], r$seed[i], ":", substr(r$warn[i], 1, 180), "\n")
}
a <- list.files(out, "-away[.]tsv$", full.names = TRUE)
aw <- do.call(rbind, lapply(a, read.delim, stringsAsFactors = FALSE))
cat("\nAWAY points:", nrow(aw), "\n")
cat("max ad_vs_fd:", format(max(aw$ad_vs_fd, na.rm = TRUE)), "\n")
print(format(aw[order(-aw$ad_vs_fd), ][1:10, ], digits = 3), row.names = FALSE)
cat("\nworst z below the exact (fn under exact, false improvement):\n")
print(format(aw[order(aw$z), ][1:15, ], digits = 3), row.names = FALSE)
cat("\nworst z above the exact:\n")
print(format(aw[order(-aw$z), ][1:15, ], digits = 3), row.names = FALSE)
cat("\nnon-finite fn:\n")
print(format(aw[!is.finite(aw$fn), ], digits = 3), row.names = FALSE)
