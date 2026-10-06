# Reviewer round 2: summarise dev/gpby-rev2-nugget-*.tsv, paired by
# design and seed against the 1e-6 arm.
wt <- "C:/Users/adf44/source/r/frmtmb-wt-gpby"
fs <- list.files(file.path(wt, "dev"), "^gpby-rev2-nugget-.*[.]tsv$",
                 full.names = TRUE)
x <- do.call(rbind, lapply(fs, function(f) {
  z <- utils::read.delim(f, comment.char = "", colClasses = "character")
  z <- z[z$design != "DONE", ]
  z[, -1] <- lapply(z[, -1], as.numeric)
  z
}))
ref <- x[x$nugget == 1e-6, c("design", "seed", "logLik", "se_past", "se_in")]
m <- merge(x, ref, by = c("design", "seed"), suffixes = c("", ".r"))
for (dz in unique(m$design)) for (nu in sort(unique(m$nugget),
                                            decreasing = TRUE)) {
  q <- m[m$design == dz & m$nugget == nu, ]
  cat(sprintf(paste0("%-5s nugget %-6s fits %2d | conv 0: %2d | pdHess: %2d",
                     " | max_grad median %.1e max %.1e | logLik - 1e-6's: ",
                     "median %+.1e range [%+.1e, %+.1e] | se_past / 1e-6's",
                     " range [%.4f, %.4f] | se_in ratio range [%.4f, %.4f]\n"),
              dz, format(nu), nrow(q), sum(q$conv == 0, na.rm = TRUE),
              sum(q$pdHess == 1, na.rm = TRUE), median(q$max_grad),
              max(q$max_grad), median(q$logLik - q$logLik.r),
              min(q$logLik - q$logLik.r), max(q$logLik - q$logLik.r),
              min(q$se_past / q$se_past.r), max(q$se_past / q$se_past.r),
              min(q$se_in / q$se_in.r), max(q$se_in / q$se_in.r)))
}
