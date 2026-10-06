# Summarise dev/gpby-cov/<arm>-*.tsv into the block the findings paste.
for (arm in c("base", "lane")) {
  fs <- list.files("dev/gpby-cov", paste0("^", arm, "-.*[.]tsv$"),
                   full.names = TRUE)
  d <- do.call(rbind, lapply(fs, utils::read.delim))
  d <- d[order(d$seed), ]
  libs <- unique(unlist(lapply(paste0(fs, ".lib"), function(f) {
    grep("^lib:", readLines(f), value = TRUE)
  })))
  ok <- d$ok == 1
  k <- sum(d$sim_cover[ok])
  n <- sum(ok)
  p <- k / n
  cat(sprintf("== arm %s (%s)\n", arm, paste(libs, collapse = "; ")))
  cat(sprintf("  seeds %d to %d, %d distinct, %d fitted, %d failed\n",
              min(d$seed), max(d$seed), length(unique(d$seed)), n,
              sum(!ok)))
  cat(sprintf("  simultaneous whole-curve coverage %d / %d = %.4f, ",
              k, n, p))
  cat(sprintf("binomial mcse %.4f, Wilson 95%% (%.4f, %.4f)\n",
              sqrt(p * (1 - p) / n),
              (p + 1.96^2 / (2 * n) - 1.96 * sqrt(p * (1 - p) / n +
                1.96^2 / (4 * n^2))) / (1 + 1.96^2 / n),
              (p + 1.96^2 / (2 * n) + 1.96 * sqrt(p * (1 - p) / n +
                1.96^2 / (4 * n^2))) / (1 + 1.96^2 / n)))
  cat(sprintf("  pointwise per-point coverage %.4f\n",
              mean(d$pt_cover[ok])))
  cat(sprintf("  critical value: median %.4f, range (%.4f, %.4f)\n",
              stats::median(d$crit[ok]), min(d$crit[ok]), max(d$crit[ok])))
}
b <- do.call(rbind, lapply(list.files("dev/gpby-cov", "^base-.*[.]tsv$",
                                      full.names = TRUE), utils::read.delim))
l <- do.call(rbind, lapply(list.files("dev/gpby-cov", "^lane-.*[.]tsv$",
                                      full.names = TRUE), utils::read.delim))
m <- merge(b, l, by = "seed", suffixes = c(".base", ".lane"))
cat(sprintf(paste0("== paired on %d seeds: crit lane / base median %.4f, ",
                   "range (%.4f, %.4f); covered by lane only %d, by base ",
                   "only %d; pointwise identical on %d\n"),
            nrow(m), stats::median(m$crit.lane / m$crit.base),
            min(m$crit.lane / m$crit.base), max(m$crit.lane / m$crit.base),
            sum(m$sim_cover.lane == 1 & m$sim_cover.base == 0),
            sum(m$sim_cover.lane == 0 & m$sim_cover.base == 1),
            sum(m$pt_cover.lane == m$pt_cover.base)))
