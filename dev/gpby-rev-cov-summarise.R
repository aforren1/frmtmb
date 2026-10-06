# Reviewer: summarise dev/gpby-rev-cov/*.tsv and compare the plug-in
# arm's rows with the lane's own dev/gpby-cov/*-1-100.tsv seed by seed.
wt <- "C:/Users/adf44/source/r/frmtmb-wt-gpby"
rd <- function(p) utils::read.delim(file.path(wt, p))
cat_arm <- function(arm, mode) {
  fs <- list.files(file.path(wt, "dev/gpby-rev-cov"),
                   paste0("^", arm, "-", mode, "-.*[.]tsv$"))
  do.call(rbind, lapply(file.path("dev/gpby-rev-cov", fs), rd))
}
wil <- function(k, n) {
  p <- k / n; z <- 1.959964
  c((p + z^2 / (2 * n) - z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) /
      (1 + z^2 / n),
    (p + z^2 / (2 * n) + z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) /
      (1 + z^2 / n))
}
for (mode in c("plugin", "oracle")) {
  for (arm in c("base", "lane")) {
    x <- cat_arm(arm, mode)
    if (is.null(x)) next
    k <- sum(x$sim_cover, na.rm = TRUE); n <- sum(x$ok == 1)
    w <- wil(k, n)
    cat(sprintf(paste0("%s %s: seeds %d to %d, %d distinct, %d fitted | ",
                       "whole-curve %d / %d = %.4f, Wilson (%.4f, %.4f), ",
                       "binomial p vs 0.95 %.3g | pointwise %.4f | crit ",
                       "median %.4f%s\n"),
                mode, arm, min(x$seed), max(x$seed), length(unique(x$seed)),
                n, k, n, k / n, w[1], w[2],
                stats::binom.test(k, n, 0.95)$p.value,
                mean(x$pt_cover, na.rm = TRUE), stats::median(x$crit),
                if (mode == "oracle") {
                  sprintf(" | Sigma vs closed form max rel %.3e",
                          max(x$cf_rel, na.rm = TRUE))
                } else ""))
  }
  b <- cat_arm("base", mode); l <- cat_arm("lane", mode)
  if (!is.null(b) && !is.null(l)) {
    m <- merge(b, l, by = "seed", suffixes = c(".b", ".l"))
    cat(sprintf(paste0("%s paired on %d seeds: lane only %d, base only %d, ",
                       "pointwise identical on %d\n"), mode, nrow(m),
                sum(m$sim_cover.l == 1 & m$sim_cover.b == 0),
                sum(m$sim_cover.l == 0 & m$sim_cover.b == 1),
                sum(m$pt_cover.l == m$pt_cover.b)))
  }
}
for (arm in c("base", "lane")) {
  mine <- cat_arm(arm, "plugin")
  theirs <- rd(sprintf("dev/gpby-cov/%s-1-100.tsv", arm))
  m <- merge(mine, theirs, by = "seed", suffixes = c(".rev", ".lane"))
  cat(sprintf(paste0("REPRO %s seeds 1-100: sim_cover equal on %d, crit ",
                       "identical on %d, max |crit diff| %.4g, max |se_max ",
                       "rel diff| %.3g\n"), arm, sum(m$sim_cover.rev ==
                                                  m$sim_cover.lane),
              sum(m$crit.rev == m$crit.lane),
              max(abs(m$crit.rev - m$crit.lane)),
              max(abs(m$se_max.rev / m$se_max.lane - 1))))
}
