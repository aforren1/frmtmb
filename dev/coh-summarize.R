## Summarize the item 2.6 sweep: coverage AND width, per rung, per arm.
##
## Coverage alone cannot say WHY a rung misses, so this splits the two
## explanations the row is about:
##
##   bias   (mean estimate - truth) / sd(estimate), the "wrong place"
##          reading, in units of the estimator's own spread;
##   calib  mean(reported se) / sd(estimate), the "too narrow" reading.
##          A correctly calibrated interval reports 1.0. The sd here is
##          the Monte Carlo estimate of what the estimator actually
##          does, so it IS the reference the run measures for itself.
##
## Predicted coverage from those two alone, under normality, is
##   P(|Z + bias/calib| < 1.96 / calib) after standardizing by the
## reported se, and it is printed beside the observed coverage: if the
## two agree, nothing but the width and the location is at work.
##
## Run: Rscript dev/coh-summarize.R dev/coh-recovery-main.tsv ...

args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) args <- c("dev/coh-recovery-main.tsv")

read_tsv_kv <- function(path) {
  ln <- readLines(path, warn = FALSE)
  ln <- ln[nzchar(trimws(ln))]
  rows <- lapply(ln, function(l) {
    kv <- strsplit(strsplit(l, "\t", fixed = TRUE)[[1L]], "=",
                   fixed = TRUE)
    setNames(trimws(vapply(kv, function(x) paste(x[-1L], collapse = "="),
                           character(1))),
             trimws(vapply(kv, `[`, character(1), 1L)))
  })
  nms <- unique(unlist(lapply(rows, names)))
  m <- do.call(rbind, lapply(rows, function(r) r[nms]))
  colnames(m) <- nms
  d <- as.data.frame(m, stringsAsFactors = FALSE)
  for (cl in c("rep", "seed", "est", "se", "lo", "hi", "width", "sd_id",
               "sd_idcond", "loglik", "conv", "maxgrad", "nbadse",
               "secs")) {
    if (cl %in% names(d)) d[[cl]] <- as.numeric(d[[cl]])
  }
  for (cl in c("ok", "covers", "pdhess")) {
    if (cl %in% names(d)) d[[cl]] <- d[[cl]] == "TRUE"
  }
  d
}

## Wilson score interval: the replicate count is small enough that the
## normal approximation on a proportion near 1 is not honest.
wilson <- function(k, n, conf = 0.95) {
  if (n == 0) return(c(NA, NA))
  z <- stats::qnorm(1 - (1 - conf) / 2)
  p <- k / n
  cen <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
  hw <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
  c(cen - hw, cen + hw)
}

truth <- 0.5
z95 <- stats::qnorm(0.975)

for (path in args) {
  d <- read_tsv_kv(path)
  ## COMPLETE replicates only. A run stopped part way leaves its last
  ## replicate with some rungs and not others, and a rung table built
  ## from different replicate sets is not a paired comparison.
  n_rung <- length(unique(d$rung))
  keep <- names(which(table(d$seed) == n_rung))
  dropped <- length(unique(d$seed)) - length(keep)
  d <- d[as.character(d$seed) %in% keep, ]
  cat("\n==== ", path, ": ", nrow(d), " rows, ", length(keep),
      " complete replicates (", dropped, " partial dropped) ====\n",
      sep = "")
  bad <- d[!d$ok | d$conv != 0 | !d$pdhess, ]
  cat("fits that failed, did not converge or lost the Hessian: ",
      nrow(bad), "\n", sep = "")
  if (nrow(bad)) print(bad[, c("rung", "seed", "conv", "pdhess")])
  cat(sprintf("%-8s %4s %7s %-16s %8s %8s %7s %7s %8s %7s %7s\n",
              "rung", "n", "cover", "wilson95", "mean_w", "sd_est",
              "mean_se", "calib", "bias_sd", "sd_z", "mean_z"))
  full_w <- NA_real_
  for (rg in unique(d$rung)) {
    s <- d[d$rung == rg & d$ok, ]
    n <- nrow(s)
    k <- sum(s$covers)
    ci <- wilson(k, n)
    sd_est <- stats::sd(s$est)
    mean_se <- mean(s$se)
    calib <- mean_se / sd_est
    bias <- (mean(s$est) - truth) / sd_est
    if (identical(rg, "full")) full_w <- mean(s$width)
    ## The standardized error, which is what the interval is made of:
    ## sd(z) = 1 is a calibrated interval and coverage IS P(|z| < 1.96).
    z <- (s$est - truth) / s$se
    cat(sprintf(paste0("%-8s %4d %7.3f (%.3f, %.3f) %8.4f %8.4f %7.4f",
                       " %7.3f %8.3f %7.3f %7.3f\n"),
                rg, n, k / n, ci[1L], ci[2L], mean(s$width), sd_est,
                mean_se, calib, bias, stats::sd(z), mean(z)))
  }
  cat("\nwidth relative to the correct model (full = 1):\n")
  for (rg in unique(d$rung)) {
    s <- d[d$rung == rg & d$ok, ]
    cat(sprintf("  %-8s %6.3f\n", rg, mean(s$width) / full_w))
  }
  cat("\npredicted coverage from location and width alone:\n")
  for (rg in unique(d$rung)) {
    s <- d[d$rung == rg & d$ok, ]
    sd_est <- stats::sd(s$est)
    calib <- mean(s$se) / sd_est
    mu <- (mean(s$est) - truth) / sd_est
    pred <- stats::pnorm(z95 * calib - mu) - stats::pnorm(-z95 * calib - mu)
    cat(sprintf("  %-8s predicted %.3f  observed %.3f\n", rg, pred,
                mean(s$covers)))
  }
  ## Paired, because every rung of a replicate saw the same data.
  if (all(c("id", "full") %in% d$rung)) {
    a <- d[d$rung == "id" & d$ok, ]
    b <- d[d$rung == "full" & d$ok, ]
    m <- merge(a[, c("seed", "width")], b[, c("seed", "width")],
               by = "seed", suffixes = c("_id", "_full"))
    r <- m$width_full / m$width_id
    cat(sprintf("\npaired width ratio full/id over %d replicates: ",
                nrow(m)))
    cat(sprintf("min %.3f median %.3f max %.3f\n", min(r),
                stats::median(r), max(r)))
  }
  ## The variance components, which say whether a rung collapsed into
  ## the rung below it.
  for (rg in unique(d$rung)) {
    s <- d[d$rung == rg & d$ok, ]
    if (all(is.na(s$sd_id)) && all(is.na(s$sd_idcond))) next
    cat(sprintf(paste0("  %-8s sd_id %6.3f  sd_idcond %6.3f",
                       "  (n collapsed < 0.01: %d)\n"),
                rg, mean(s$sd_id, na.rm = TRUE),
                mean(s$sd_idcond, na.rm = TRUE),
                sum(c(s$sd_id, s$sd_idcond) < 0.01, na.rm = TRUE)))
  }
}
