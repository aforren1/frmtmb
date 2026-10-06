# Identifiability of ordinal mixtures: label switching, which component
# the start values give which class, and how often the default start
# reaches the best of 12 starts. Replicate seeds 1..R, n = 500.
# Usage: Rscript dev/ordmix-ident.R <none|mu> [R]
# Output on stdout; the driver writes dev/ordmix-log-ident-<mode>.txt.
args <- commandArgs(TRUE)
mode <- if (length(args)) args[1] else "none"
R <- if (length(args) >= 2) as.integer(args[2]) else 20L
.libPaths(c("C:/Users/adf44/source/r/wt-ordmix-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("mode", mode, "R", R, "lib", find.package("frmtmb"), "\n")
gen <- function(seed, n = 500) {
  set.seed(seed)
  x <- rnorm(n)
  cls <- rbinom(n, 1, 0.4)
  lat <- if (mode == "none") {
    # class 1 sits higher on the latent scale: thresholds of its own
    ifelse(cls == 1, 1.5 * x + 1, -0.8 * x - 1) + rlogis(n)
  } else {
    # one set of thresholds, the classes apart in their slopes
    ifelse(cls == 1, 2 * x, -1.5 * x) + rlogis(n)
  }
  data.frame(x = x, cls = cls,
             y = 1L + (lat > -1.5) + (lat > 0) + (lat > 1.5))
}
fam <- mixture(cumulative(), cumulative(), order = mode)
quiet <- function(e) suppressWarnings(suppressMessages(e))
out <- NULL
for (s in seq_len(R)) {
  d <- gen(s)
  fit <- quiet(frm(bf(y ~ x), family = fam, data = d))
  ll0 <- as.numeric(logLik(fit))
  est <- fit$estimates
  # label switching is an identity of the likelihood: the components
  # swapped, and theta1 negated, give the same value
  p <- fit$opt$par
  nm <- names(p)
  sw <- p
  sw[nm == "beta"] <- rev(p[nm == "beta"])
  sw[nm == "betad"] <- -p[nm == "betad"]
  if (mode == "none") {
    i1 <- which(nm == "tau_raw1")
    i2 <- which(nm == "tau_raw2")
    sw[i1] <- p[i2]
    sw[i2] <- p[i1]
  }
  swap_diff <- abs(fit$obj$fn(sw) - fit$obj$fn(p))
  # which class component 1 took: the one whose slope it is nearer
  b <- est$beta
  truth <- if (mode == "none") c(1.5, -0.8) else c(2, -1.5)
  comp1_is_class1 <- abs(b[[1]] - truth[1]) + abs(b[[2]] - truth[2]) <
    abs(b[[1]] - truth[2]) + abs(b[[2]] - truth[1])
  # 11 more starts around the template, seed 1000 + s
  set.seed(1000 + s)
  best <- ll0
  for (j in 1:11) {
    st <- lapply(fit$frame$par_template, function(v) {
      if (!length(v)) return(v)
      v + rnorm(length(v), 0, 1)
    })
    st <- st[c("beta", "betad", grep("^tau_raw", names(st), value = TRUE))]
    # the mapped disc entries keep their value: start is on the template
    f2 <- tryCatch(quiet(frm(bf(y ~ x), family = fam, data = d,
                             start = st)),
                   error = function(e) NULL)
    if (!is.null(f2)) best <- max(best, as.numeric(logLik(f2)))
  }
  row <- data.frame(seed = s, ll_default = ll0, ll_best12 = best,
                    gap = best - ll0, swap_diff = swap_diff,
                    comp1_is_class1 = comp1_is_class1,
                    b1 = b[[1]], b2 = b[[2]],
                    theta1 = stats::plogis(est$betad[["theta1_(Intercept)"]]))
  out <- rbind(out, row)
  cat(sprintf("REP seed=%d ll=%.6f best12=%.6f gap=%.3g swap=%.3g comp1_class1=%s b=(%.3f, %.3f) theta1=%.3f\n",
              s, ll0, best, best - ll0, swap_diff, comp1_is_class1,
              b[[1]], b[[2]], row$theta1))
}
cat(sprintf("SUMMARY mode=%s R=%d default_at_best(gap<1e-4)=%d max_gap=%.4g comp1_is_class1=%d max_swap_diff=%.3g\n",
            mode, nrow(out), sum(out$gap < 1e-4), max(out$gap),
            sum(out$comp1_is_class1), max(out$swap_diff)))
