# Lane alapse, backlog 5490: every row right censored, as
# dev/wcheck2-allcens.R (seed 20261021, 80 rows at min(ya)), with and
# without (1 | g4), asym_laplace() against gaussian(). Per fit: code,
# logLik, the largest gradient, every warning and message, and what the
# standard-error analysis and the boundary-stop test read.
#   Rscript dev/alapse-allcens.R <base|lane>
lib_arg <- commandArgs(TRUE)[1]
source("dev/alapse-common.R")
set.seed(20261021)
n <- 80
d <- data.frame(x = rnorm(n), obs = factor(seq_len(n)),
                g4 = factor(rep(1:20, each = 4)))
d$y <- 1 + 0.5 * d$x + rnorm(n)
d$ya <- 1 + 0.5 * d$x + rnorm(n)
d$ca <- "right"
d$yac <- rep(min(d$ya), n)
fos <- list(re = yac | cens(ca) ~ x + (1 | g4), fixed = yac | cens(ca) ~ x)
for (nm in names(fos)) {
  for (fam in list(asym_laplace(), gaussian())) {
    r <- cap(frm(fos[[nm]], family = fam, data = d))
    cat(sprintf("== %s %s\n", fam$family, nm))
    if (inherits(r$v, "error")) {
      cat("  ERROR", conditionMessage(r$v), "\n")
      next
    }
    f <- r$v
    g <- tryCatch(max(abs(f$obj$gr(f$opt$par))), error = function(e) NA)
    cat(sprintf("  code %d (%s) | logLik %.4g | max|grad| %.3g\n",
                f$opt$convergence, f$opt$message, as.numeric(logLik(f)), g))
    cat("  par:", paste(names(f$opt$par), format(f$opt$par, digits = 4),
                        collapse = ", "), "\n")
    for (i in seq_along(r$w)) {
      cat("  WARN [", r$wc[i], "]: ", substr(r$w[i], 1, 400), "\n", sep = "")
    }
    for (i in seq_along(r$m)) {
      cat("  MSG [", r$mc[i], "]: ", substr(r$m[i], 1, 300), "\n", sep = "")
    }
    an <- tryCatch(ns$se_fit_analysis(f, force = TRUE), error = function(e) e)
    if (inherits(an, "error")) {
      cat("  analysis ERROR", conditionMessage(an), "\n")
    } else {
      cat("  analysis lost:", paste(names(an$lost), an$lost, sep = "=",
                                    collapse = ", "), "\n")
    }
    cat("  se_explained:", format(f$cache$se_explained), "| edge probed:",
        tryCatch(ns$se_edge_probed(f), error = function(e) NA), "\n")
    s <- cap(summary(f))
    out <- capture.output(print(s$v))
    i <- grep("without a standard error", out)
    if (length(i)) cat(paste0("  | ", out[i:length(out)]), sep = "\n")
  }
}
cat("DONE\n")
