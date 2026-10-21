# After the release review, finding 1: where does a censored model with
# an observation-level effect run away? The check skips cens() on the
# build this runs on, so every fit here is the unguarded objective.
#   Rscript dev/rel073-wcens-probe.R <lib>   (log dev/rel073-log/wcens-probe-<tag>.txt)
# Gaussian, 80 rows, (1 | obs), the top 16 rows right censored at the
# 80th percentile with weight 1, and the 64 uncensored rows at weight
# 1 + E / 64, so the uncensored rows' excess sum(w - 1) is E. Laplace
# and quadrature (the exact integral of the one block). Then
# asym_laplace, unweighted, with 0, 8, 40 and 72 of the 80 rows right
# censored, under Laplace.
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib: frmtmb", format(packageVersion("frmtmb")), dirname(find.package("frmtmb")),
    "\n")
run <- function(label, expr) {
  w <- character()
  r <- tryCatch(withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage")),
  error = function(e) structure(conditionMessage(e), class = "err"))
  if (inherits(r, "err")) {
    cat(sprintf("%-44s ERROR: %s\n", label,
                substr(gsub("[[:space:]]+", " ", r), 1, 110)))
  } else {
    cat(sprintf("%-44s logLik %12.4f  log sigma %8.3f  conv %d  warnings %d\n",
                label, as.numeric(logLik(r)), r$estimates$betad[1],
                r$opt$convergence, length(w)))
  }
}
set.seed(20261021)
n <- 80
d <- data.frame(x = rnorm(n), obs = factor(seq_len(n)))
d$y <- 1 + 0.5 * d$x + rnorm(n)
q <- quantile(d$y, 0.8)
d$cc <- ifelse(d$y > q, "right", "none")
d$yc <- pmin(d$y, q)
nu <- sum(d$cc == "none")
for (E in c(0, 2, 4, 8, 16, 32)) {
  d$w <- ifelse(d$cc == "none", 1 + E / nu, 1)
  for (quad in c(FALSE, TRUE)) {
    run(sprintf("gaussian E %2g %s", E, if (quad) "quadrature" else "Laplace"),
        frm(bf(yc | weights(w) + cens(cc) ~ x + (1 | obs)),
            family = gaussian(), data = d, quadrature = quad))
  }
}
d$ya <- 1 + 0.5 * d$x + rnorm(n)
for (k in c(0, 8, 40, 72)) {
  o <- order(d$ya, decreasing = TRUE)
  ca <- rep("none", n)
  if (k) ca[o[seq_len(k)]] <- "right"
  qa <- if (k) d$ya[o[k]] else Inf
  d$ca <- ca
  d$yac <- pmin(d$ya, qa)
  run(sprintf("asym_laplace, %2d of 80 right censored", k),
      frm(yac | cens(ca) ~ x + (1 | obs), family = asym_laplace(), data = d))
}
cat("END\n")
