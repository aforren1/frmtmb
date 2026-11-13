# Lane wcheck: do the guard-absent fits of the new tests in
# test-weights-unbounded.R end the same way on every rounding? Each fit
# is repeated on its test data with x (and the gaussian or Gamma
# response) scaled by 1 + 2e-15 * k, k = 0..24, and its nlminb code,
# warnings and messages are counted. The refusals are decided from the
# weights and the design alone, so they are checked once per draw too.
#   Rscript dev/wcheck-perturb.R <lib>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
# wu_data() of the test file
wu_data <- function(seed = 20261020) {
  set.seed(seed)
  ng <- 12; nt <- 6; n <- ng * nt
  d <- data.frame(g = factor(rep(seq_len(ng), each = nt)),
                  tt = factor(rep(seq_len(nt), ng)),
                  x = rnorm(n), xs = rnorm(n), w = runif(n, 0.5, 2))
  d$obs <- factor(seq_len(n))
  eta <- 1 + 0.4 * d$x + rnorm(n, 0, 0.5)
  d$y <- rgamma(n, shape = 3, rate = 3 / exp(eta))
  d$ly <- log(d$y)
  d$k <- rpois(n, exp(eta))
  d$lsig <- rep(log(0.7), n)
  d$fac <- factor(rep(c("a", "b"), length.out = n))
  d
}
fit <- function(...) {
  w <- character(); m <- character()
  f <- tryCatch(withCallingHandlers(frm(...), warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  }, message = function(x) {
    m <<- c(m, conditionMessage(x)); invokeRestart("muffleMessage")
  }), error = function(e) e)
  olre <- "gives every observation its own random effect"
  if (inherits(f, "error")) {
    return(if (grepl("has no maximum|Laplace approximation",
                     conditionMessage(f))) "REFUSED" else "ERROR")
  }
  sprintf("code %d, %d other warnings, %d messages", f$opt$convergence,
          sum(!grepl(olre, w, fixed = TRUE)), length(m))
}
d0 <- wu_data()
set.seed(20261111)
d0$ly2 <- d0$ly + stats::rnorm(nrow(d0), 0, 0.5)
d0$k2 <- stats::rpois(nrow(d0), exp(0.3 + 0.2 * d0$x +
                                     stats::rnorm(nrow(d0), 0, 0.6)))
P <- function(...) set_prior(..., dpar = "sigma")
out <- list()
for (k in 0:24) {
  d <- d0
  s <- 1 + 2e-15 * k
  d$x <- d$x * s; d$ly <- d$ly * s; d$ly2 <- d$ly2 * s; d$y <- d$y * s
  dw <- d; dw$w <- d$w * 0.7
  r <- c(
    mv_pois = fit(mvbf(bf(k ~ x + (1 | obs)), bf(k2 ~ x + (1 | obs))),
                  family = poisson(), data = d),
    mv_pois_w = fit(mvbf(bf(k | weights(w) ~ x + (1 | obs)),
                         bf(k2 | weights(w) ~ x + (1 | obs))),
                    family = poisson(), data = d),
    mv_gauss_w07 = fit(mvbf(bf(ly | weights(w) ~ x + (1 | obs)),
                            bf(ly2 | weights(w) ~ x + (1 | obs))),
                       family = gaussian(), data = dw),
    cells_xs = fit(bf(ly | weights(w) ~ x + (1 | obs), sigma ~ 0 + fac + xs),
                   family = gaussian(), data = d,
                   prior = c(P("normal(0, 1)", class = "b", coef = "faca"),
                             P("normal(0, 1)", class = "b", coef = "facb"))),
    gamma_t30 = fit(bf(y | weights(w) ~ x + (1 | obs)),
                    family = Gamma(link = "log"), data = d,
                    prior = set_prior("student_t(30, 0, 2.5)",
                                      class = "shape")),
    mv_refused = fit(mvbf(bf(ly ~ x + (1 | obs)),
                          bf(ly2 | weights(w) ~ x + (1 | obs))),
                     family = gaussian(), data = d),
    cell_refused = fit(bf(ly | weights(w) ~ x + (1 | obs), sigma ~ 0 + fac),
                       family = gaussian(), data = d,
                       prior = P("normal(0, 1)", class = "b",
                                 coef = "faca")))
  out[[k + 1L]] <- r
}
tab <- do.call(rbind, out)
for (j in colnames(tab)) {
  cat(sprintf("%-14s %s\n", j, paste(names(table(tab[, j])),
                                      table(tab[, j]), sep = " x",
                                      collapse = "; ")))
}
cat("DONE\n")
