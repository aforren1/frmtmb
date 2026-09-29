# Lane wt-arcovsample: localize the student + set_rescor(TRUE)
# disagreement dev/arcovsample-rescor.R found. Is it the ARMA shift, or
# is rescor_row_loglik() already off the objective for a Student-t
# rescor model with NO autocorrelation at all?
#
# ARCOVSAMPLE_REF=true runs it against the base commit instead.
#   Rscript dev/arcovsample-rescor-t.R

LIB <- if (identical(Sys.getenv("ARCOVSAMPLE_REF"), "true")) {
  "C:/Users/adf44/source/r/rellib-r3"
} else {
  "C:/Users/adf44/source/r/wt-arcovsample-lib"
}
.libPaths(c(LIB, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
suppressMessages(library(frmtmb.sample))
cat("library:", LIB, "\n")

set.seed(4023)
ng <- 8L
nt <- 9L
dd <- expand.grid(t = seq_len(nt), g = factor(seq_len(ng)))
dd <- dd[-c(4L, 20L, 33L), ]
dd$x <- stats::rnorm(nrow(dd))
dd$y <- 1 + 0.5 * dd$x + stats::rnorm(nrow(dd))
dd$y2 <- 0.4 + 0.8 * dd$y + stats::rnorm(nrow(dd), 0, 0.7)
dd <- dd[sample(nrow(dd)), ]
rownames(dd) <- NULL

one <- function(label, form, family) {
  ds <- suppressWarnings(suppressMessages(
    frm_sample(form, family = family, data = dd, chains = 1, iter = 400,
               refresh = 0, seed = 9)))
  ll <- log_lik(ds)
  idx <- frmtmb.sample:::draws_par_index(ds$fit)
  nd <- nrow(ds$draws)
  res <- numeric(nd)
  for (d in seq_len(nd)) {
    sh <- frmtmb.sample:::draws_fit_at(ds, d, idx)
    res[d] <- sum(ll[d, ]) +
      as.numeric(frmtmb::build_objective(sh$frame)(sh$estimates))
  }
  nu <- if ("nu" %in% colnames(ds$draws)) ds$draws[, "nu"] else NA_real_
  cat("\n---- ", label, " ----\n", sep = "")
  cat("  max |row sum + nll| : ", format(max(abs(res)), digits = 10),
      "\n", sep = "")
  cat("  draws with |res| > 1e-8 : ", sum(abs(res) > 1e-8), " of ", nd,
      "\n", sep = "")
  cat("  nu range : ", format(range(nu), digits = 6), "\n", sep = " ")
  # is the residual a function of nu?
  if (!all(is.na(nu))) {
    big <- abs(res) > 1e-8
    if (any(big)) {
      cat("  nu of the disagreeing draws : ",
          format(range(nu[big]), digits = 6), "\n", sep = " ")
      cat("  nu of the agreeing draws    : ",
          format(range(nu[!big]), digits = 6), "\n", sep = " ")
    }
  }
  invisible(NULL)
}

one("gaussian rescor, no autocor",
    bf(y ~ x) + bf(y2 ~ x) + set_rescor(TRUE), gaussian())
one("student rescor, NO autocor",
    bf(y ~ x) + bf(y2 ~ x) + set_rescor(TRUE), student())
one("student rescor + ar(1) on both",
    bf(y ~ x + ar(t, g)) + bf(y2 ~ x + ar(t, g)) + set_rescor(TRUE),
    student())
cat("\nDONE\n")
