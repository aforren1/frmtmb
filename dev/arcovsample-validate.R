# Lane wt-arcovsample, validation (a) and (c): log_lik() row sums
# against the taped objective at the SAME draw, and loo() on the result.
#
# The row sum against the objective is an IDENTITY by construction (both
# sides compose row_lpdf() over arma_cond_dpars()'s mu), so what is
# reported is the residual at full precision over every draw, not a
# confirmation. The interesting part is the RANDOM-EFFECT correction:
# build_objective() is the bare likelihood plus each block's own prior,
# and only the row product is supposed to be in log_lik().
#
#   Rscript dev/arcovsample-validate.R > dev/arcovsample-log/validate.txt

LIB <- "C:/Users/adf44/source/r/wt-arcovsample-lib"
.libPaths(c(LIB, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))

mk <- "C:/Users/adf44/Documents/.R/Makevars.win"
if (!nzchar(Sys.getenv("R_MAKEVARS_USER")) && file.exists(mk)) {
  Sys.setenv(R_MAKEVARS_USER = mk)
}
stopifnot(utils::packageVersion("tmbstan") >= "1.2.1",
          any(grepl("-std=gnu++17",
                    readLines(tools::makevars_user(), warn = FALSE),
                    fixed = TRUE)))
suppressMessages(library(frmtmb))
suppressMessages(library(frmtmb.sample))
cat("frmtmb        ", format(packageVersion("frmtmb")), "\n")
cat("frmtmb.sample ", format(packageVersion("frmtmb.sample")), "\n")

SEED <- 4021L

# ragged groups with interior gaps, rows shuffled: the sort order and
# the row-counted lags are both exercised
mk_data <- function(seed = SEED, ng = 10L, nt = 10L) {
  set.seed(seed)
  dd <- expand.grid(t = seq_len(nt), g = factor(seq_len(ng)))
  dd <- dd[-c(3L, 17L, 41L, 55L, 72L), ]
  dd$x <- stats::rnorm(nrow(dd))
  u <- stats::rnorm(ng, 0, 0.6)
  dd$y <- 1 + 0.5 * dd$x + u[as.integer(dd$g)] +
    unlist(lapply(split(seq_len(nrow(dd)), dd$g), function(r) {
      as.numeric(stats::arima.sim(list(ar = 0.6, ma = 0.3), length(r),
                                  sd = 0.6))
    }))[order(order(dd$g))]
  dd$w <- stats::runif(nrow(dd), 0.5, 2)
  dd$cc <- sample(c(0, 0, 0, 1, -1), nrow(dd), TRUE)
  dd$yp <- dd$y + 5
  dd[sample(nrow(dd)), ]
}

dd <- mk_data()
cat("N =", nrow(dd), " groups =", length(unique(dd$g)), " seed =", SEED,
    "\n")

# Each random-effect block's own prior contribution at one parameter
# vector: build_objective() carries it and log_lik() must not.
re_prior_of <- function(sh) {
  bks <- sh$frame[["re_blocks"]] %||% list()
  tot <- 0
  for (bk in bks) {
    f <- frmtmb:::covstruct_registry[[bk$covstruct]]$nll
    tot <- tot + as.numeric(f(sh$estimates[["b"]][bk$b_idx],
                              sh$estimates[["theta"]][bk$theta_idx], bk))
  }
  tot
}
`%||%` <- function(a, b) if (is.null(a)) b else a

one <- function(label, form, family, data = dd) {
  ds <- suppressWarnings(suppressMessages(
    frm_sample(form, family = family, data = data,
               chains = 1, iter = 500, refresh = 0, seed = 7)))
  ll <- log_lik(ds)
  idx <- frmtmb.sample:::draws_par_index(ds$fit)
  nd <- nrow(ds$draws)
  res <- numeric(nd)
  for (d in seq_len(nd)) {
    sh <- frmtmb.sample:::draws_fit_at(ds, d, idx)
    nll <- as.numeric(frmtmb::build_objective(sh$frame)(sh$estimates))
    res[d] <- sum(ll[d, ]) - (-nll - re_prior_of(sh))
  }
  lo <- suppressWarnings(loo(ds))
  wa <- suppressWarnings(waic(ds))
  cat("\n---- ", label, " ----\n", sep = "")
  cat("  draws x cols      : ", nd, " x ", ncol(ll), " (N = ",
      nrow(data), ")\n", sep = "")
  cat("  max |row sum - (-nll - re_prior)| : ",
      format(max(abs(res)), digits = 12), "\n", sep = "")
  cat("  max relative residual             : ",
      format(max(abs(res) / abs(rowSums(ll))), digits = 12), "\n",
      sep = "")
  cat("  any non-finite log_lik cell       : ", anyNA(ll) ||
        any(!is.finite(ll)), "\n", sep = "")
  cat("  elpd_loo / se                     : ",
      format(lo$estimates["elpd_loo", "Estimate"], digits = 10), " / ",
      format(lo$estimates["elpd_loo", "SE"], digits = 6), "\n", sep = "")
  cat("  p_loo                             : ",
      format(lo$estimates["p_loo", "Estimate"], digits = 8), "\n",
      sep = "")
  cat("  elpd_waic                         : ",
      format(wa$estimates["elpd_waic", "Estimate"], digits = 10), "\n",
      sep = "")
  cat("  max pareto k                      : ",
      format(max(lo$diagnostics$pareto_k), digits = 4), "\n", sep = "")
  invisible(list(ds = ds, ll = ll, res = res, loo = lo))
}

r <- list()
r$ar1 <- one("gaussian ar(1), no random effect",
             bf(y ~ x + ar(t, g)), gaussian())
r$arma_re <- one("student arma(1,1) + (1 | g)",
                 bf(y ~ x + arma(t, g, p = 1, q = 1) + (1 | g)),
                 student())
r$sigma <- one("gaussian ar(1), sigma ~ x",
               bf(y ~ x + ar(t, g), sigma ~ x), gaussian())
r$weights <- one("gaussian ma(1), weights(w)",
                 bf(y | weights(w) ~ x + ma(t, g)), gaussian())
r$cens <- one("gaussian arma(1,1), cens(cc)",
              bf(y | cens(cc) ~ x + arma(t, g, p = 1, q = 1)), gaussian())
r$trunc <- one("gaussian ma(1), trunc(lb = 0)",
               bf(yp | trunc(lb = 0) ~ x + ma(t, g)), gaussian())
r$ar2 <- one("gaussian ar(2) + (1 | g)",
             bf(y ~ x + ar(t, g, p = 2) + (1 | g)), gaussian())

# The one-step mean is what posterior_epred() already reported, so the
# two surfaces must agree cell for cell. A MEASUREMENT, not an identity:
# posterior_epred goes through frm_linpred() and log_lik through
# arma_cond_dpars(), two different code paths to the same mu.
cat("\n---- log_lik()'s mu against posterior_epred()'s ----\n")
# trunc() is left out on purpose: posterior_epred() of a truncated
# response is the truncated mean, not mu, so the two are not the same
# quantity there and a difference would say nothing about the shift.
for (nm in setdiff(names(r), "trunc")) {
  ds <- r[[nm]]$ds
  fit <- frmtmb.sample:::draws_base_fit(ds)
  rn <- names(fit$spec$responses)[1L]
  idx <- frmtmb.sample:::draws_par_index(fit)
  ep <- posterior_epred(ds)
  mu <- t(vapply(seq_len(nrow(ds$draws)), function(d) {
    sh <- frmtmb.sample:::draws_fit_at(ds, d, idx)
    dpv <- frmtmb::arma_cond_dpars(sh, frmtmb::eval_dpars(sh))
    as.numeric(dpv[[rn]][["mu"]])
  }, numeric(ncol(ep))))
  cat("  ", nm, ": max |mu - epred| = ",
      format(max(abs(mu - ep)), digits = 6), " (epred sd ",
      format(stats::sd(ep), digits = 4), ")\n", sep = "")
}

# cov = TRUE must still refuse, with a message true of it
cat("\n---- cov = TRUE is still refused ----\n")
uc <- frm(bf(y ~ x + ar(t, g, cov = TRUE)), family = gaussian(),
          data = dd, dry_run = "objective")
lab <- c(frmtmb::brms_par_labels(uc), "lp__")
fd <- structure(list(stanfit = NULL,
                     draws = matrix(0, 4L, length(lab),
                                    dimnames = list(NULL, lab)),
                     fit = uc), class = "frmtmb_draws")
cat(conditionMessage(tryCatch(log_lik(fd), error = identity)), "\n")
uu <- frm(bf(y ~ x + unstr(t, g)), family = gaussian(), data = dd,
          dry_run = "objective")
lab <- c(frmtmb::brms_par_labels(uu), "lp__")
fu <- structure(list(stanfit = NULL,
                     draws = matrix(0, 4L, length(lab),
                                    dimnames = list(NULL, lab)),
                     fit = uu), class = "frmtmb_draws")
cat(conditionMessage(tryCatch(log_lik(fu), error = identity)), "\n")

cat("\nDONE\n")
