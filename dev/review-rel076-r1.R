# Review: is release fix R1's ?set_prior sentence true? Under
# bf(y ~ x + ar(t, g, cov = TRUE), sigma ~ x), evaluate frmtmb's prior
# penalty (resolve_prior_input() then neg_log_prior_fn(), the functions
# frm() and frm_sample() use) at internal parameter vectors that differ
# only in the ar parameter, for a prior on class Intercept dpar sigma
# and on class b dpar sigma. If the intercept prior is read at the log
# marginal sd, its value does not move with ar; brms's (at the log
# innovation sd, log marginal + log c) would.
#   Rscript dev/review-rel076-r1.R > dev/review-rel076-out/r1.txt
.libPaths(c("C:/Users/adf44/source/r/rellib-r14",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(6)
ng <- 12; k <- 6
d <- data.frame(t = rep(seq_len(k), ng), g = rep(seq_len(ng), each = k),
                x = rnorm(ng * k))
d$y <- 1 + 0.5 * d$x + rnorm(nrow(d))
fo <- bf(y ~ x + ar(t, g, cov = TRUE), sigma ~ x)
probe <- function(pr, lab) {
  fit <- frm(fo, data = d, prior = pr, dry_run = "objective")
  ri <- frmtmb:::resolve_prior_input(fit, pr)
  f <- frmtmb:::neg_log_prior_fn(ri$entries)
  e <- ri$entries[[1]]
  cat("\n--", lab, ": entries", length(ri$entries), "comp", e$comp, "idx",
      e$idx, "innov", if (is.null(e$innov)) "no" else "yes", "\n")
  pars <- fit$par_list %||% fit$obj$env$parList()
  cat("  par blocks:", paste(names(pars), lengths(pars), collapse = "; "),
      "\n")
  out <- NULL
  for (th in c(-1.5, 0, 1.5)) {
    p <- pars
    p$thetaac[] <- th
    for (v in c(-0.5, 0.4)) {
      p[[e$comp]][e$idx] <- v
      out <- rbind(out, data.frame(theta_ar = th, ar = th / sqrt(1 + th^2),
                                   value = v, nlp = f(p)))
    }
  }
  print(out, digits = 10)
  invisible(out)
}
`%||%` <- function(a, b) if (is.null(a)) b else a
a <- probe(set_prior("normal(0, 1)", class = "Intercept", dpar = "sigma"),
           "Intercept dpar sigma")
b <- probe(set_prior("normal(0, 1)", class = "b", dpar = "sigma"),
           "b dpar sigma")
# what brms's density would be at the same internal point: the intercept
# shifted by log c = -0.5 log(1 + theta^2) (innov_logc()), read at the
# centered intercept as frmtmb's entry is
lc <- -0.5 * log(1 + a$theta_ar^2)
cat("\nIntercept prior: nlp spread over ar at a fixed value:",
    max(tapply(a$nlp, a$value, function(z) diff(range(z)))), "\n")
cat("b prior: nlp spread over ar at a fixed value:",
    max(tapply(b$nlp, b$value, function(z) diff(range(z)))), "\n")
