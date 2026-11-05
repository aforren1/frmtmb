# car() and sar() (frmtmb.spatial) on families with no residual, against brms
# 2.23.0: what brms builds, what frm() does, L1 where both build, and what each
# does with the field under re_formula = NA. Seed 167.
setwd("C:/Users/adf44/source/r/frmtmb-wt-latres")
.libPaths(c("C:/Users/adf44/source/r/wt-latres-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-latres/dev/stan-cache")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.spatial)
  library(rstan)
})
hel <- new.env(parent = asNamespace("frmtmb.spatial"))
hel$`%||%` <- function(a, b) if (is.null(a)) b else a
sys.source("dev/equiv-lib.R", envir = hel)
# car() rows of the translator live with the extension
sys.source("extensions/frmtmb.spatial/tests/testthat/helper-brms.R", envir = hel)
set.seed(167)
nl <- 12L
M <- matrix(0, nl, nl)
for (i in seq_len(nl - 1L)) M[i, i + 1L] <- M[i + 1L, i] <- 1
M[1, nl] <- M[nl, 1] <- 1
locs <- paste0("L", seq_len(nl))
dimnames(M) <- list(locs, locs)
d <- data.frame(g = factor(rep(locs, each = 8L), levels = locs))
n <- nrow(d)
d$x <- rnorm(n)
u <- rnorm(nl, 0, 0.5)
eta <- 0.8 + 0.4 * d$x + u[as.integer(d$g)]
d$count <- rpois(n, exp(eta))
d$bin <- rbinom(n, 1, plogis(eta - 0.8))
d$pos <- rgamma(n, 4, 4 / exp(eta))
d$w <- M[1, ] [as.integer(d$g)]
Wn <- diag(n) * 0
for (i in seq_len(n - 1L)) Wn[i, i + 1L] <- Wn[i + 1L, i] <- 1
Wn <- Wn / rowSums(Wn)

cases <- list(
  list("count ~ x + car(M, gr = g)", "poisson()"),
  list("count ~ x + car(M, gr = g, type = 'icar')", "poisson()"),
  list("count ~ x + car(M, gr = g, type = 'esicar')", "poisson()"),
  list("count ~ x + car(M, gr = g, type = 'bym2')", "poisson()"),
  list("bin ~ x + car(M, gr = g)", "bernoulli()"),
  list("pos ~ x + car(M, gr = g)", "Gamma(link = 'log')"),
  list("count ~ x + sar(Wn)", "poisson()"),
  list("count ~ x + sar(Wn, type = 'error')", "poisson()"),
  list("bin ~ x + sar(Wn)", "bernoulli()"),
  list("pos ~ x + sar(Wn)", "Gamma(link = 'log')"))
d2 <- list(M = M, Wn = Wn)
for (cs in cases) {
  cat("\n=====", cs[[1]], cs[[2]], "\n")
  f <- as.formula(cs[[1]])
  fb <- hel$eq_eval(cs[[2]], "brms")
  ff <- hel$eq_eval(cs[[2]], "frmtmb")
  bform <- brms::bf(f)
  code <- tryCatch({
    prior <- hel$brms_flat_prior(bform, data = d, family = fb, data2 = d2)
    brms::make_stancode(bform, data = d, family = fb, prior = prior, data2 = d2)
  }, error = function(e) e)
  cat("brms:", if (inherits(code, "error")) conditionMessage(code) else "builds",
      "\n")
  fit <- tryCatch(suppressWarnings(frm(frmtmb::bf(f), family = ff, data = d,
                                       data2 = d2)),
                  error = function(e) e)
  cat("frm: ", if (inherits(fit, "error")) conditionMessage(fit) else {
    sprintf("fits logLik %.6f conv %d", logLik(fit), fit$opt$convergence)
  }, "\n")
  if (inherits(code, "error") || inherits(fit, "error")) next
  l1 <- tryCatch(hel$eq_l1(list(data = d, brms_args = list(data2 = d2)), bform,
                           fb, d, fit, prior, joint = TRUE),
                 error = function(e) e)
  if (inherits(l1, "error")) {
    cat("L1 could not be read:", conditionMessage(l1), "\n")
    next
  }
  cat(sprintf("L1 const %.4g grad %.3g\n", l1$at$const, l1$at$grad))
  bfit <- tryCatch(suppressWarnings(hel$brms_fixed_fit(bform, fb, d, fit,
                                                       data2 = d2)),
                   error = function(e) e)
  if (inherits(bfit, "error")) next
  a <- as.numeric(fitted(fit, re_formula = NA)[, 1L])
  b <- as.numeric(brms::posterior_epred(bfit, re_formula = NA)[1L, ])
  cat("fitted(re_formula = NA) max rel diff from brms",
      max(abs(a - b) / pmax(1, abs(a), abs(b))), "\n")
  a <- as.numeric(fitted(fit)[, 1L])
  b <- as.numeric(brms::posterior_epred(bfit)[1L, ])
  cat("fitted() max rel diff from brms",
      max(abs(a - b) / pmax(1, abs(a), abs(b))), "\n")
}
