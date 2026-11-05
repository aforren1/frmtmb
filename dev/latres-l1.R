# The L1 identity of every latent-residual term x family class against brms
# 2.23.0's log_prob, through the survey's instrument (dev/equiv-lib.R,
# eq_l1(): the joint at the conditional modes plus the map's log-Jacobian at
# the estimate and at two points stepped off it), and the Laplace value
# against one computed from brms's own Hessian in z (a second, chart-free
# measure: dev/latres-laplace below).
#
# Usage: Rscript dev/latres-l1.R <worker> <workers> <outfile>
# Data seed 167. Each row: term, family, const at the three points, inner
# gradient, logLik, and the difference from the Laplace value on brms's
# program.
args <- commandArgs(TRUE)
worker <- as.integer(args[[1]])
workers <- as.integer(args[[2]])
outf <- args[[3]]
setwd("C:/Users/adf44/source/r/frmtmb-wt-latres")
.libPaths(c("C:/Users/adf44/source/r/wt-latres-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-latres/dev/stan-cache")
suppressPackageStartupMessages({
  library(frmtmb)
  library(rstan)
})
# The instrument reads the objective, which frm_sample() reads too, on the forms
# frm() refuses for its Laplace value.
assignInNamespace("latres_refuse_laplace", function(spec, frame) NULL, "frmtmb")
assignInNamespace("latres_refuse_margin", function(fit) NULL, "frmtmb")
`%||%` <- function(a, b) if (is.null(a)) b else a
# the helpers read core's internals, as the test environment gives them
hel <- new.env(parent = asNamespace("frmtmb"))
hel$`%||%` <- `%||%`
sys.source("dev/equiv-lib.R", envir = hel)
eq_l1 <- hel$eq_l1
eq_step_optimizer <- hel$eq_step_optimizer
eq_eval <- hel$eq_eval
eq_try <- hel$eq_try
brms_flat_prior <- hel$brms_flat_prior
brms_stan_model <- hel$brms_stan_model

terms <- c(ar1t = "ar(t, g, p = 1, cov = TRUE)", ar1f = "ar(t, g, p = 1)",
           ar2f = "ar(t, g, p = 2)", ma1t = "ma(t, g, q = 1, cov = TRUE)",
           arma1t = "arma(t, g, p = 1, q = 1, cov = TRUE)",
           cosy = "cosy(t, g)", unstr = "unstr(t, g)")
if (length(args) > 3L) terms <- terms[strsplit(args[[4]], ",")[[1]]]
fams <- list(
  c("poisson()", "count"), c("negbinomial()", "count"),
  c("zero_inflated_poisson()", "count"), c("hurdle_poisson()", "count"),
  c("bernoulli()", "bin"), c("binomial()", "trials"),
  c("Gamma(link = 'log')", "pos"), c("Gamma()", "pos"),
  c("lognormal()", "pos"), c("Beta()", "unit"),
  c("cumulative()", "ord"))
if (length(args) > 4L) {
  keep <- strsplit(args[[5]], ",")[[1]]
  fams <- Filter(function(f) f[[1]] %in% keep, fams)
}

mk <- function() {
  set.seed(167)
  d <- expand.grid(t = 1:5, g = factor(1:30))
  d <- d[order(d$g, d$t), ]
  n <- nrow(d)
  d$x <- rnorm(n)
  e <- as.vector(replicate(30, arima.sim(list(ar = 0.5), 5, sd = 0.5)))
  eta <- 0.8 + 0.4 * d$x + e
  d$count <- rpois(n, exp(eta))
  d$bin <- rbinom(n, 1, plogis(eta - 0.8))
  d$trials <- rbinom(n, 10, plogis(eta - 0.8))
  d$tr <- 10L
  d$pos <- rgamma(n, 4, 4 / exp(eta))
  d$unit <- rbeta(n, 8 * plogis(eta - 0.8), 8 * (1 - plogis(eta - 0.8)))
  d$ord <- factor(findInterval(eta + rlogis(n), c(0, 1, 2)) + 1L,
                  ordered = TRUE)
  d
}
d <- mk()
lhs <- function(resp) if (resp == "trials") "trials | trials(tr)" else resp

grid <- expand.grid(term = names(terms), fam = seq_along(fams),
                    stringsAsFactors = FALSE)
grid <- grid[(seq_len(nrow(grid)) - 1L) %% workers == worker - 1L, ]
if (!file.exists(outf)) {
  cat(paste(c("term", "family", "const", "const_up", "const_dn", "grad",
              "loglik", "laplace_gap", "conv", "note"), collapse = "\t"),
      "\n", sep = "", file = outf)
}

# The Laplace value on brms's program: the joint at z* plus the normal
# integral over the inner block, by a finite difference of Stan's gradient
# in z. It reads no frmtmb chart at all.
stan_laplace <- function(sf, pars, idx, lp) {
  u <- rstan::unconstrain_pars(sf, pars)
  k <- length(idx)
  H <- matrix(0, k, k)
  h <- 1e-4
  for (i in seq_len(k)) {
    up <- u
    dn <- u
    up[idx[i]] <- up[idx[i]] + h
    dn[idx[i]] <- dn[idx[i]] - h
    H[i, ] <- (rstan::grad_log_prob(sf, up, adjust_transform = FALSE)[idx] -
                 rstan::grad_log_prob(sf, dn, adjust_transform = FALSE)[idx]) /
      (2 * h)
  }
  H <- (H + t(H)) / 2
  lp + 0.5 * k * log(2 * pi) -
    0.5 * as.numeric(determinant(-H, logarithm = TRUE)$modulus)
}

for (i in seq_len(nrow(grid))) {
  tm <- grid$term[i]
  fm <- fams[[grid$fam[i]]]
  ftxt <- paste0(lhs(fm[[2]]), " ~ x + ", terms[[tm]])
  rec <- c(tm, fm[[1]], rep(NA, 7), "")
  res <- tryCatch({
    fam_f <- eq_eval(fm[[1]], "frmtmb")
    fam_b <- eq_eval(fm[[1]], "brms")
    fit <- suppressWarnings(frm(frmtmb::bf(as.formula(ftxt)), family = fam_f,
                                data = d))
    bform <- brms::bf(as.formula(ftxt))
    prior <- brms_flat_prior(bform, data = d, family = fam_b)
    refit <- function(k) {
      suppressWarnings(frm(frmtmb::bf(as.formula(ftxt)), family = fam_f,
                           data = d,
                           control = frmtmb_control(optimizer = eq_step_optimizer(c(1, -1)[k]),
                               restarts = 0, check_se = "ignore")))
    }
    spec <- list(data = d)
    l1 <- eq_l1(spec, bform, fam_b, d, fit, prior, joint = TRUE, refit = refit)
    pars <- hel$stan_pars_from_fit(fit, l1$sdat, l1$code, l1$rtab)
    idx <- hel$brms_inner_index(l1$sf, pars)
    lap <- stan_laplace(l1$sf, pars, idx, l1$at$lp)
    ex <- function(k) {
      if (length(l1$extra) >= k && !is.null(l1$extra[[k]]$const)) {
        l1$extra[[k]]$const
      } else NA_real_
    }
    list(l1$at$const, ex(1L), ex(2L), l1$at$grad, as.numeric(logLik(fit)),
         lap - as.numeric(logLik(fit)), fit$opt$convergence,
         paste(vapply(l1$extra, function(e) e$error %||% "", ""),
               collapse = ";"))
  }, error = function(e) e)
  if (inherits(res, "error")) {
    rec[10] <- gsub("[\t\r\n]+", " ", conditionMessage(res))
  } else {
    rec[3:10] <- vapply(res, function(v) {
      if (is.numeric(v)) sprintf("%.10g", v) else as.character(v)
    }, "")
  }
  cat(paste(rec, collapse = "\t"), "\n", sep = "", file = outf, append = TRUE)
}
