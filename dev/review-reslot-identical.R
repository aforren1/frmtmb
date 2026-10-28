# Reviewer copy of dev/reslot-identical.R, writing to dev/review-reslot-out/ident.
# Lane reslot: models the change should not touch, fitted on rellib-r14
# and on the lane library: estimates, standard errors, logLik, MAP fits
# under class priors, both default_prior() tables, and frm_sample()'s
# default stack (-log prior at the estimate). Saves one rds per arm;
# `compare` reads both and prints the largest relative difference.
# Usage: Rscript dev/reslot-identical.R before|after|compare. Seed 1.
arm <- commandArgs(trailingOnly = TRUE)[1L]
out <- "C:/Users/adf44/source/r/frmtmb-wt-reslot/dev/review-reslot-out/ident"
if (identical(arm, "compare")) {
  a <- readRDS(file.path(out, "identical-before.rds"))
  b <- readRDS(file.path(out, "identical-after.rds"))
  rel <- function(x, y) {
    if (is.data.frame(x) || is.character(x)) {
      return(if (identical(x, y)) 0 else NA_real_)
    }
    x <- unlist(x)
    y <- unlist(y)
    if (length(x) != length(y) || !identical(names(x), names(y))) {
      return(NA_real_)
    }
    if (!length(x)) return(0)
    ok <- is.finite(x) | is.finite(y)
    if (!identical(is.finite(x), is.finite(y))) return(NA_real_)
    max(c(0, abs(x[ok] - y[ok]) / pmax(1, abs(x[ok]))))
  }
  for (m in names(a)) {
    for (k in names(a[[m]])) {
      r <- rel(a[[m]][[k]], b[[m]][[k]])
      cat(sprintf("%-14s %-12s %s\n", m, k,
                  if (is.na(r)) "DIFFERS" else sprintf("%.3g", r)))
      if (is.na(r) && is.data.frame(a[[m]][[k]])) {
        x <- a[[m]][[k]]
        y <- b[[m]][[k]]
        kx <- do.call(paste, x[, c("class", "coef", "group", "prior")])
        ky <- do.call(paste, y[, c("class", "coef", "group", "prior")])
        cat("    only before:", paste(setdiff(kx, ky), collapse = " | "),
            "\n    only after: ", paste(setdiff(ky, kx), collapse = " | "),
            "\n")
      }
    }
  }
  quit(save = "no")
}
source("C:/Users/adf44/source/r/frmtmb-wt-reslot/dev/reslot-common.R")
suppressPackageStartupMessages(library(frmtmb.sample))
dd <- re_data(1L)
dd$ord <- factor(cut(dd$y, 4, labels = FALSE), ordered = TRUE)
dd$cnt <- rpois(nrow(dd), exp(0.3 + 0.2 * dd$x))
dd$yp <- abs(dd$y) + 0.2
models <- list(
  gauss = list(f = bf(y ~ x), fam = gaussian(),
               pr = set_prior("normal(0, 1)", class = "b") +
                 set_prior("student_t(3, 0, 2)", class = "sigma")),
  pois = list(f = bf(cnt ~ x), fam = poisson(),
              pr = set_prior("normal(0, 2)", class = "Intercept")),
  ar_cov = list(f = bf(y ~ x + ar(t, g, cov = TRUE)), fam = gaussian(),
                pr = set_prior("normal(0.3, 0.3)", class = "ar")),
  cosy = list(f = bf(y ~ x + cosy(t, g)), fam = gaussian(),
              pr = set_prior("beta(2, 2)", class = "cosy")),
  ordflex = list(f = bf(ord ~ x), fam = cumulative(),
                 pr = set_prior("normal(0, 3)", class = "Intercept")),
  student = list(f = bf(y ~ x), fam = student(),
                 pr = set_prior("gamma(2, 0.1)", class = "nu")),
  gammas = list(f = bf(yp ~ x), fam = Gamma("log"),
                pr = set_prior("gamma(0.01, 0.01)", class = "shape")),
  # the random-effect and rescor models: the fits and every class prior
  # are unchanged; the tables lose their class theta rows
  re1 = list(f = bf(y ~ x + (1 | g)), fam = gaussian(),
             pr = set_prior("exponential(1)", class = "sd")),
  re2 = list(f = bf(y ~ x + (1 + x | g)), fam = gaussian(),
             pr = set_prior("lkj(2)", class = "cor") +
               set_prior("exponential(1)", class = "sd", group = "g",
                         coef = "x")),
  rescor = list(f = bf(y ~ x) + bf(y2 ~ x) + set_rescor(TRUE), fam = NULL,
                pr = set_prior("lkj(2)", class = "rescor")),
  ordequi = list(f = bf(ord ~ x), fam = cumulative(threshold = "equidistant"),
                 pr = set_prior("normal(0, 3)", class = "Intercept")))
res <- list()
for (m in names(models)) {
  md <- models[[m]]
  q <- function(e) suppressWarnings(suppressMessages(e))
  fit <- q(frm(md$f, family = md$fam, data = dd))
  map <- q(frm(md$f, family = md$fam, data = dd, prior = md$pr))
  se <- function(x) sqrt(diag(as.matrix(vcov(x, full = TRUE))))
  rp <- q(frmtmb.sample:::sample_resolve_priors(fit, NULL))
  pars <- fit$obj$env$parList(fit$obj$env$last.par.best)
  res[[m]] <- list(
    est = fit$opt$par, se = tryCatch(se(fit), error = function(e) NA),
    loglik = as.numeric(logLik(fit)),
    map_est = map$opt$par, map_obj = map$opt$objective,
    tab_fit = as.data.frame(default_prior(md$f, family = md$fam,
                                          data = dd)),
    tab_sample = as.data.frame(default_prior(md$f, family = md$fam,
                                             data = dd, route = "sample")),
    dflt_nlp = if (!is.null(rp$ri)) {
      frmtmb:::neg_log_prior_fn(rp$ri$entries)(pars)
    } else 0,
    summary_prior = as.data.frame(rp$effective))
  cat(m, "done\n")
}
saveRDS(res, file.path(out, paste0("identical-", arm, ".rds")))
