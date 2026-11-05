# Item 9: is p069's L4 max z of 0.627 a correlation boundary effect? The same
# design (eq_dat("poisson"), y ~ x + (1 + x | g), 10 groups of 20), brms flat
# posterior, Laplace and importance read by the instrument's own eq_l4, at
# seed 69 and at seeds whose ML correlation is off the boundary. If the reading
# is above 0.5 at seeds off the boundary too, "boundary" is not the cause.
# Usage (from the worktree root):
#   Rscript dev/review-impsamp-p069.R <seed> [seed ...]
.libPaths(c("C:/Users/adf44/source/r/wt-impsamp-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache", winslash = "/"),
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({
  library(frmtmb)
  library(testthat)
})
attach(list(me = brms::me), name = "equiv_brms_me", warn.conflicts = FALSE)
E <- new.env(parent = asNamespace("frmtmb"))
for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R", "dev/equiv-flow.R")) {
  sys.source(f, envir = E)
}
seeds <- as.integer(commandArgs(TRUE))
ff <- brms::bf(y ~ x + (1 + x | g))
spec <- list(method = NULL, margs = list())
out <- "dev/review-impsamp-out/p069-flat.tsv"
for (s in seeds) {
  d <- E$eq_dat("poisson", seed = s)
  f2 <- frmtmb::bf(y ~ x + (1 + x | g))
  lap <- suppressWarnings(frm(f2, family = poisson(), data = d))
  imp <- suppressWarnings(tryCatch(
    frm(f2, family = poisson(), data = d, importance = 200L),
    error = function(e) e))
  vc <- frmtmb:::varcorr_matrices(lap)[[1L]]
  cors <- vc[1, 2] / sqrt(vc[1, 1] * vc[2, 2])
  post <- suppressWarnings(E$eq_brms_sample(spec, ff, stats::poisson(), d, E$eq_brms_prior(spec, ff, stats::poisson(), d)))
  rhat <- suppressWarnings(max(brms::rhat(post), na.rm = TRUE))
  ep <- brms::posterior_epred(post)
  pm <- colMeans(ep)
  psd <- apply(ep, 2, stats::sd)
  zl <- (fitted(lap)[, "Estimate"] - pm) / psd
  zi <- if (inherits(imp, "error")) NA else
    (fitted(imp)[, "Estimate"] - pm) / psd
  ps <- as.matrix(post, variable = c("sd_g__Intercept", "sd_g__x",
                                     "cor_g__Intercept__x"))
  worst <- which.max(abs(zl))
  row <- data.frame(seed = s, ml_cor = cors, sd1 = sqrt(vc[1, 1]),
                    sd2 = sqrt(vc[2, 2]), post_cor = mean(ps[, 3]),
                    post_sd1 = mean(ps[, 1]), post_sd2 = mean(ps[, 2]),
                    rhat = rhat, z_lap_max = max(abs(zl)),
                    z_imp_max = max(abs(zi)),
                    z_lap_med = median(abs(zl)), z_imp_med = median(abs(zi)),
                    dfit = max(abs(fitted(lap)[, "Estimate"] -
                                     if (inherits(imp, "error")) NA else
                                       fitted(imp)[, "Estimate"])) ,
                    worst_group = as.character(d$g[worst]),
                    worst_x = d$x[worst],
                    worst_y = d$y[worst],
                    fa = fitted(lap)[worst, "Estimate"],
                    pm = pm[worst], psd = psd[worst],
                    imp_ll = if (inherits(imp, "error")) NA else
                      as.numeric(logLik(imp)),
                    lap_ll = as.numeric(logLik(lap)))
  print(format(row, digits = 4))
  write.table(row, out, sep = "\t", quote = FALSE, row.names = FALSE,
              append = file.exists(out), col.names = !file.exists(out))
}
