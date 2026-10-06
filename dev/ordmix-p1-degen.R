# Punch round 1, B2: the measurements a degenerate-component criterion
# is chosen from. For every fit, the reviewer's label (any |estimate| >
# 30 or a non-finite standard error, dev/ordmix-rev-degen.R) beside the
# per-component statistics of mixture_ord_degeneracy(): the largest
# latent distance |disc (tau - cs - eta)|, a collapsed gap, the smallest
# mean weight.
# Usage: Rscript dev/ordmix-p1-degen.R <config>
#   rev_cum_cum_300 rev_cum_cum_500 rev_cum_sr_300 rev_cum_sr_500
#     (the reviewer's 160 fits: seeds 1..40, its data-generating process)
#   id_<model> (the reviewer's identified designs of
#     dev/ordmix-rev-falsealarm.R, seeds 1..20, n = 500)
# Output: dev/ordmix-p1-degen-log/<config>.txt (the driver redirects)
args <- commandArgs(TRUE)
cfg <- args[1]
.libPaths(c("C:/Users/adf44/source/r/wt-ordmix-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("config", cfg, "lib", find.package("frmtmb"), "\n")
cut4 <- function(lat) 1L + (lat > -1.5) + (lat > 0) + (lat > 1.5)
# the reviewer's identified designs, as dev/ordmix-rev-falsealarm.R
# generates them
gen_id <- function(seed, n = 500) {
  set.seed(seed)
  x <- rnorm(n)
  z <- rnorm(n)
  cls <- rbinom(n, 1, 0.4)
  lat_none <- ifelse(cls == 1, 1.5 * x + 1, -0.8 * x - 1) + rlogis(n)
  lat_mu <- ifelse(cls == 1, 2 * x, -1.5 * x) + rlogis(n)
  lat_pr <- ifelse(cls == 1, 1.5 * x + 1, -0.8 * x - 1) + rnorm(n)
  y <- cut4(lat_none)
  hu <- plogis(ifelse(cls == 1, -2, -0.4) + 0.8 * z)
  yh <- ifelse(runif(n) < hu, 0L, y)
  data.frame(x, z, y, yh, ymu = cut4(lat_mu), ypr = cut4(lat_pr),
             gr = factor(sample(letters[1:2], n, TRUE)))
}
id_spec <- list(
  cum2 = list(f = bf(y ~ x), fam = mixture(cumulative(), cumulative())),
  probit2 = list(f = bf(ypr ~ x),
                 fam = mixture(cumulative("probit"), cumulative("probit"))),
  cum_sr_clog = list(f = bf(y ~ x),
                     fam = mixture(cumulative(), sratio("cloglog"))),
  mu2 = list(f = bf(ymu ~ x),
             fam = mixture(cumulative(), cumulative(), order = "mu")),
  hurdle_huz = list(f = bf(yh ~ x, hu1 ~ z, hu2 ~ z),
                    fam = mixture(hurdle_cumulative(), hurdle_cumulative())),
  sr_acat_th = list(f = bf(y ~ x, theta1 ~ z),
                    fam = mixture(sratio(), acat())),
  gr_cum2 = list(f = bf(y | thres(gr = gr) ~ x),
                 fam = mixture(cumulative(), cumulative())))
rev_label <- function(f) {
  fx <- suppressWarnings(fixef(f))
  any(abs(fx[, "Estimate"]) > 30) || any(!is.finite(fx[, "Est.Error"]))
}
one <- function(s, d, f, fam) {
  w <- character(0)
  fit <- tryCatch(withCallingHandlers(frm(f, family = fam, data = d),
    warning = function(e) {
      w <<- c(w, conditionMessage(e))
      invokeRestart("muffleWarning")
    }), error = function(e) e)
  if (inherits(fit, "error")) {
    cat(sprintf("REP cfg=%s seed=%d ERROR %s\n", cfg, s,
                substr(conditionMessage(fit), 1, 120)))
    return(invisible())
  }
  dg <- frmtmb:::mixture_ord_degeneracy(fit, names(fit$spec$responses)[1])
  gm <- max(abs(fit$obj$gr(fit$opt$par)))
  cat(sprintf(paste0("REP cfg=%s seed=%d rev_degenerate=%s ll=%.4f ",
                     "max_reach=%.4g collapsed=%s min_cover=%.3g min_weight=%.3g ",
                     "gradfinite=%s nwarn=%d\n"),
              cfg, s, rev_label(fit), as.numeric(logLik(fit)),
              max(dg$reach), any(dg$collapsed), min(dg$cover),
              min(dg$weight),
              is.finite(gm), length(w)))
}
if (startsWith(cfg, "rev_")) {
  p <- strsplit(sub("^rev_", "", cfg), "_")[[1]]
  n <- as.integer(p[length(p)])
  fam <- if (identical(p[2], "cum")) mixture(cumulative(), cumulative()) else
    mixture(cumulative(), sratio())
  for (s in 1:40) {
    set.seed(s)
    x <- rnorm(n)
    cls <- rbinom(n, 1, 0.4)
    lat <- ifelse(cls == 1, 1.5 * x + 1, -0.8 * x - 1) + rlogis(n)
    d <- data.frame(x = x, y = cut4(lat))
    one(s, d, bf(y ~ x), fam)
  }
} else {
  sp <- id_spec[[sub("^id_", "", cfg)]]
  for (s in 1:20) one(s, gen_id(s), sp$f, sp$fam)
}
