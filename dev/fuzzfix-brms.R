# brms 2.23.0 on the two fuzz models, same data as frmtmb.
#   Rscript dev/fuzzfix-brms.R <lib> <cum|gam>
# cum: cumulative + weights + (1 | ga/gb); frmtmb quadrature and Laplace.
# gam: Gamma + weights + ar1 latent errors + mo(xo) * z + shape ~ z.
#   brms spells frmtmb's ar1(tt + 0 | g) as ar(time = tt, gr = g), whose
#   latent residuals are the same observation-level AR(1) effects.
#   Run with the fuzz weights (sum(w) > n, no maximum likelihood) and
#   with the weights times 0.7 (sum(w) < n).
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r10"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(c(a[1], base, ul))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages({library(frmtmb); library(brms)})
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"),
    "brms", format(packageVersion("brms")), "\n")
q <- function(x) suppressMessages(suppressWarnings(x))
show <- function(lbl, est, se) {
  cat(sprintf("%-10s %-28s %10.5f (%.5f)\n", lbl, names(est), est, se),
      sep = "")
}
brm_row <- function(fit, pars) {
  s <- posterior::summarise_draws(
    posterior::subset_draws(posterior::as_draws_df(fit), variable = pars),
    "mean", "sd", "rhat")
  for (i in seq_len(nrow(s))) {
    cat(sprintf("%-10s %-28s %10.5f (%.5f) rhat %.3f\n", "brms",
                s$variable[i], s$mean[i], s$sd[i], s$rhat[i]))
  }
}
if (a[2] == "cum") {
  d <- readRDS("dev/fuzzfix-out/data-cum.rds")
  fo <- frmtmb::bf(y | weights(w) ~ 1 + x + (1 | ga/gb))
  fq <- q(frm(fo, family = cumulative(), data = d, quadrature = TRUE))
  fl <- q(frm(fo, family = cumulative(), data = d))
  for (f in list(quad = fq, laplace = fl)) {
    sm <- q(summary(f))
    print(sm$coefficients %||% sm)
  }
  cat("frmtmb logLik quad", as.numeric(logLik(fq)), "laplace",
      as.numeric(logLik(fl)), "\n")
  print(q(VarCorr(fq)))
  print(q(VarCorr(fl)))
  bfit <- brm(y | weights(w) ~ 1 + x + (1 | ga/gb), family = cumulative(),
              data = d, chains = 4, iter = 2000, seed = 20261020,
              refresh = 0, silent = 2)
  print(summary(bfit))
} else {
  d <- readRDS("dev/fuzzfix-out/data-gam.rds")
  G <- Gamma(link = "log")
  bfo <- brms::bf(y | weights(w) ~ 1 + x + mo(xo) * z +
                    ar(time = tt, gr = g), shape ~ 1 + z)
  d7 <- d; d7$w <- d$w * 0.7
  cat("sum(w) - n:", sum(d$w) - nrow(d), " scaled:", sum(d7$w) - nrow(d7),
      "\n")
  ff <- tryCatch(q(frm(frmtmb::bf(y | weights(w) ~ 1 + x + mo(xo) * z +
                            ar1(tt + 0 | g), shape ~ 1 + z),
                       family = G, data = d)), error = function(e) e)
  cat("frmtmb, fuzz weights:", if (inherits(ff, "error"))
    conditionMessage(ff) else "fitted", "\n")
  f7 <- q(frm(frmtmb::bf(y | weights(w) ~ 1 + x + mo(xo) * z + ar1(tt + 0 | g),
                 shape ~ 1 + z), family = G, data = d7))
  print(q(summary(f7)))
  cat("frmtmb logLik scaled:", as.numeric(logLik(f7)), "\n")
  b7 <- brm(bfo, family = G, data = d7, chains = 4, iter = 2000,
            seed = 20261020, refresh = 0, silent = 2,
            control = list(adapt_delta = 0.95))
  print(summary(b7))
  bw <- update(b7, newdata = d, seed = 20261020, refresh = 0, silent = 2)
  cat("brms with the fuzz weights:\n")
  print(summary(bw))
}
cat("END\n")
