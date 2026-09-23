# Punch round 1, the feasibility question: what would it take to give
# simulate() a `newdata` and let pp_check() pass it through, so that
# pp_check(newdata = ) answers as brms answers?
#
# This is a PROBE, not an implementation. It assembles the draw outside
# the package, from the same pieces frmtmb.sample's
# posterior_predict(newdata = ) already uses, and then asks the two
# ported assertions (brmsfit-methods:675 and :682) whether they would
# hold. Nothing here is installed.
#   Rscript dev/correct-newdata-feasibility.R
source("C:/Users/adf44/source/r/frmtmb-wt-correct/dev/correct-prelude.R")
suppressMessages(library(frmtmb))
source("C:/Users/adf44/source/r/frmtmb-wt-correct/dev/correct-data.R")

# The recipe, copied from frmtmb.sample:::posterior_predict.frmtmb_draws:
# dpars at newdata through frm_linpred(), aterms through
# aterms_for_newdata() when the response is truncated, a refusal when
# the draw is structured, and sim_draw(sim_context(...)) per replicate.
simulate_newdata <- function(fit, newdata, nsim = 10, re_formula = NA) {
  rspec <- frmtmb:::single_response(fit, "simulate()")
  av <- if (has_trunc(rspec)) aterms_for_newdata(rspec, newdata) else list()
  if (sim_is_structured(sim_context(fit, rspec, list(), aterms = av))) {
    stop("structured draw: the structure indexes the fitted rows")
  }
  out <- vector("list", nsim)
  for (s in seq_len(nsim)) {
    dp <- list()
    for (dnm in names(rspec$dpars)) {
      dp[[dnm]] <- as.vector(frm_linpred(fit, newdata = newdata,
                                         dpar = dnm, type = "response",
                                         re_formula = re_formula))
    }
    out[[s]] <- sim_draw(sim_context(fit, rspec, dp, aterms = av,
                                     n = length(dp[[1L]]),
                                     extra = fit_extras(fit)))
  }
  do.call(cbind, out)
}

d <- correct_data_gauss()
fit <- frm(bf(y ~ x + (1 | g)) + gaussian(), data = d)
nd <- d[1:12, ]

cat("== the pieces the fit method would need ==\n")
for (nm in c("frm_linpred", "aterms_for_newdata", "sim_context", "sim_draw",
             "sim_is_structured", "fit_extras", "has_trunc")) {
  where <- if (nm %in% getNamespaceExports("frmtmb")) "exported" else
    "internal"
  cat(sprintf("  %-20s %s\n", nm, where))
}

cat("\n== a draw at newdata, outside the package ==\n")
set.seed(1)
yrep <- t(simulate_newdata(fit, nd, nsim = 10))
cat("yrep is", nrow(yrep), "draws x", ncol(yrep), "rows; y has",
    nrow(nd), "\n")
p1 <- bayesplot::ppc_dens_overlay(nd$y, yrep)
cat("ported :675, expect_ggplot(pp_check(fit, newdata = 10 rows)):",
    inherits(p1, "ggplot"), "\n")
p2 <- bayesplot::ppc_violin_grouped(nd$y, yrep, group = nd$g)
cat("ported :682, violin_grouped with group read off newdata:",
    inherits(p2, "ggplot"), "\n")
invisible(ggplot2::ggplot_build(p1))
invisible(ggplot2::ggplot_build(p2))

cat("\n== where it would refuse ==\n")
# a structured family: the group-level mixture, whose draw walks the
# fitted groups
set.seed(4)
dg <- data.frame(g = factor(rep(seq_len(40), each = 5)))
cls <- rep(c(1, 2), each = 20)
dg$y <- stats::rnorm(nrow(dg), c(0, 4)[cls[as.integer(dg$g)]], 1)
fg <- frm(bf(y ~ 1) + mixture(gaussian(), gaussian(), groups = ~g),
          data = dg)
cat("group-level mixture:",
    tryCatch({simulate_newdata(fg, dg[1:10, ], nsim = 2); "drew"},
             error = function(e) conditionMessage(e)), "\n")
# a truncated response: the bounds have to follow the rows
set.seed(5)
dt <- data.frame(x = stats::rnorm(200))
dt$y <- stats::rnorm(200, 1 + dt$x, 1)
dt <- dt[dt$y > 0, ]
ft <- frm(bf(y | trunc(lb = 0) ~ x) + gaussian(), data = dt)
set.seed(2)
ynd <- simulate_newdata(ft, dt[1:20, ], nsim = 5)
cat("trunc(lb = 0) at newdata, all draws above the bound:",
    all(ynd > 0), "\n")

cat("\n== the one piece that is not a copy: re_formula = NA ==\n")
# simulate(re_formula = NA) draws FRESH group effects per replicate
# (draw_b()); frm_linpred(re_formula = NA) drops them instead. At
# newdata the two have to be composed: draw b, then evaluate the newdata
# design at that b, which is what ce_draw_new_levels() does for
# conditional_effects (it hands back a fit whose b is the draw).
sim_newdata_re <- function(fit, newdata, nsim = 10) {
  rspec <- frmtmb:::single_response(fit, "simulate()")
  out <- vector("list", nsim)
  for (s in seq_len(nsim)) {
    f2 <- fit
    f2$estimates[["b"]] <- frmtmb:::draw_b(fit)
    dp <- list()
    for (dnm in names(rspec$dpars)) {
      dp[[dnm]] <- as.vector(frm_linpred(f2, newdata = newdata, dpar = dnm,
                                         type = "response",
                                         re_formula = NULL))
    }
    out[[s]] <- sim_draw(sim_context(f2, rspec, dp, aterms = list(),
                                     n = length(dp[[1L]]),
                                     extra = fit_extras(f2)))
  }
  do.call(cbind, out)
}
# 2000 replicates, not 200: at 200 the correlation below has a standard
# error of 0.07 and cannot separate 0 from the 0.116 the shared draw
# implies, which is how this probe first read the wrong answer
set.seed(11)
a <- simulate_newdata(fit, nd, nsim = 2000, re_formula = NA)
set.seed(11)
b <- sim_newdata_re(fit, nd, nsim = 2000)
set.seed(11)
onrows <- as.matrix(simulate(fit, nsim = 2000, re_formula = NA))
cat(sprintf("sd of one row's draws: population %.3f, fresh b %.3f, fitted rows %.3f\n",
            stats::sd(a[1, ]), stats::sd(b[1, ]), stats::sd(onrows[1, ])))
cat(sprintf("the model: sigma %.3f, sd(g) %.3f\n", sigma(fit),
            sqrt(varcorr_matrices(fit)[[1]][1, 1])))

# The spread of ONE row cannot separate the two (200 draws put 1.068 and
# 1.114 inside each other's error). What does is the correlation between
# two rows of the SAME group across replicates: a shared group draw puts
# it at sd_g^2 / (sd_g^2 + sigma^2), and a population draw at 0.
sg <- varcorr_matrices(fit)[[1]][1, 1]
pairs <- do.call(rbind, lapply(levels(droplevels(nd$g)), function(lv) {
  ii <- which(nd$g == lv)
  if (length(ii) < 2L) NULL else t(utils::combn(ii, 2))
}))
cor_of <- function(m) {
  mean(vapply(seq_len(nrow(pairs)), function(k) {
    stats::cor(m[pairs[k, 1], ], m[pairs[k, 2], ])
  }, 0))
}
cat(sprintf("mean within-group correlation over %d row pairs: population %.3f, fresh b %.3f, sd_g^2/(sd_g^2+sigma^2) = %.3f\n",
            nrow(pairs), cor_of(a), cor_of(b), sg / (sg + sigma(fit)^2)))

cat("\n== what calls simulate() today ==\n")
cat("  R/bootstrap.R frm_bootstrap(), R/conditional-effects.R pp_check(),",
    "R/diagnostics.R dharma_residuals()\n")
