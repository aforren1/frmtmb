# Lane surface, item 9, third step: is the defect the exact mode, or
# the exact point the tape was recorded at?
#
#   Rscript dev/surface-gp-diag3.R lane|base > dev/surface-out/gp-diag3-<arm>.txt
arm <- commandArgs(TRUE)[1]
source("dev/surface-env.R")
surface_env(arm)
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("arm", arm, "| frmtmb", as.character(packageVersion("frmtmb")), "\n")
ns <- asNamespace("frmtmb.sample")
set.seed(5)
n <- 60
d9 <- data.frame(x = round(stats::runif(n, 0, 6), 1))
d9$y <- 0.5 + sin(d9$x) + stats::rnorm(n, 0, 0.3)
fit <- frm(bf(y ~ gp(x)), family = gaussian(), data = d9)
rp <- suppressMessages(ns$sample_resolve_priors(fit, NULL, base = fit$prior,
                                                defaults = TRUE))
nll <- frmtmb::build_objective(fit$frame)
nlp <- frmtmb::neg_log_prior_fn(rp$ri$entries)
tape_at <- function(pars) {
  RTMB::MakeADFun(function(p) nll(p) + nlp(p), pars, random = "b",
                  map = fit$frame[["map"]], silent = TRUE)
}
run <- function(tag, obj, x0) {
  sf <- suppressWarnings(tmbstan::tmbstan(
    obj, chains = 1, iter = 300, refresh = 0, seed = 4,
    init = list(array(x0, dim = length(x0)))))
  sp <- rstan::get_sampler_params(sf, inc_warmup = TRUE)[[1]]
  cat(sprintf("%-46s post-warmup accept %.3f, divergent %d\n", tag,
              mean(sp[151:300, "accept_stat__"]),
              sum(sp[151:300, "divergent__"])))
}
est <- fit$estimates
lpb <- unlist(est[names(fit$frame$par_template)])
jit <- est
set.seed(2)
for (k in names(jit)) jit[[k]] <- jit[[k]] + rnorm(length(jit[[k]]), 0, 1e-3)
run("taped at mode, init at mode", tape_at(est), lpb)
run("taped at mode, init at mode + 1e-12", tape_at(est), lpb + 1e-12)
run("taped at mode + N(0, 1e-3), init at mode", tape_at(jit), lpb)
# which coordinate: zero the field at the tape, keep the rest
z <- est
z$b[] <- 0
run("taped with field 0, init at mode", tape_at(z), lpb)
o <- tape_at(est)
cat("f at mode, order 0 twice:", o$env$f(lpb, order = 0),
    o$env$f(lpb, order = 0), "\n")
g1 <- o$env$f(lpb, order = 1)
g2 <- tape_at(jit)$env$f(lpb, order = 1)
cat("gradient at the mode from the two tapes: max abs diff",
    format(max(abs(g1 - g2))), "| max abs g", format(max(abs(g1))), "\n")
