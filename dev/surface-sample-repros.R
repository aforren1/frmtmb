# Lane surface, items 5 and 9: the sampling-route defects.
#
#   Rscript dev/surface-sample-repros.R lane|base > dev/surface-out/sample-<arm>.txt
#
# 5: frm_sample() repeats a warning the frame build gives (brms's cs()
#    experimental warning on cumulative(), 6 times at 0.68.0), on both
#    routes.
# 9: frm_sample(fit) on an exact y ~ gp(x) fit does not move (the gpby
#    review's construction: 60 points, data seed 5, chains = 1,
#    iter = 600, seed = 4).
arm <- commandArgs(TRUE)[1]
source("dev/surface-env.R")
surface_env(arm)
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("arm", arm, "| frmtmb", as.character(packageVersion("frmtmb")), "from",
    find.package("frmtmb"), "| frmtmb.sample",
    as.character(packageVersion("frmtmb.sample")), "from",
    find.package("frmtmb.sample"), "\n")
options(mc.cores = 1)
count_warn <- function(tag, expr) {
  w <- character()
  r <- tryCatch(withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage")),
  error = function(e) structure(conditionMessage(e), class = "err"))
  tab <- table(substr(w, 1, 60))
  cat(sprintf("%-40s %s; %d warning(s)\n", tag,
              if (inherits(r, "err")) paste("ERROR:", substr(r, 1, 120))
              else "ok", length(w)))
  for (k in names(tab)) cat(sprintf("    %2d x %s\n", tab[[k]], k))
  invisible(r)
}

## 5
set.seed(20261006)
n <- 400
d <- data.frame(x = rnorm(n), z = rnorm(n))
u <- rlogis(n) + 0.5 * d$x + 0.4 * d$z
d$y <- 1L + (u > -1 + 0.15 * d$x) + (u > 0.3) + (u > 1.5 - 0.15 * d$x)
fit5 <- count_warn("5 frm(y ~ z + cs(x), cumulative)",
                   frm(y ~ z + cs(x), family = cumulative(), data = d))
count_warn("5 frm_sample(formula route)",
           frm_sample(y ~ z + cs(x), family = cumulative(), data = d,
                      chains = 1, iter = 300, refresh = 0, seed = 3,
                      prior = set_prior("normal(0, 2)", class = "b")))
count_warn("5 frm_sample(fit route)",
           frm_sample(fit5, chains = 1, iter = 300, refresh = 0, seed = 3))
count_warn("5 brms::make_standata (brms's own count)",
           brms::make_standata(y ~ z + cs(x), family = brms::cumulative(),
                               data = d))

## 9
set.seed(5)
n <- 60
d9 <- data.frame(x = round(stats::runif(n, 0, 6), 1))
d9$y <- 0.5 + sin(d9$x) + stats::rnorm(n, 0, 0.3)
fit9 <- frm(bf(y ~ gp(x)), family = gaussian(), data = d9)
for (route in c("fit", "formula")) {
  ds <- count_warn(paste("9 frm_sample,", route, "route"),
    if (route == "fit") {
      frm_sample(fit9, chains = 1, iter = 600, refresh = 0, seed = 4)
    } else {
      frm_sample(bf(y ~ gp(x)), data = d9, family = gaussian(), chains = 1,
                 iter = 600, refresh = 0, seed = 4)
    })
  if (inherits(ds, "err")) next
  m <- as.matrix(ds)
  s <- apply(m, 2, stats::sd)
  sp <- rstan::get_sampler_params(ds$stanfit, inc_warmup = FALSE)[[1]]
  cat(sprintf(paste0("    %d columns, %d with sd 0; accept %.3f, ",
                     "stepsize %.3g, divergent %d\n"),
              ncol(m), sum(s == 0), mean(sp[, "accept_stat__"]),
              sp[1, "stepsize__"], sum(sp[, "divergent__"])))
}
cat("SAMPLE REPROS DONE\n")
