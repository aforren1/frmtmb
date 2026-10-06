# Reviewer of lane ordmix: plain ordinal fits, base against lane, to the
# bit, over every family x threshold structure x disc case, the links,
# thres(gr = ), cs(), weights and a random effect. Per model: opt$par,
# logLik, obj$fn and obj$gr at a fixed perturbed point (the tape itself),
# fitted() in sample and on newdata, simulate(seed = 1), fixef(),
# variables(), default_prior() and print(). Data seed 20261096.
# Usage: Rscript dev/ordmix-rev-bitwise.R <base|lane|compare>
args <- commandArgs(TRUE)
arm <- args[1]
root <- "C:/Users/adf44/source/r/frmtmb-wt-ordmix/dev"
if (arm == "compare") {
  a <- readRDS(file.path(root, "ordmix-rev-bitwise-base.rds"))
  b <- readRDS(file.path(root, "ordmix-rev-bitwise-lane.rds"))
  n <- 0L
  same <- 0L
  stopifnot(identical(names(a), names(b)))
  for (m in names(a)) {
    for (o in union(names(a[[m]]), names(b[[m]]))) {
      n <- n + 1L
      ok <- identical(a[[m]][[o]], b[[m]][[o]])
      same <- same + ok
      if (!ok) cat("DIFFERS", m, o, "\n")
    }
  }
  cat(sprintf("BITWISE %d of %d outputs identical over %d models\n", same,
              n, length(a)))
  quit(save = "no")
}
libs <- c("C:/Users/adf44/source/r/rellib-r5",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "lane") libs <- c("C:/Users/adf44/source/r/wt-ordmix-lib", libs)
.libPaths(libs)
suppressPackageStartupMessages(library(frmtmb))
cat("arm", arm, find.package("frmtmb"), "\n")
set.seed(20261096)
n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n),
                g = factor(sample(c("a", "b"), n, TRUE)),
                h = factor(sample(1:10, n, TRUE)),
                w = runif(n, 0.5, 2))
u <- rlogis(n) / exp(0.3 * d$z) + 0.8 * d$x
d$y <- 1L + (u > -1) + (u > 0.3) + (u > 1.5)
d$yh <- ifelse(runif(n) < 0.2, 0L, d$y)
nd <- d[c(1, 5, 9), ]
mods <- list()
for (fam in c("cumulative", "sratio", "cratio", "acat")) {
  for (th in c("flexible", "equidistant", "sum_to_zero")) {
    for (dc in c("held", "modeled")) {
      nm <- paste(fam, th, dc, sep = "_")
      mods[[nm]] <- local({
        fam_ <- fam; th_ <- th; dc_ <- dc
        function() {
          f <- if (dc_ == "held") bf(y ~ x) else bf(y ~ x, disc ~ 0 + z)
          frm(f, family = get(fam_)(threshold = th_), data = d)
        }
      })
    }
  }
}
for (lk in c("probit", "cloglog", "cauchit")) {
  for (fam in c("cumulative", "sratio", "acat")) {
    nm <- paste(fam, lk, sep = "_")
    mods[[nm]] <- local({
      fam_ <- fam; lk_ <- lk
      function() frm(bf(y ~ x), family = get(fam_)(lk_), data = d)
    })
  }
}
mods$cum_gr <- function() frm(bf(y | thres(gr = g) ~ x),
                              family = cumulative(), data = d)
mods$sratio_gr_equi <- function() {
  frm(bf(y | thres(gr = g) ~ x), family = sratio(threshold = "equidistant"),
      data = d)
}
mods$cratio_cs <- function() frm(bf(y ~ cs(x) + z), family = cratio(),
                                 data = d)
mods$acat_cs_disc <- function() frm(bf(y ~ cs(x), disc ~ 0 + z),
                                    family = acat(), data = d)
mods$cum_weights <- function() frm(bf(y | weights(w) ~ x),
                                   family = cumulative(), data = d)
mods$cum_re <- function() frm(bf(y ~ x + (1 | h)), family = cumulative(),
                              data = d)
mods$cum_thres5 <- function() frm(bf(y | thres(4) ~ x),
                                  family = cumulative(), data = d)
for (th in c("flexible", "equidistant", "sum_to_zero")) {
  mods[[paste0("hurdle_", th)]] <- local({
    th_ <- th
    function() frm(bf(yh ~ x, hu ~ z),
                   family = hurdle_cumulative(threshold = th_), data = d)
  })
}
mods$hurdle_disc_probit <- function() {
  frm(bf(yh ~ x, disc ~ 0 + z), family = hurdle_cumulative("probit"),
      data = d)
}
out <- list()
for (m in names(mods)) {
  fit <- suppressWarnings(mods[[m]]())
  set.seed(5)
  q <- fit$opt$par + rnorm(length(fit$opt$par), 0, 0.2)
  out[[m]] <- list(
    par = fit$opt$par, logLik = as.numeric(logLik(fit)),
    fn = fit$obj$fn(q), gr = fit$obj$gr(q),
    fitted = suppressWarnings(fitted(fit)),
    fitted_nd = suppressWarnings(fitted(fit, newdata = nd)),
    sim = simulate(fit, nsim = 2, seed = 1),
    fixef = suppressWarnings(fixef(fit)), variables = variables(fit),
    prior = tryCatch(default_prior(fit), error = function(e) "error"),
    print = capture.output(suppressWarnings(print(fit))))
}
cat("models", length(out), "\n")
saveRDS(out, file.path(root, paste0("ordmix-rev-bitwise-", arm, ".rds")))
