# Short check of lane priortab's fix round: the final-stop guard, the
# narrowed support box, and 4925 on frm_sample().
# PRIORTAB_LIB unset = lane library; PRIORTAB_LIB=none = rellib-r20.
lib <- Sys.getenv("PRIORTAB_LIB", "C:/Users/adf44/source/r/wt-priortab-lib")
if (identical(lib, "none")) lib <- ""
.libPaths(c(if (nzchar(lib)) lib, "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("frmtmb", dirname(find.package("frmtmb")), "\n")
options(width = 150, digits = 10)
show <- function(label, expr) {
  cat("\n=== ", label, "\n", sep = "")
  r <- tryCatch(withCallingHandlers(force(expr),
         warning = function(w) {
           cat("WARN: ", substr(gsub("\n", " ", conditionMessage(w)), 1, 150),
               "\n")
           invokeRestart("muffleWarning")
         },
         message = function(m) invokeRestart("muffleMessage")),
       error = function(e) {
         cat("ERROR: ", substr(conditionMessage(e), 1, 200), "\n")
         NULL
       })
  if (!is.null(r)) print(r)
  invisible(r)
}
set.seed(1)
n <- 120
d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:10, each = 12)))
d$y <- 1 + d$z * 0.3 + rnorm(10)[d$g] * 0.3 + rnorm(n, sd = 0.3)
d$yneg <- d$y - 0.8 * d$z
s <- function(f) c(conv = f$opt$convergence, obj = f$opt$objective,
                   z = unname(f$opt$par[names(f$opt$par) == "beta"][2]),
                   msg = f$opt$message)
sp <- function(dn, ...) set_prior(dn, class = "b", coef = "z", ...)
show("A1 beta(2,2) no bound, yneg", s(frm(yneg ~ z, data = d, prior = sp("beta(2, 2)"))))
show("A2 beta(2,2) lb = 0, yneg", s(frm(yneg ~ z, data = d, prior = sp("beta(2, 2)", lb = 0))))
show("A3 lognormal lb = 0, yneg + (1 | g)",
     s(frm(yneg ~ z + (1 | g), data = d, prior = sp("lognormal(0, 1)", lb = 0))))
show("A4 gamma(0.5, 1) lb = 0, yneg",
     s(frm(yneg ~ z, data = d, prior = sp("gamma(0.5, 1)", lb = 0))))
show("A5 lognormal lb = 0, yneg, optim",
     s(frm(yneg ~ z, data = d, prior = sp("lognormal(0, 1)", lb = 0),
           control = frmtmb_control(optimizer = "optim"))))
show("A6 uniform(0, 0.2) no bound, y + 0.5 z (valid stop on the bound)",
     s(frm(I(y + 0.5 * z) ~ z, data = d, prior = sp("uniform(0, 0.2)"))))
show("A7 exponential(1) no bound, yneg (boxed at 0, finite)",
     s(frm(yneg ~ z, data = d, prior = sp("exponential(1)"))))
show("A8 beta(2,2) no bound, y (interior, base 59.80059518)",
     s(frm(y ~ z, data = d, prior = sp("beta(2, 2)"))))
show("A9 sd lb = 0.05 box (valid stop on a bound? finite)",
     s(frm(y ~ z + (1 | g), data = d, prior = set_prior("", class = "sd",
                                                        lb = 0.5))))
show("A10 update() of the A2 fit", {
  f <- suppressWarnings(frm(yneg ~ z, data = d, prior = sp("beta(2, 2)", lb = 0)))
  s(update(f, data = d))
})
show("A11 zi uniform(0, 0.3) dpar class box", {
  dd <- d; dd$cnt <- rpois(n, 2) * rbinom(n, 1, 0.5)
  f <- frm(cnt ~ z, data = dd, family = zero_inflated_poisson(),
           prior = set_prior("uniform(0, 0.3)", class = "zi"))
  c(conv = f$opt$convergence, obj = f$opt$objective,
    zi_internal = unname(f$opt$par[grep("zi", names(f$opt$par))]),
    logit0.3 = qlogis(0.3))
})
smp <- function(fit, pr) {
  dr <- frm_sample(fit, prior = pr, chains = 1, iter = 400, warmup = 200,
                   seed = 3, refresh = 0)
  zz <- as.matrix(dr)[, "b_z"]
  c(n = length(zz), min = min(zz), max = max(zz))
}
show("S1 frm_sample beta(2,2) on z, no bound, data slope 0.26",
     smp(frm(y ~ z, data = d), sp("beta(2, 2)")))
show("S2 frm_sample lognormal(0,1) on z, no bound, data slope 0.26",
     smp(frm(y ~ z, data = d), sp("lognormal(0, 1)")))
show("S3 frm_sample beta(2,2) on z, no bound, data slope -0.5",
     smp(frm(yneg ~ z, data = d), sp("beta(2, 2)")))
show("S4 frm_sample uniform(0, 0.2), no bound (4925 as filed)",
     smp(frm(I(y + 0.5 * z) ~ z, data = d), sp("uniform(0, 0.2)")))
cat("\ndone\n")
