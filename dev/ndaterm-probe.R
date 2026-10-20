# Lane ndaterm: does every newdata path read the addition terms from
# newdata? Core paths only (no sampler). Seed 7 throughout.
# Usage: Rscript dev/ndaterm-probe.R <arm>, arm in mine, r9, r8, r7, r6.
# Logs: dev/ndaterm-out/probe-<arm>.txt
args <- commandArgs(TRUE)
arm <- if (length(args)) args[1] else "r9"
base <- "C:/Users/adf44/source/r"
libs <- switch(arm,
  mine = c(file.path(base, "wt-ndaterm-lib"), file.path(base, "rellib-r9")),
  file.path(base, paste0("rellib-", arm)))
.libPaths(c(libs, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.eam)})
cat("arm", arm, "frmtmb", format(packageVersion("frmtmb")), "from",
    find.package("frmtmb"), " eam", format(packageVersion("frmtmb.eam")),
    "\n")
try_it <- function(label, expr) {
  r <- tryCatch(expr, error = function(e) {
    structure(paste("ERROR:", conditionMessage(e)), class = "err")
  })
  if (inherits(r, "err")) cat(label, ":", r, "\n")
  invisible(r)
}
ok <- function(r) !inherits(r, "err")
show <- function(label, x) {
  cat(sprintf("%-46s %s\n", label,
              paste(format(x, digits = 5), collapse = "  ")))
}

cat("\n== binomial trials(N), newdata N = 1000, 2\n")
set.seed(7)
n <- 200
db <- data.frame(x = rnorm(n), N = sample(5:15, n, TRUE))
db$y <- rbinom(n, db$N, plogis(0.3 + 0.5 * db$x))
fb <- frm(bf(y | trials(N) ~ x), family = binomial(), data = db)
ndb <- data.frame(x = c(0, 0), N = c(1000, 2))
set.seed(7)
p <- try_it("predict(summary = FALSE)",
            predict(fb, newdata = ndb, summary = FALSE, ndraws = 400))
if (ok(p)) {
  show("predict(summary = FALSE) col means", colMeans(p))
  show("predict(summary = FALSE) col max", apply(p, 2, max))
}
ps <- try_it("predict() summary", predict(fb, newdata = ndb, ndraws = 400))
if (ok(ps)) show("predict() Estimate", ps[, "Estimate"])
s <- try_it("simulate", simulate(fb, nsim = 400, newdata = ndb, seed = 7))
if (ok(s)) show("simulate row means", rowMeans(as.matrix(s)))
f <- try_it("fitted", fitted(fb, newdata = ndb))
if (ok(f)) show("fitted Estimate", f[, "Estimate"])
pr <- try_it("predict(newdata without N)",
             predict(fb, newdata = data.frame(x = 0), ndraws = 20))
if (ok(pr)) show("predict(no N) Estimate", pr[, "Estimate"])

cat("\n== zero_inflated_binomial and beta_binomial trials(N)\n")
set.seed(7)
dz <- data.frame(x = rnorm(n), N = sample(5:15, n, TRUE))
dz$y <- rbinom(n, dz$N, plogis(0.2 + 0.4 * dz$x)) *
  rbinom(n, 1, 0.8)
for (fm in list(zero_inflated_binomial(), beta_binomial())) {
  fz <- try_it(fm$family, frm(bf(y | trials(N) ~ x), family = fm,
                              data = dz))
  if (!ok(fz)) next
  set.seed(7)
  pz <- try_it("predict", predict(fz, newdata = ndb, summary = FALSE,
                                  ndraws = 200))
  if (ok(pz)) show(paste(fm$family, "predict col means"), colMeans(pz))
  fz2 <- fitted(fz, newdata = ndb)
  show(paste(fm$family, "fitted"), fz2[, "Estimate"])
}

cat("\n== gaussian se(s), newdata s = 0.05, 5\n")
set.seed(7)
dg <- data.frame(x = rnorm(n), s = runif(n, 0.2, 1))
dg$y <- 1 + dg$x + rnorm(n, 0, dg$s)
ndg <- data.frame(x = c(0, 0), s = c(0.05, 5))
for (sg in c(FALSE, TRUE)) {
  fg <- if (sg) {
    frm(bf(y | se(s, sigma = TRUE) ~ x), family = gaussian(), data = dg)
  } else frm(bf(y | se(s) ~ x), family = gaussian(), data = dg)
  lab <- if (sg) "se(s, sigma = TRUE)" else "se(s)"
  set.seed(7)
  pg <- try_it(lab, predict(fg, newdata = ndg, summary = FALSE,
                            ndraws = 400))
  if (ok(pg)) show(paste(lab, "predict col sd"), apply(pg, 2, sd, na.rm = TRUE))
  sgd <- simulate(fg, nsim = 400, newdata = ndg, seed = 7)
  show(paste(lab, "simulate row sd"), apply(as.matrix(sgd), 1, sd))
}

cat("\n== gaussian ARMA cov = FALSE with se(): fit-time refusal\n")
set.seed(7)
da <- data.frame(t = 1:120, s = runif(120, 0.2, 1))
da$y <- as.numeric(arima.sim(list(ar = 0.5), 120)) + rnorm(120, 0, da$s)
try_it("frm(y | se(s) ~ 1 + arma(cov = FALSE))",
       frm(bf(y | se(s) ~ 1 + ar(time = t, p = 1, cov = FALSE)),
           family = gaussian(), data = da))
cat("\n== gaussian ARMA cov = FALSE with trunc(): newdata fill\n")
da$y <- pmax(da$y, -2.5)
fa <- try_it("fit", frm(bf(y | trunc(lb = -3) ~ 1 + ar(time = t, p = 1,
                                                    cov = FALSE)),
                        family = gaussian(), data = da))
if (ok(fa)) {
  nda <- da[100:110, ]
  nda$y[3:6] <- NA
  set.seed(7)
  pa <- try_it("predict", predict(fa, newdata = nda, summary = FALSE,
                                  ndraws = 200))
  if (ok(pa)) show("ARMA trunc predict min over draws", min(pa))
}

cat("\n== custom family with vint() trials and vreal() known sd\n")
cbin <- custom_family(
  "cbin", dpars = "mu", links = list(mu = "logit"),
  lpdf = function(y, dpars, aterms) {
    RTMB::dbinom(y, aterms[["vint1"]], dpars[["mu"]], log = TRUE)
  },
  required_aterms = "vint1", type = "discrete",
  post = list(mean_fn = function(dpars, aterms) {
    dpars[["mu"]] * aterms[["vint1"]]
  }),
  sim = function(dpars, aterms, n) {
    stats::rbinom(n, aterms[["vint1"]], dpars[["mu"]])
  })
fc <- frm(bf(y | vint(N) ~ x), family = cbin, data = db)
set.seed(7)
pc <- try_it("vint predict", predict(fc, newdata = ndb, summary = FALSE,
                                     ndraws = 400))
if (ok(pc)) show("vint predict col means", colMeans(pc))
show("vint simulate row means",
     rowMeans(as.matrix(simulate(fc, nsim = 400, newdata = ndb, seed = 7))))
show("vint fitted", fitted(fc, newdata = ndb)[, "Estimate"])
cnorm <- custom_family(
  "cnorm", dpars = "mu", links = list(mu = "identity"),
  lpdf = function(y, dpars, aterms) {
    RTMB::dnorm(y, dpars[["mu"]], aterms[["vreal1"]], log = TRUE)
  },
  required_aterms = "vreal1",
  post = list(mean_fn = function(dpars, aterms) dpars[["mu"]]),
  sim = function(dpars, aterms, n) {
    stats::rnorm(n, dpars[["mu"]], aterms[["vreal1"]])
  })
fr <- frm(bf(y | vreal(s) ~ x), family = cnorm, data = dg)
set.seed(7)
pv <- try_it("vreal predict", predict(fr, newdata = ndg, summary = FALSE,
                                      ndraws = 400))
if (ok(pv)) show("vreal predict col sd", apply(pv, 2, sd))
show("vreal simulate row sd",
     apply(as.matrix(simulate(fr, nsim = 400, newdata = ndg, seed = 7)), 1,
           sd))

cat("\n== multinomial trials(N)\n")
set.seed(7)
dm <- data.frame(x = rnorm(n), N = sample(5:15, n, TRUE))
dm$Y <- t(vapply(seq_len(n), function(i) {
  rmultinom(1, dm$N[i], c(0.5, 0.3, 0.2))
}, numeric(3)))
fm3 <- try_it("multinomial fit",
              frm(bf(Y | trials(N) ~ x), family = multinomial(K = 3), data = dm))
if (ok(fm3)) {
  set.seed(7)
  pm <- try_it("multinomial predict", predict(fm3, newdata = ndb,
                                               summary = FALSE, ndraws = 50))
  if (ok(pm)) show("multinomial predict row sums, draw 1",
                   rowSums(pm[1, , ]))
}

cat("\n== wiener dec(upper)\n")
set.seed(7)
dat <- ddm_simulate(300, mu = 0.8, bs = 1.4, ndt = 0.3, bias = 0.3)
dat$x <- rnorm(300)
fw <- frm(bf(rt | dec(upper) ~ x, bias = 0.3), family = wiener(),
          data = dat)
ndw <- data.frame(upper = c(1, 0, 1), x = 0)
set.seed(7)
pw <- try_it("wiener predict(summary = FALSE)",
             predict(fw, newdata = ndw, summary = FALSE, ndraws = 400))
if (ok(pw)) show("wiener predict col means", colMeans(pw))
sw <- try_it("wiener simulate", simulate(fw, nsim = 400, newdata = ndw,
                                         seed = 7))
if (ok(sw)) show("wiener simulate row means", rowMeans(as.matrix(sw)))
fwf <- try_it("wiener fitted", fitted(fw, newdata = ndw))
if (ok(fwf)) show("wiener fitted Estimate", fwf[, "Estimate"])

cat("\n== conditional_effects(method = 'predict')\n")
set.seed(7)
ce <- try_it("binomial ce predict",
             conditional_effects(fb, "x", method = "predict",
                                 conditions = data.frame(N = 1000)))
if (ok(ce)) show("binomial ce estimate range", range(ce[[1]]$estimate__))
ce2 <- try_it("wiener ce predict",
              conditional_effects(fw, "x", method = "predict"))
if (ok(ce2)) show("wiener ce estimate range", range(ce2[[1]]$estimate__))
ce3 <- try_it("wiener ce predict, upper = 0",
              conditional_effects(fw, "x", method = "predict",
                                  conditions = data.frame(upper = 0)))
if (ok(ce3)) show("wiener ce estimate range, upper = 0",
                  range(ce3[[1]]$estimate__))
cat("done\n")
