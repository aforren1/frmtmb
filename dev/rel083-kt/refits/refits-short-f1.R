# Short check of the refits fix round, F1: refit_bernoulli_codes() and
# boot_refits(). Seed 9, n = 40 (probe5 data) unless stated.
source("C:/Users/adf44/source/r/frmtmb-wt-rel083/dev/rel083-kt/refits/refits-review-lib.R")
cat("build has fix:",
    any(grepl("reads two ways", deparse(frmtmb:::refit_bernoulli_codes))),
    any(grepl("as.logical", deparse(frmtmb:::boot_refits))), "\n")
set.seed(9)
n <- 40
x <- rnorm(n)
s <- rbinom(n, 1, plogis(0.3 + x))
mk <- function(y) frm(y ~ x, data = data.frame(x = x, y = y),
                      family = bernoulli())
cf <- function(f) fixef(f, flatten = TRUE)
same <- function(a, b) identical(cf(a), cf(b))
q <- function(e) suppressWarnings(e)
f01 <- mk(s)
f12 <- mk(s + 1)
f02 <- mk(2 * s)
ff <- mk(factor(s + 1, levels = c("1", "2")))
fl <- mk(s == 1)
for (nm in c("f01", "f12", "f02", "ff", "fl")) {
  cat(nm, "bin_levels:", get(nm)$spec$responses[[1]]$family$bin_levels, "\n")
}
r01_1 <- q(refit(f01, rep(1, n)))
r01_0 <- q(refit(f01, rep(0, n)))

cat("\n== A. fit {0,1} ==\n")
go("A1 refit(f01, rep(1, n)) not refused, intercept sign",
   sign(cf(r01_1)[1]))
go("A2 refit(f01, rep(0, n)) intercept sign", sign(cf(r01_0)[1]))
go("A3 refit(f01, s) = f01", same(refit(f01, s), f01))
go("A4 refit(f01, as.character(s)) = f01", same(refit(f01, as.character(s)), f01))
go("A5 refit(f01, s == 1) = f01", same(refit(f01, s == 1), f01))
go("A6 refit(f01, c(s[-1], 2)) refused", refit(f01, c(s[-1], 2)))

cat("\n== B. fit {0,2} ==\n")
go("B1 refit(f02, s) (0/1, 1 not a value: codes) = refit(f01, s)",
   same(refit(f02, s), refit(f01, s)))
go("B2 refit(f02, 2 * s) (values) = refit(f01, s)",
   same(refit(f02, 2 * s), refit(f01, s)))
go("B3 refit(f02, rep(1, n)) = refit(f01, rep(1, n)) (1 a code: success)",
   same(q(refit(f02, rep(1, n))), r01_1))
go("B4 refit(f02, rep(2, n)) = refit(f01, rep(1, n))",
   same(q(refit(f02, rep(2, n))), r01_1))
go("B5 refit(f02, rep(0, n)) = refit(f01, rep(0, n))",
   same(q(refit(f02, rep(0, n))), r01_0))
go("B6 refit(f02, c(0,1,2,...)) mix refused",
   refit(f02, c(0, 1, 2, s[-(1:3)] * 2)))
go("B7 residuals(f02, newdata y = 1) (residuals reads values only)",
   residuals(f02, newdata = data.frame(x = x, y = 1)))

cat("\n== C. fit factor levels c('1','2') ==\n")
go("C1 refit(ff, s + 1) numeric values = refit(f01, s)",
   same(refit(ff, s + 1), refit(f01, s)))
go("C2 refit(ff, s) numeric codes = refit(f01, s)",
   same(refit(ff, s), refit(f01, s)))
go("C3 refit(ff, rep(1, n)) refused (ambiguous)", refit(ff, rep(1, n)))
go("C4 refit(ff, factor(rep('1', n))) = refit(f01, rep(0, n))",
   same(q(refit(ff, factor(rep("1", n)))), r01_0))
go("C5 refit(ff, rep(TRUE, n)) = refit(f01, rep(1, n))",
   same(q(refit(ff, rep(TRUE, n))), r01_1))
go("C6 refit(ff, factor(s, levels = 0:1)) refused (factor states values)",
   refit(ff, factor(s, levels = 0:1)))
go("C7 refit(ff, factor(s + 1, levels = c('2','1'))) = refit(f01, s)",
   same(refit(ff, factor(s + 1, levels = c("2", "1"))), refit(f01, s)))
go("C8 refit(f12, factor(rep('1', n))) = refit(f01, rep(0, n))",
   same(q(refit(f12, factor(rep("1", n)))), r01_0))
go("C9 refit(f12, rep(TRUE, n)) = refit(f01, rep(1, n))",
   same(q(refit(f12, rep(TRUE, n))), r01_1))
go("C10 refit(f12, rep(2, n)) = refit(f01, rep(1, n))",
   same(q(refit(f12, rep(2, n))), r01_1))
go("C11 refit(f12, rep(0, n)) = refit(f01, rep(0, n)) (0 only a code)",
   same(q(refit(f12, rep(0, n))), r01_0))
go("C12 refit(fl, rep(1, n)) logical fit: 1 = TRUE = success",
   same(q(refit(fl, rep(1, n))), r01_1))
go("C13 refit(fl, rep(FALSE, n))", same(q(refit(fl, rep(FALSE, n))), r01_0))

cat("\n== D. NA ==\n")
go("D1 refit(f12, c(NA, s[-1] + 1))", refit(f12, c(NA, s[-1] + 1)))
go("D2 refit(f01, c(NA, s[-1]))", refit(f01, c(NA, s[-1])))
go("D3 refit(ff, factor(c(NA, s[-1] + 1)))", refit(ff, factor(c(NA, s[-1] + 1))))
go("D4 refit(f12, c(NA, rep(1, n - 1)))", refit(f12, c(NA, rep(1, n - 1))))

cat("\n== E. -1/0 fit ==\n")
fm <- mk(s - 1)
go("E1 refit(fm, rep(0, n)) refused (0 = value success, code failure)",
   refit(fm, rep(0, n)))
go("E2 refit(fm, rep(-1, n)) = refit(f01, rep(0, n))",
   same(q(refit(fm, rep(-1, n))), r01_0))
go("E3 refit(fm, s - 1) = refit(f01, s)", same(refit(fm, s - 1), refit(f01, s)))

cat("\n== F. bootstrap, all-one-class draws ==\n")
set.seed(21)
n2 <- 12
x2 <- rnorm(n2)
mk2 <- function(y) frm(y ~ x, data = data.frame(x = x2, y = y),
                       family = bernoulli())
for (dir in c("succ", "fail")) {
  s2 <- if (dir == "succ") c(rep(1, 11), 0) else c(rep(0, 11), 1)
  codings <- list(c01 = s2, c12 = s2 + 1, c02 = 2 * s2, cm10 = s2 - 1,
                  fac = factor(ifelse(s2 == 1, "yes", "no")),
                  lgl = s2 == 1)
  fits <- lapply(codings, mk2)
  sims <- simulate(fits$c01, nsim = 40, seed = 5)
  cat(dir, ": draws all 1:", sum(vapply(sims, function(v) all(v == 1), NA)),
      " all 0:", sum(vapply(sims, function(v) all(v == 0), NA)), "of 40\n")
  bs <- lapply(fits, function(f) q(frm_bootstrap(f, nsim = 40, seed = 5)))
  for (k in names(bs)[-1]) {
    cat(sprintf("  %-5s vs c01: identical t %s, failed %d vs %d\n", k,
                identical(bs[[k]]$t, bs$c01$t), sum(!bs[[k]]$converged),
                sum(!bs$c01$converged)))
  }
  cat("  NA rows in c01 t:", sum(!complete.cases(bs$c01$t)), "\n")
}

cat("\n== G. update(newdata =) and influence() ==\n")
nd <- function(y) data.frame(x = x, y = y)
go("G1 update(f12, newdata y = s + 1) = f01",
   same(q(update(f12, newdata = nd(s + 1))), f01))
go("G2 update(f12, newdata y = s) (0/1) = f01",
   same(q(update(f12, newdata = nd(s))), f01))
# constant responses separate: the optimizer stops at different points
# from a warm and a cold start, so compare the sign of the intercept
sg <- function(f) unname(sign(cf(f)[1]))
go("G3 update(f12, newdata y = 1s): intercept sign (+1 success, -1 failure)",
   sg(q(update(f12, newdata = nd(rep(1, n))))))
go("G4 update(f12, newdata y = 2s): intercept sign",
   sg(q(update(f12, newdata = nd(rep(2, n))))))
go("G5 update(f12, newdata y = 0s): intercept sign",
   sg(q(update(f12, newdata = nd(rep(0, n))))))
go("G6 update(ff, newdata y = factor('1') levels 1,2): intercept sign",
   sg(q(update(ff, newdata = nd(factor(rep("1", n), levels = c("1", "2")))))))
go("G6b brms 2.23.0 fresh standata(y = 1s)$Y and update(newdata) rebuilds basis", {
  suppressPackageStartupMessages(library(brms))
  unique(standata(y ~ x, data = nd(rep(1, n)), family = bernoulli())$Y) })
# influence: one deletion leaves only failures (the value 1 on 1/2)
set.seed(31)
n3 <- 30
x3 <- rnorm(n3)
s3 <- rep(0, n3)
s3[which.min(abs(x3))] <- 1
dd3 <- data.frame(x = x3, s = s3)
i12 <- q(influence(frm(y ~ x, data = transform(dd3, y = s + 1),
                       family = bernoulli()), force = TRUE))
i01 <- q(influence(frm(y ~ x, data = transform(dd3, y = s),
                       family = bernoulli()), force = TRUE))
k <- which(s3 == 1)
go("G7 influence 1/2 vs 0/1: identical fixed tables",
   identical(unname(as.matrix(i01$fixed)), unname(as.matrix(i12$fixed))))
go("G7 deletion of the only success: row k of each",
   list(k = k, f01 = as.matrix(i01$fixed)[k, ], f12 = as.matrix(i12$fixed)[k, ]))
go("G7 max|diff| other rows",
   max(abs(as.matrix(i01$fixed)[-k, ] - as.matrix(i12$fixed)[-k, ])))
