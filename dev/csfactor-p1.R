# Punch round 1 fixes, measured. Seeds are the reviewer's.
#   Rscript dev/csfactor-p1.R <lib> > dev/csfactor-log/p1.txt
args <- commandArgs(TRUE)
LIB <- if (length(args)) args[[1L]] else
  "C:/Users/adf44/source/r/wt-csfactor-lib"
.libPaths(c(LIB, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb from", dirname(system.file(package = "frmtmb")), "\n")
go <- function(lab, expr) {
  cat("\n--", lab, "--\n")
  print(tryCatch(expr, error = function(e) paste("ERROR:",
                                                 conditionMessage(e))))
}

## B2: the mo() hole. Reviewer's dev/csfactor-rev-gap.R, seed 1907.
set.seed(1907)
n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n))
d$m <- factor(sample(1:4, n, TRUE), ordered = TRUE)
d$f <- factor(sample(c("a", "b", "c"), n, TRUE))
d$yo <- sample(1:3, n, TRUE)
go("mo(m) + cs(m) must be REFUSED",
   frm(bf(yo ~ mo(m) + cs(m)), family = sratio(), data = d))
go("cs(m) alone must FIT",
   as.numeric(logLik(frm(bf(yo ~ cs(m)), family = sratio(), data = d))))
go("mo(m) alone must FIT",
   as.numeric(logLik(frm(bf(yo ~ mo(m)), family = sratio(), data = d))))
go("mo(m) + cs(x) must FIT (no false alarm)",
   as.numeric(logLik(frm(bf(yo ~ mo(m) + cs(x)), family = sratio(),
                         data = d))))
go("mo(m) + cs(f) must FIT: f is not a coarsening of m here",
   as.numeric(logLik(frm(bf(yo ~ mo(m) + cs(f)), family = sratio(),
                         data = d))))
# a COARSENING of m: mo(m) is generically not inside span(1, fc2 dummies)
d$mc <- factor(ifelse(as.integer(d$m) <= 2L, "lo", "hi"))
go("mo(m) + cs(mc), mc a 2-level coarsening of m, must FIT",
   as.numeric(logLik(frm(bf(yo ~ mo(m) + cs(mc)), family = sratio(),
                         data = d))))
go("mo(m):z + cs(m) must FIT: z times a function of m is not in cs(m) span",
   frm(bf(yo ~ mo(m):z + cs(m)), family = sratio(), data = d))
go("me(x, sdx) + cs(x): accepted, standard errors finite?", {
  d$sdx <- 0.3
  fm <- frm(bf(yo ~ me(x, sdx) + cs(x)), family = sratio(), data = d)
  c(logLik = as.numeric(logLik(fm)),
    nan_se = sum(is.na(fixef(fm)[, "Est.Error"])),
    npar = nrow(fixef(fm)))
})

## Nit 5: the advice must match what the formula holds.
go("poly(x, 2) + cs(x) message",
   frm(bf(yo ~ poly(x, 2) + cs(x)), family = sratio(), data = d))
go("s(x) + cs(x) message",
   frm(bf(yo ~ s(x) + cs(x)), family = sratio(), data = d))
go("x + cs(x) message (unchanged advice)",
   frm(bf(yo ~ x + cs(x)), family = sratio(), data = d))
go("x + z + cs(I(x + z)) message names both columns",
   frm(bf(yo ~ x + z + cs(I(x + z))), family = sratio(), data = d))

## B1: multivariate cs() fits. Reviewer's dev/csfactor-rev-docs.R, seed 77.
set.seed(77)
n2 <- 300
dm <- data.frame(x = rnorm(n2), z = rnorm(n2))
dm$fc <- factor(sample(c("a", "b", "c"), n2, TRUE))
dm$yo <- sample(1:3, n2, TRUE)
dm$yo2 <- sample(1:3, n2, TRUE)
go("multivariate cs() in both responses", {
  mv <- frm(bf(yo ~ x + cs(fc)) + bf(yo2 ~ z + cs(fc)), family = sratio(),
            data = dm)
  u1 <- frm(bf(yo ~ x + cs(fc)), family = sratio(), data = dm)
  u2 <- frm(bf(yo2 ~ z + cs(fc)), family = sratio(), data = dm)
  list(npar = nrow(fixef(mv)), rows = rownames(fixef(mv)),
       vars = grep("^bcs", variables(mv), value = TRUE),
       logLik = as.numeric(logLik(mv)),
       sum_univariate = as.numeric(logLik(u1)) + as.numeric(logLik(u2)),
       difference = as.numeric(logLik(mv)) -
         (as.numeric(logLik(u1)) + as.numeric(logLik(u2))))
})
go("multivariate with ONE bad predictor is still refused",
   frm(bf(yo ~ x + cs(x)) + bf(yo2 ~ z + cs(fc)), family = sratio(),
       data = dm))

## B3: unused level, and the influence() warning. Reviewer's seeds.
set.seed(405)
n3 <- 240
du <- data.frame(x = rnorm(n3))
du$f <- factor(c(rep("a", 114), rep("b", 126)), levels = c("a", "b", "c"))
du$yo <- sample(1:3, n3, TRUE)
go("an unused level gives ONE cs column, no refusal", {
  fu <- frm(bf(yo ~ x + cs(f)), family = sratio(), data = du)
  list(rows = rownames(fixef(fu)),
       cs_cols = length(fu$frame$linpreds[["yo.mu"]][["cs"]]),
       logLik = as.numeric(logLik(fu)))
})

set.seed(90291)
n4 <- 120
di <- data.frame(x = rnorm(n4))
di$f <- factor(c(rep("a", 60), rep("b", 59), "c"),
               levels = c("a", "b", "c"))
di$yo <- sample(1:3, n4, TRUE)
go("influence(): the dropped-coefficient warning fires, with the unit", {
  fi <- frm(bf(yo ~ x + cs(f)), family = sratio(), data = di)
  w <- NULL
  inf <- withCallingHandlers(
    influence(fi, force = TRUE),
    warning = function(cnd) {
      w <<- c(w, conditionMessage(cnd))
      invokeRestart("muffleWarning")
    })
  cd <- cooks.distance(inf)
  list(warnings = w, n_na_cooks = sum(is.na(cd)),
       which_na = which(is.na(cd)))
})
cat("\ndone\n")
