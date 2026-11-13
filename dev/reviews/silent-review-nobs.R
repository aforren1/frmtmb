# Review of item 1562. Usage: Rscript silent-review-nobs.R lane|base
arm <- commandArgs(TRUE)[1]
libs <- c("C:/Users/adf44/source/r/rellib-r20",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (arm == "lane") libs <- c("C:/Users/adf44/source/r/wt-silent-lib", libs)
.libPaths(libs)
suppressPackageStartupMessages(library(frmtmb))
cat("arm", arm, "frmtmb from", find.package("frmtmb"), "\n")
lln <- function(f) attr(logLik(f), "nobs")
try_s <- function(e) tryCatch({force(e); "ran"},
                              error = function(err) paste("ERROR:",
                                substr(conditionMessage(err), 1, 160)))

set.seed(5)
d6 <- data.frame(x = rnorm(80))
d6$y1 <- d6$x + rnorm(80)
d6$y2 <- d6$x + rnorm(80)
d6$s1 <- seq_len(80) <= 72
d6$s9 <- seq_len(80) > 8
f <- frm(bf(y1 | subset(s1) ~ x) + bf(y2 | subset(s1) ~ x) +
           set_rescor(FALSE), data = d6)
cat("A: logLik nobs", lln(f), " nobs()", nobs(f), " df.residual",
    df.residual(f), "\n")
# readers of n
ins <- tryCatch(insight::n_obs(f), error = function(e) conditionMessage(e))
cat("  insight::n_obs", ins, "\n")
cat("  stats::BIC", BIC(f), " AIC(k = log(80))", AIC(f, k = log(80)),
    " AIC(k = log(72))", AIC(f, k = log(72)), "\n")
bb <- tryCatch(bbmle::AICctab(f, nobs = NULL), error = function(e)
  conditionMessage(e))
cat("  bbmle::AICc", tryCatch(bbmle::AICc(f), error = function(e)
  conditionMessage(e)), "\n")
# rows 1..72 vs rows 9..80: equal count, different data
f9 <- frm(bf(y1 | subset(s9) ~ x) + bf(y2 | subset(s9) ~ 1) +
            set_rescor(FALSE), data = d6)
cat("E: anova(rows 1..72 fit, rows 9..80 fit):", try_s(anova(f, f9)), "\n")

# mi(): A reads mi(x) in y's formula, B does not; some rows miss both
set.seed(8)
n <- 120
dm <- data.frame(x = rnorm(n))
dm$y <- 0.5 * dm$x + rnorm(n)
dm$x[c(3, 10, 20, 30, 40)] <- NA
dm$y[c(3, 10, 50, 60)] <- NA
mA <- frm(bf(y | mi() ~ mi(x)) + bf(x | mi() ~ 1) + set_rescor(FALSE),
          data = dm)
mB <- frm(bf(y | mi() ~ 1) + bf(x | mi() ~ 1) + set_rescor(FALSE),
          data = dm)
cat("M: logLik nobs A", lln(mA), " B", lln(mB), " nobs() A", nobs(mA),
    " B", nobs(mB), "\n")
cat("   anova(B, A):", try_s(anova(mB, mA)), "\n")
cat("   logLik A", as.numeric(logLik(mA)), " B", as.numeric(logLik(mB)),
    "\n")
# univariate mi(): y missing on 4 rows, nested models
set.seed(9)
du <- data.frame(x = rnorm(n), z = rnorm(n))
du$y <- du$x + rnorm(n)
du$y[c(2, 4, 6, 8)] <- NA
u0 <- frm(bf(y | mi() ~ x), data = du)
u1 <- frm(bf(y | mi() ~ x + z), data = du)
cat("U: logLik nobs", lln(u0), lln(u1), " nobs()", nobs(u0),
    " df.residual", df.residual(u0), " anova:", try_s(anova(u0, u1)), "\n")
# mv where y2 reads mi(y1) and y1 missing: the latent row is read
set.seed(10)
dl <- data.frame(x = rnorm(n))
dl$y1 <- dl$x + rnorm(n)
dl$y2 <- 0.5 * dl$y1 + rnorm(n)
dl$y1[1:6] <- NA
dl$y2[1:3] <- NA
l1 <- frm(bf(y1 | mi() ~ x) + bf(y2 | mi() ~ mi(y1)) + set_rescor(FALSE),
          data = dl)
l0 <- frm(bf(y1 | mi() ~ x) + bf(y2 | mi() ~ 1) + set_rescor(FALSE),
          data = dl)
cat("L: logLik nobs with mi(y1) read", lln(l1), " without", lln(l0),
    " anova(l0, l1):", try_s(anova(l0, l1)), "\n")
