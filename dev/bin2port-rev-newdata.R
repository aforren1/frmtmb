# Reviewer: lane bin2port's wiener newdata defect, minimal, and its
# scope. Does predict(newdata =) (core predict_simulate()) read an
# addition term other than trunc(), rate() and thres(gr = ) from
# newdata? Seed 7. Log dev/bin2port-rev-out/newdata.txt.
.libPaths(c("C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.eam)})
cat("frmtmb", format(packageVersion("frmtmb")), "frmtmb.eam",
    format(packageVersion("frmtmb.eam")), "\n")
set.seed(7)
try_it <- function(label, expr) {
  r <- tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(label, ":", if (is.character(r)) r else
    paste("ran, dim", paste(dim(as.matrix(r)), collapse = " x ")), "\n")
  invisible(r)
}
cat("== wiener, dec() in newdata\n")
dat <- ddm_simulate(150, mu = 0.8, bs = 1.4, ndt = 0.3, bias = 0.5)
fit <- frm(bf(rt | dec(upper) ~ 1, bias = 0.5), family = wiener(),
           data = dat)
nd <- dat[1:3, , drop = FALSE]
try_it("fitted(newdata)", fitted(fit, newdata = nd))
try_it("predict(newdata, summary = FALSE)",
       predict(fit, newdata = nd, summary = FALSE, ndraws = 20))
try_it("predict(no newdata, summary = FALSE)",
       predict(fit, summary = FALSE, ndraws = 20))
try_it("simulate(newdata)", simulate(fit, nsim = 2, newdata = nd))

cat("== binomial, trials() in newdata: are the draws bounded by the",
    "newdata trials?\n")
n <- 200
db <- data.frame(x = rnorm(n), N = sample(5:15, n, TRUE))
db$y <- rbinom(n, db$N, plogis(0.3 + 0.5 * db$x))
fb <- frm(bf(y | trials(N) ~ x), family = binomial(), data = db)
ndb <- data.frame(x = c(0, 0), N = c(1000, 2))
pb <- try_it("predict(newdata N = 1000, 2, summary = FALSE)",
             predict(fb, newdata = ndb, summary = FALSE, ndraws = 200))
if (is.matrix(pb)) {
  cat("  column means of the draws:", format(colMeans(pb), digits = 4),
      "  max:", apply(pb, 2, max), "\n")
}
sb <- try_it("simulate(newdata)", simulate(fb, nsim = 200, newdata = ndb))
if (is.data.frame(sb)) cat("  simulate row means:",
                           format(rowMeans(as.matrix(sb)), digits = 4), "\n")
fbn <- try_it("fitted(newdata)", fitted(fb, newdata = ndb))
