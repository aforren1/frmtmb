# Reviewer, claim 2: factor levels under subset(). The frame is built on
# every row with drop.unused.levels, then cut to a response's rows, so a
# level present only outside the subset may survive into that
# response's design. Compared with the fit on the subsetted data and
# with brms's standata(). Seed 202. Log: dev/aterms2-rev-log-02a.txt
.libPaths(c("C:/Users/adf44/source/r/wt-aterms2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
try_ <- function(expr) tryCatch(expr, error = function(e) {
  cat("  ERROR:", conditionMessage(e), "\n"); NULL
})
set.seed(202)
n <- 90
d <- data.frame(x = rnorm(n), z = rnorm(n),
                f = factor(rep(c("a", "b", "c"), length.out = n)),
                g = factor(rep(letters[1:9], length.out = n)))
d$s <- d$f != "c"
d$s2 <- d$g != "i"
d$y <- 1 + d$x + as.integer(d$f) * 0.3 + rnorm(n)
d$y2 <- 0.5 * d$z + rnorm(n)
d$o <- factor(cut(d$x + rnorm(n), c(-Inf, -1, 0, 1, Inf), labels = FALSE),
              ordered = TRUE)
d$s3 <- d$o != "4"
# every level of g keeps rows under s5
d$s5 <- d$x > -0.3
stopifnot(all(table(d$g[d$s5]) > 0))

cat("== A. univariate, factor predictor, level c only outside subset\n")
fa <- try_(frm(y | subset(s) ~ x + f, data = d))
fb <- frm(y ~ x + f, data = d[d$s, ])
if (!is.null(fa)) {
  cat("  joint fixef:", format(fixef(fa)[, 1]), "\n")
  cat("  X colnames:", colnames(fa$frame$linpreds[[1]]$X), "\n")
  cat("  logLik subset", format(logLik(fa), digits = 15),
      " data[s,]", format(logLik(fb), digits = 15), "\n")
}
sd1 <- brms::standata(brms::bf(y | subset(s) ~ x + f), data = d)
cat("  brms X colnames:", colnames(sd1$X), " N", sd1$N, "\n")

cat("== B. multivariate, factor predictor on the subset response\n")
fm <- try_(frm(bf(y | subset(s) ~ x + f) + bf(y2 ~ z), data = d,
               family = gaussian()))
if (!is.null(fm)) {
  f1 <- frm(y ~ x + f, data = d[d$s, ])
  f2 <- frm(y2 ~ z, data = d)
  cat("  joint", format(logLik(fm), digits = 15), " separate sum",
      format(as.numeric(logLik(f1)) + as.numeric(logLik(f2)), digits = 15),
      "\n")
  print(fixef(fm)[, 1])
}
sd2 <- brms::standata(brms::bf(y | subset(s) ~ x + f) +
                        brms::bf(y2 ~ z) + brms::set_rescor(FALSE), data = d)
cat("  brms X_y colnames:", colnames(sd2$X_y), " N_y", sd2$N_y, "\n")

cat("== C. ordinal response, category 4 only outside subset\n")
fo <- try_(frm(o | subset(s3) ~ x, data = d, family = cumulative()))
fo2 <- frm(o ~ x, data = d[d$s3, ], family = cumulative())
if (!is.null(fo)) {
  cat("  thresholds joint:", length(fo$frame$par_template$tau_raw %||%
                                      numeric()), " logLik",
      format(logLik(fo), digits = 15), "\n")
}
cat("  thresholds data[s3,]:", length(fo2$frame$par_template$tau_raw %||%
                                        numeric()), " logLik",
    format(logLik(fo2), digits = 15), "\n")
sd3 <- brms::standata(brms::bf(o | subset(s3) ~ x), data = d,
                      family = brms::cumulative())
cat("  brms nthres", sd3$nthres, "\n")

cat("== D. mv, grouping level i only outside subset, (1 | g) unshared\n")
fg <- try_(frm(bf(y | subset(s2) ~ x + (1 | g)) + bf(y2 ~ z + (1 | g)),
               data = d, family = gaussian()))
if (!is.null(fg)) {
  f1 <- frm(y ~ x + (1 | g), data = d[d$s2, ])
  f2 <- frm(y2 ~ z + (1 | g), data = d)
  cat("  joint", format(logLik(fg), digits = 15), " separate sum",
      format(as.numeric(logLik(f1)) + as.numeric(logLik(f2)), digits = 15),
      "\n")
  cat("  ranef levels y:", vapply(ranef(fg), function(b) paste(dimnames(b)[[1]], collapse = ""), ""), "\n")
}
cat("== E. mv, same, |ID| shared: refusal expected\n")
fid <- try_(frm(bf(y | subset(s2) ~ x + (1 | p | g)) +
                  bf(y2 ~ z + (1 | p | g)), data = d, family = gaussian()))
cat("== F. mv, |ID| shared, subset keeps every level: no refusal\n")
fid2 <- try_(frm(bf(y | subset(s5) ~ x + (1 | p | g)) +
                   bf(y2 ~ z + (1 | p | g)), data = d, family = gaussian()))
if (!is.null(fid2)) cat("  fitted, logLik", format(logLik(fid2),
                                                   digits = 15), "\n")
cat("== G. mv, |ID| where BOTH responses subset to the same level set\n")
d$s4 <- d$g != "i" & d$x > -1
fid3 <- try_(frm(bf(y | subset(s2) ~ x + (1 | p | g)) +
                   bf(y2 | subset(s4) ~ z + (1 | p | g)), data = d,
                 family = gaussian()))
if (!is.null(fid3)) cat("  fitted, logLik", format(logLik(fid3),
                                                   digits = 15), "\n")
cat("== H. factor g with an unused level in data (droplevels issue)\n")
d$g2 <- factor(as.character(d$g), levels = c(letters[1:9], "zz"))
fh <- try_(frm(bf(y | subset(s5) ~ x + (1 | p | g2)) +
                 bf(y2 ~ z + (1 | p | g2)), data = d, family = gaussian()))
if (!is.null(fh)) cat("  fitted, logLik", format(logLik(fh), digits = 15),
                      "\n")
