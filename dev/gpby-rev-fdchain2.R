# Reviewer: where the chain route and the per-coefficient route of
# fit_fd_se() disagree past truncation error (me() + gp(), a new level),
# vary the step of each and see which one is stable. Also: disc ~ gp()
# fitted() error on both arms, and a categorical per-category gp().
.libPaths(c("C:/Users/adf44/source/r/wt-gpby-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
set.seed(11)
n <- 120
d <- data.frame(x = round(runif(n, 0, 5), 1), z = rnorm(n),
                g = factor(sample(letters[1:8], n, TRUE)),
                m = factor(sample(1:4, n, TRUE), ordered = TRUE))
d$sdx <- 0.2
d$xo <- d$z + rnorm(n, 0, 0.2)
lat <- sin(d$x) + 0.5 * d$z + rnorm(8)[as.integer(d$g)] * 0.5 +
  0.3 * as.integer(d$m) + rlogis(n)
d$y <- cut(lat, quantile(lat, c(0, 0.3, 0.6, 1)), include.lowest = TRUE,
           labels = FALSE)
d$y <- factor(d$y, ordered = TRUE)
d$yc <- factor(sample(c("A", "B", "C"), n, TRUE))
nd_new <- data.frame(x = c(5.6, 6.3, 2.55), z = c(0, 1, -1),
                     g = factor(c("a", "b", "c"), levels = levels(d$g)),
                     m = factor(c(1, 2, 3), levels = 1:4, ordered = TRUE),
                     xo = c(0, 1, -1), sdx = 0.2)
nd_lev <- nd_new
nd_lev$g <- factor(c("zz", "zz", "yy"))
se_of <- function(fit, newdata, ranl, chain, eps) {
  f <- function(fx) ns$fitted_point(fx, newdata, NULL, "response", NULL,
                                    NULL, ranl)
  want <- sort(unique(c(ns$smooth_b_idx(fit), ns$re_governed_b(fit))))
  used <- ns$re_used_b(fit, newdata, NULL, ranl)
  b_idx <- if (is.null(used)) want else intersect(want, used)
  bt <- ns$re_b_batches(fit, newdata, NULL, ranl, b_idx)
  ch <- ns$fd_eta_chain(fit, newdata, NULL, ranl, TRUE, b_idx)
  as.vector(ns$fit_fd_se(fit, f, eps = eps, b_idx = b_idx, b_batch = bt,
                         b_chain = if (chain) ch$chain))
}
for (cs in list(list("me(xo, sdx) + gp(x)", y ~ me(xo, sdx) + gp(x), nd_new,
                     FALSE),
                list("gp(x) + (1|g) new level", y ~ gp(x) + (1 | g), nd_lev,
                     TRUE),
                list("gp(x) + (1|g) seen level", y ~ gp(x) + (1 | g), nd_new,
                     FALSE))) {
  fit <- suppressWarnings(frm(bf(cs[[2]]), data = d, family = cumulative()))
  cat("==", cs[[1]], "\n")
  for (eps in c(1e-3, 1e-4, 1e-5, 1e-6)) {
    a <- se_of(fit, cs[[3]], cs[[4]], TRUE, eps)
    b <- se_of(fit, cs[[3]], cs[[4]], FALSE, eps)
    cat(sprintf("eps %.0e | chain %s | per-coef %s | max rel diff %.3e\n",
                eps, paste(sprintf("%.10f", a[1:3]), collapse = " "),
                paste(sprintf("%.10f", b[1:3]), collapse = " "),
                max(abs(a / b - 1))))
  }
}
cat("== disc ~ gp(x, k = 6)\n")
fd <- suppressWarnings(frm(bf(y ~ z, disc ~ 0 + gp(x, k = 6)), data = d,
                           family = cumulative()))
r <- tryCatch(fitted(fd, newdata = nd_new), error = function(e) e)
if (inherits(r, "error")) {
  cat("fitted(newdata) error:", conditionMessage(r), "\n")
  print(sys.calls)
  tr <- tryCatch(withCallingHandlers(fitted(fd, newdata = nd_new),
                                     error = function(e) {
                                       print(utils::head(sys.calls(), 40))
                                     }), error = function(e) NULL)
}
r2 <- tryCatch(fitted(fd), error = function(e) e)
cat("fitted() in sample:", if (inherits(r2, "error")) conditionMessage(r2)
    else "ok", "\n")
fz <- suppressWarnings(frm(bf(y ~ z, disc ~ 0 + x), data = d,
                           family = cumulative()))
r3 <- tryCatch(fitted(fz, newdata = nd_new), error = function(e) e)
cat("disc ~ 0 + x fitted(newdata):", if (inherits(r3, "error"))
  conditionMessage(r3) else "ok", "\n")
cat("== categorical muC ~ gp(x, k = 6)\n")
fc <- tryCatch(suppressWarnings(frm(bf(yc ~ z, muC ~ gp(x, k = 6)), data = d,
                                    family = categorical())),
               error = function(e) e)
if (inherits(fc, "error")) cat("refused:", conditionMessage(fc), "\n") else {
  for (eps in c(1e-4, 1e-5)) {
    a <- se_of(fc, nd_new, FALSE, TRUE, eps)
    b <- se_of(fc, nd_new, FALSE, FALSE, eps)
    cat(sprintf("eps %.0e | max rel diff chain vs per-coef %.3e\n", eps,
                max(abs(a / b - 1))))
  }
}
cat("DONE\n")
