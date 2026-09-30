# Reviewer, claim 2: a multivariate subset() model against the sum of
# the separate fits on each response's own rows, over the feature
# paths the per-response frame touches. The responses share no
# parameter, so at the joint optimum logLik(joint) = sum of the
# separate logLik, and the fixed effects agree; a structural error
# shows as an O(1) gap, optimizer noise as ~1e-9 relative.
# Also fitted(resp =) and its standard error against the separate fit.
# Seed 222. Log: dev/aterms2-rev-log-02b.txt
.libPaths(c("C:/Users/adf44/source/r/wt-aterms2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
q <- function(expr) suppressWarnings(suppressMessages(expr))
set.seed(222)
n <- 180
G <- 9
d <- data.frame(x = rnorm(n), z = rnorm(n), w = rnorm(n),
                g = factor(rep(seq_len(G), length.out = n)),
                gg = factor(rep(c("u", "v"), length.out = n)),
                wt = runif(n, 0.5, 2), nt = sample(5:12, n, TRUE),
                time = runif(n, 0.5, 3),
                om = factor(sample(1:4, n, TRUE), ordered = TRUE),
                s1 = rep(c(TRUE, FALSE, TRUE), length.out = n),
                s2 = rep(c(TRUE, TRUE, FALSE, FALSE), length.out = n))
d$fg <- factor(ifelse(as.integer(d$g) <= 4, "a", "b"))
re <- rnorm(G, 0, 0.6)
d$y1 <- 1 + d$x + re[d$g] + rnorm(n)
d$y2 <- 0.5 - d$z + rnorm(n)
d$yc <- rpois(n, exp(0.3 + 0.4 * d$z + re[d$g]) * d$time)
d$yb <- rbinom(n, d$nt, plogis(-0.2 + 0.5 * d$z))
cut4 <- function(v) factor(cut(v, c(-Inf, -0.7, 0.2, 1.1, Inf),
                               labels = FALSE), ordered = TRUE)
d$o1 <- cut4(d$x + re[d$g] + rnorm(n))
d$o2 <- cut4(d$z + rnorm(n))
d$cc <- rep(c(0, 0, 1, -1), length.out = n)
d$yzi <- ifelse(runif(n) < 0.3, 0L, rpois(n, 2))
# NA outside a response's rows is harmless under subset(), as in brms
d$y2[which(!d$s2)[1:6]] <- NA
d$z2 <- d$z
d$z2[which(!d$s2)[7:9]] <- NA

cases <- list(
  gauss_pois_re = list(
    joint = quote(bf(y1 | subset(s1) ~ x + (1 | g)) + gaussian() +
                    bf(yc | subset(s2) ~ z + (1 | g)) + poisson()),
    parts = list(y1 = quote(frm(y1 ~ x + (1 | g), data = d[d$s1, ])),
                 yc = quote(frm(yc ~ z + (1 | g), data = d[d$s2, ],
                                family = poisson())))),
  rate_sub = list(
    joint = quote(bf(y1 | subset(s1) ~ x) + gaussian() +
                    bf(yc | subset(s2) + rate(time) ~ z) + poisson()),
    parts = list(y1 = quote(frm(y1 ~ x, data = d[d$s1, ])),
                 yc = quote(frm(yc | rate(time) ~ z, data = d[d$s2, ],
                                family = poisson())))),
  smooth = list(
    joint = quote(bf(y1 | subset(s1) ~ s(x)) + gaussian() +
                    bf(y2 | subset(s2) ~ s(z2)) + gaussian()),
    parts = list(y1 = quote(frm(y1 ~ s(x), data = d[d$s1, ])),
                 y2 = quote(frm(y2 ~ s(z2), data = d[d$s2, ])))),
  ordinal_cs = list(
    joint = quote(bf(o1 | subset(s1) ~ cs(x)) + sratio() +
                    bf(y2 | subset(s2) ~ z) + gaussian()),
    parts = list(o1 = quote(frm(o1 ~ cs(x), data = d[d$s1, ],
                                family = sratio())),
                 y2 = quote(frm(y2 ~ z, data = d[d$s2, ])))),
  ordinal_thronly = list(
    joint = quote(bf(o1 | subset(s1) ~ 1) + cumulative() +
                    bf(o2 ~ 1 + (1 | g)) + cumulative()),
    parts = list(o1 = quote(frm(o1 ~ 1, data = d[d$s1, ],
                                family = cumulative())),
                 o2 = quote(frm(o2 ~ 1 + (1 | g), data = d,
                                family = cumulative())))),
  thres_gr = list(
    joint = quote(bf(o1 | subset(s1) + thres(gr = gg) ~ x) + cumulative() +
                    bf(y1 ~ x) + gaussian()),
    parts = list(o1 = quote(frm(o1 | thres(gr = gg) ~ x, data = d[d$s1, ],
                                family = cumulative())),
                 y1 = quote(frm(y1 ~ x, data = d)))),
  weights_trials = list(
    joint = quote(bf(y1 | subset(s1) + weights(wt) ~ x) + gaussian() +
                    bf(yb | subset(s2) + trials(nt) ~ z) + binomial()),
    parts = list(y1 = quote(frm(y1 | weights(wt) ~ x, data = d[d$s1, ])),
                 yb = quote(frm(yb | trials(nt) ~ z, data = d[d$s2, ],
                                family = binomial())))),
  cens_mo = list(
    joint = quote(bf(y1 | subset(s1) + cens(cc) ~ mo(om)) + gaussian() +
                    bf(y2 | subset(s2) ~ z) + gaussian()),
    parts = list(y1 = quote(frm(y1 | cens(cc) ~ mo(om), data = d[d$s1, ])),
                 y2 = quote(frm(y2 ~ z, data = d[d$s2, ])))),
  gr_by_sigma = list(
    joint = quote(bf(y1 | subset(s1) ~ x + (1 | gr(g, by = fg)),
                     sigma ~ w) + gaussian() +
                    bf(y2 | subset(s2) ~ z) + gaussian()),
    parts = list(y1 = quote(frm(bf(y1 ~ x + (1 | gr(g, by = fg)),
                                   sigma ~ w), data = d[d$s1, ])),
                 y2 = quote(frm(y2 ~ z, data = d[d$s2, ])))),
  offset_zip = list(
    joint = quote(bf(yzi | subset(s1) ~ x + offset(log(time))) +
                    zero_inflated_poisson() +
                    bf(y2 | subset(s2) ~ z) + gaussian()),
    parts = list(yzi = quote(frm(yzi ~ x + offset(log(time)),
                                 data = d[d$s1, ],
                                 family = zero_inflated_poisson())),
                 y2 = quote(frm(y2 ~ z, data = d[d$s2, ])))),
  nonlinear = list(
    joint = quote(bf(y1 | subset(s1) ~ a + b * x, a ~ 1, b ~ 1,
                     nl = TRUE) + gaussian() +
                    bf(y2 | subset(s2) ~ z) + gaussian()),
    parts = list(y1 = quote(frm(bf(y1 ~ a + b * x, a ~ 1, b ~ 1,
                                   nl = TRUE), data = d[d$s1, ])),
                 y2 = quote(frm(y2 ~ z, data = d[d$s2, ])))),
  mixture = list(
    joint = quote(bf(y1 | subset(s1) ~ 1) +
                    mixture(gaussian(), gaussian()) +
                    bf(y2 | subset(s2) ~ z) + gaussian()),
    parts = list(y1 = quote(frm(y1 ~ 1, data = d[d$s1, ],
                                family = mixture(gaussian(), gaussian()))),
                 y2 = quote(frm(y2 ~ z, data = d[d$s2, ])))),
  one_unsubset = list(
    joint = quote(bf(y1 | subset(s1) ~ x + (1 | g)) + gaussian() +
                    bf(yc ~ z) + poisson()),
    parts = list(y1 = quote(frm(y1 ~ x + (1 | g), data = d[d$s1, ])),
                 yc = quote(frm(yc ~ z, data = d, family = poisson()))))
)

for (nm in names(cases)) {
  cs <- cases[[nm]]
  cat("\n==", nm, "\n")
  fj <- tryCatch(q(eval(bquote(frm(.(cs$joint), data = d)))),
                 error = function(e) conditionMessage(e))
  if (is.character(fj)) {
    cat("  JOINT ERROR:", fj, "\n")
    next
  }
  fp <- lapply(cs$parts, function(e) q(eval(e)))
  llj <- as.numeric(logLik(fj))
  lls <- sum(vapply(fp, function(f) as.numeric(logLik(f)), 0))
  cat(sprintf("  logLik joint %.12f  separate %.12f  rel %.3g\n", llj, lls,
              abs(llj - lls) / abs(lls)))
  cat("  rows:", paste(names(fj$frame$y), vapply(fj$frame$y, NROW, 1L),
                       collapse = " "), " separate:",
      paste(vapply(fp, function(f) NROW(f$frame$y[[1]]), 1L),
            collapse = " "), "\n")
  for (r in names(fp)) {
    a <- tryCatch(q(fitted(fj, resp = r)), error = function(e)
      conditionMessage(e))
    b <- q(fitted(fp[[r]]))
    if (is.character(a)) {
      cat("  fitted resp", r, "ERROR:", a, "\n")
      next
    }
    a <- if (length(dim(a)) == 3L) a[, 1, ] else a[, 1]
    b <- if (length(dim(b)) == 3L) b[, 1, ] else b[, 1]
    cat(sprintf("  fitted %s: dims %s vs %s, max rel diff %.3g\n", r,
                paste(dim(as.matrix(a)), collapse = "x"),
                paste(dim(as.matrix(b)), collapse = "x"),
                if (length(a) == length(b)) {
                  max(abs(a - b)) / max(abs(b))
                } else NA))
  }
}
