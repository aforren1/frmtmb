# Review lane osa, short check: weights on the fix build.
#   Rscript dev/review-osa-weights2.R
# 1. gaussian (1 | g), weight 2 on rows t <= 3, against the rows written
#    twice: a continuous residual has no uniform, so rows compare directly.
# 2. weight 0 under (1 | g): the other rows against the exact one-step
#    predictive with those rows left out (quadrature, stats:: only).
# 3. fractional weights (0.5 and 1.5) on a poisson and a Gamma with (1 | g).
.libPaths(c("C:/Users/adf44/source/r/wt-osa-lib",
            "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib: frmtmb", as.character(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "\n")
cat("BEGIN review-osa-weights2\n")
catch <- function(expr) {
  wn <- character(0)
  v <- withCallingHandlers(expr, warning = function(w) {
    wn <<- c(wn, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  attr(v, "warn") <- wn
  v
}
G <- 20L
K <- 10L
set.seed(801)
d <- data.frame(g = factor(rep(seq_len(G), each = K)), t = rep(1:K, G),
                x = rnorm(G * K))
u <- rnorm(G, 0, 1)[d$g]
d$yg <- 1 + 0.5 * d$x + u + rnorm(G * K)
d$w <- ifelse(d$t <= 3, 2, 1)
dd <- d[rep(seq_len(nrow(d)), d$w), ]
dd$last <- !duplicated(rep(seq_len(nrow(d)), d$w), fromLast = TRUE)
fw <- frm(yg | weights(w) ~ x + (1 | g), family = gaussian(), data = d)
fd <- frm(yg ~ x + (1 | g), family = gaussian(), data = dd)
rw <- catch(residuals(fw, type = "osa")[, "Estimate"])
rd <- residuals(fd, type = "osa")[, "Estimate"][dd$last]
cat(sprintf(paste("1. gaussian (1 | g) weight 2 vs duplicated (last copy):",
                  "max|diff| all rows %.3g, weight-2 rows %.3g; estimates",
                  "differ by %.2g\n"),
            max(abs(rw - rd)), max(abs(rw - rd)[d$w == 2]),
            max(abs(c(fw$estimates$beta, fw$estimates$theta) -
                      c(fd$estimates$beta, fd$estimates$theta)))))

# 2. weight 0, poisson
d$yp <- rpois(G * K, exp(0.2 + 0.5 * d$x + u))
d$w0 <- ifelse(d$t %in% c(2, 5), 0, 1)
d$yp[d$w0 == 0] <- d$yp[d$w0 == 0] + 6L
f0 <- frm(yp | weights(w0) ~ x + (1 | g), family = poisson(), data = d)
r0 <- catch(residuals(f0, type = "osa")[, "Estimate"])
U <- local({ set.seed(123); runif(nrow(d)) })
eta <- f0$estimates$beta[[1]] + f0$estimates$beta[[2]] * d$x
sdu <- exp(f0$estimates$theta[[1]])
s <- seq(-8, 8, length.out = 1601) * sdu
ez <- rep(NA_real_, nrow(d))
for (gi in split(seq_len(nrow(d)), d$g)) {
  w <- dnorm(s, 0, sdu)
  w <- w / sum(w)
  for (i in gi) {
    if (d$w0[i] == 0) next
    m <- exp(eta[i] + s)
    p <- dpois(d$yp[i], m)
    Fh <- ppois(d$yp[i], m)
    ez[i] <- qnorm(sum(w * Fh) - U[i] * sum(w * p))
    w <- w * p
    w <- w / sum(w)
  }
}
keep <- d$w0 == 1
cat(sprintf(paste("2. poisson (1 | g) weight 0 on t = 2, 5: NA on %d of %d",
                  "weight-0 rows; other rows vs exact without them: max %.3g",
                  "median %.3g; warnings: %s\n"),
            sum(is.na(r0[!keep])), sum(!keep), max(abs(r0 - ez)[keep]),
            median(abs(r0 - ez)[keep]),
            if (length(attr(r0, "warn"))) paste(attr(r0, "warn"), collapse = " | ") else "none"))

# 3. fractional weights
for (wv in c(0.5, 1.5)) {
  d$wf <- ifelse(d$t <= 3, wv, 1)
  fp <- frm(yp ~ x + (1 | g), family = poisson(), data = transform(d, yp = rpois(nrow(d), exp(0.2 + 0.5 * x + u))))
  dp <- transform(d, yp = rpois(nrow(d), exp(0.2 + 0.5 * x + u)),
                  yga = rgamma(nrow(d), 3, 3 / exp(0.3 + 0.3 * x + 0.5 * u)))
  f1 <- frm(yp | weights(wf) ~ x + (1 | g), family = poisson(), data = dp)
  f2 <- frm(yga | weights(wf) ~ x + (1 | g), family = Gamma(link = "log"),
            data = dp)
  for (ff in list(list("poisson", f1), list("Gamma", f2))) {
    r <- tryCatch(catch(residuals(ff[[2]], type = "osa")[, "Estimate"]),
                  error = function(e) e)
    if (inherits(r, "error")) {
      cat(sprintf("3. %-7s weight %.1f: ERROR %s\n", ff[[1]], wv,
                  conditionMessage(r)))
    } else {
      cat(sprintf(paste("3. %-7s weight %.1f on t <= 3: NA %d, nonfinite",
                        "non-NA %d, weighted rows mean %.3f sd %.3f,",
                        "others mean %.3f sd %.3f; warn: %s\n"),
                  ff[[1]], wv, sum(is.na(r)), sum(!is.finite(r) & !is.na(r)),
                  mean(r[d$wf != 1], na.rm = TRUE), sd(r[d$wf != 1], na.rm = TRUE),
                  mean(r[d$wf == 1]), sd(r[d$wf == 1]),
                  if (length(attr(r, "warn"))) paste(unique(substr(attr(r, "warn"), 1, 70)), collapse = " || ") else "none"))
    }
  }
}
cat("END review-osa-weights2\n")
