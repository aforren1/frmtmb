# Review lane osa: edges of the "cdf" route. References from stats::
# functions at mu = exp(X b) from fixef(), nothing from frmtmb internals.
#
#   Rscript dev/review-osa-edge.R [lib]
args <- commandArgs(TRUE)
lib1 <- if (length(args) >= 1) args[1] else "C:/Users/adf44/source/r/wt-osa-lib"
.libPaths(c(if (lib1 != "none") lib1, "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib: frmtmb", as.character(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "\n")
cat("BEGIN review-osa-edge\n")
osa <- function(fit, ...) {
  wn <- character(0)
  r <- tryCatch(withCallingHandlers(
    residuals(fit, type = "osa", ...)[, "Estimate"],
    warning = function(w) {
      wn <<- c(wn, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) e)
  attr(r, "warn") <- wn
  r
}
show <- function(label, r, ref_lo, ref_hi, rows) {
  if (inherits(r, "error")) {
    cat(sprintf("%-44s ERROR %s\n", label, substr(conditionMessage(r), 1, 150)))
    return(invisible())
  }
  for (i in rows) {
    cat(sprintf("%-44s row %3d r %10.4f  exact step [%10.4f, %10.4f]%s\n",
                label, i, r[i], qnorm(ref_lo[i]), qnorm(ref_hi[i]),
                if (is.finite(r[i]) && (r[i] < qnorm(ref_lo[i]) - 1e-6 ||
                                         r[i] > qnorm(ref_hi[i]) + 1e-6))
                  "  OUTSIDE" else ""))
  }
  if (length(attr(r, "warn"))) {
    cat(sprintf("%-44s WARN %s\n", label, substr(attr(r, "warn")[1], 1, 200)))
  }
}
U <- function(n) local({ set.seed(123); runif(n) })
mu_of <- function(fit, x) {
  b <- fixef(fit)[, "Estimate"]
  exp(b[[1]] + b[[2]] * x)
}

# 1. a large positive count outlier: true residual qnorm(F(y-1) + U p(y))
set.seed(401)
d <- data.frame(x = rnorm(100))
d$y <- rpois(100, exp(1 + 0.3 * d$x))
d$y[c(10, 20, 30)] <- c(14L, 22L, 60L)
f <- frm(y ~ x, family = poisson(), data = d)
m <- mu_of(f, d$x)
r <- osa(f)
show("poisson outliers y = 14, 22, 60", r, ppois(d$y - 1, m), ppois(d$y, m),
     c(10, 20, 30))
cat(sprintf("%-44s upper tails S(y - 1): %s\n", "",
            paste(format(ppois(d$y[c(10, 20, 30)] - 1, m[c(10, 20, 30)],
                               lower.tail = FALSE), digits = 3),
                  collapse = " ")))
# TMB's residual is qnorm(F(y) - U p); its upper tail is S(y - 1) - (1 - U) p
tru <- function(i) -qnorm(ppois(d$y[i] - 1, m[i], lower.tail = FALSE) -
                            (1 - U(100)[i]) * dpois(d$y[i], m[i]))
cat(sprintf("%-44s true residual -qnorm(S(y-1) - (1-U) p): row 10 %.4f, row 20: %.4f, row 30: %.3f\n",
            "", tru(10), tru(20), tru(30)))

# 2. a continuous outlier: Gamma, and gaussian under trunc() ("cdf") vs
#    without (fullGaussian)
set.seed(402)
dg <- data.frame(x = rnorm(100))
dg$y <- rgamma(100, 4, 4 / exp(0.3 + 0.2 * dg$x))
dg$y[5] <- 30
fgm <- frm(y ~ x, family = Gamma(link = "log"), data = dg)
rg <- osa(fgm)
shp <- as.numeric(frmtmb:::eval_dpars(fgm)[[1]]$shape)[1]
sg <- pgamma(30, shape = shp, rate = shp / mu_of(fgm, dg$x)[5],
             lower.tail = FALSE)
cat(sprintf("%-44s row 5 r %s; shape %.3f, pgamma upper tail %.3g -> true %.3f\n",
            "Gamma outlier y = 30", format(rg[5]), shp, sg, -qnorm(sg)))
# where a continuous residual stops being finite: y moved up row 5
for (yy in c(40, 50, 60, 80)) {
  dg2 <- dg
  dg2$y[5] <- yy
  fg2 <- frm(y ~ x, family = Gamma(link = "log"), data = dg2)
  s2 <- as.numeric(frmtmb:::eval_dpars(fg2)[[1]]$shape)[1]
  r2g <- osa(fg2)
  st <- pgamma(yy, shape = s2, rate = s2 / mu_of(fg2, dg2$x)[5],
               lower.tail = FALSE)
  cat(sprintf("%-44s y %g: r %s, true %.3f (S %.3g)%s\n", "Gamma outlier",
              yy, format(r2g[5], digits = 6), -qnorm(st), st,
              if (length(attr(r2g, "warn"))) "  WARN not finite" else ""))
}
cat(sprintf("%-44s WARN %s\n", "", paste(substr(attr(rg, "warn"), 1, 200),
                                          collapse = " | ")))
dn <- data.frame(x = rnorm(100))
dn$y <- 1 + 0.5 * dn$x + rnorm(100)
dn$y[7] <- 12
fn1 <- frm(y ~ x, family = gaussian(), data = dn)
fn2 <- frm(y | trunc(lb = -50) ~ x, family = gaussian(), data = dn)
r1 <- osa(fn1)
r2 <- osa(fn2)
cat(sprintf("%-44s row 7: fullGaussian %.3f, trunc(lb = -50) via cdf %s\n",
            "gaussian outlier at about 9.5 sd", r1[7], format(r2[7])))

# 3. lower tail below 1e-300: poisson mu about 800 and y = 0
set.seed(403)
d3 <- data.frame(x = seq(0, 1, length.out = 60))
d3$y <- rpois(60, exp(2 + 4.7 * d3$x))
d3$y[60] <- 0L
f3 <- frm(y ~ x, family = poisson(), data = d3)
m3 <- mu_of(f3, d3$x)
r3 <- osa(f3)
cat(sprintf("%-44s mu %.1f, log F(0) = %.1f; r %s; true in [-Inf, %.2f]\n",
            "poisson y = 0 at large mu", m3[60], ppois(0, m3[60], log.p = TRUE),
            format(r3[60]), qnorm(ppois(0, m3[60], log.p = TRUE), log.p = TRUE)))

# 4. support edges: binomial y = 0 and y = trials; trunc(ub) at y = ub
set.seed(404)
d4 <- data.frame(x = rnorm(120), k = 6L)
d4$y <- rbinom(120, 6, plogis(1.5 + 1.5 * d4$x))
f4 <- frm(y | trials(k) ~ x, family = binomial(), data = d4)
b4 <- fixef(f4)[, "Estimate"]
p4 <- plogis(b4[[1]] + b4[[2]] * d4$x)
r4 <- osa(f4)
edge <- which(d4$y %in% c(0, 6))
ok <- r4 >= qnorm(pbinom(d4$y - 1, 6, p4)) - 1e-8 &
  r4 <= qnorm(pbinom(d4$y, 6, p4)) + 1e-8
cat(sprintf("%-44s rows at 0 or 6: %d, outside their step: %d, nonfinite %d\n",
            "binomial support edges", length(edge), sum(!ok[edge]),
            sum(!is.finite(r4))))
set.seed(405)
d5 <- data.frame(x = rnorm(150))
y5 <- rpois(150, exp(1.6 + 0.3 * d5$x))
d5 <- d5[y5 >= 3 & y5 <= 7, , drop = FALSE]
d5$y <- y5[y5 >= 3 & y5 <= 7]
f5 <- frm(y | trunc(lb = 3, ub = 7) ~ x, family = poisson(), data = d5)
m5 <- mu_of(f5, d5$x)
Z <- ppois(7, m5) - ppois(2, m5)
lo5 <- (ppois(d5$y - 1, m5) - ppois(2, m5)) / Z
hi5 <- (ppois(d5$y, m5) - ppois(2, m5)) / Z
r5 <- osa(f5)
ok5 <- r5 >= qnorm(lo5) - 1e-8 & r5 <= qnorm(pmin(hi5, 1)) + 1e-8
cat(sprintf("%-44s rows %d (at 3: %d, at 7: %d), outside step %d, nonfinite %d\n",
            "poisson trunc(lb = 3, ub = 7)", nrow(d5), sum(d5$y == 3),
            sum(d5$y == 7), sum(!ok5), sum(!is.finite(r5))))

# 5. a truncation window of tiny mass: poisson mu about 0.05, lb = 3
set.seed(406)
d6 <- data.frame(x = rnorm(80, 0, 0.3))
m6 <- exp(-3 + 0.2 * d6$x)
d6$y <- qpois(runif(80, ppois(2, m6), 1), m6)
f6 <- frm(y | trunc(lb = 3) ~ x, family = poisson(), data = d6)
mm6 <- mu_of(f6, d6$x)
Z6 <- ppois(2, mm6, lower.tail = FALSE)
lo6 <- (ppois(d6$y - 1, mm6) - ppois(2, mm6)) / Z6
hi6 <- (ppois(d6$y, mm6) - ppois(2, mm6)) / Z6
# a cancellation-free reference: sums of the mass from 3
lo6b <- vapply(seq_len(80), function(i) {
  if (d6$y[i] <= 3) return(0)
  sum(dpois(3:(d6$y[i] - 1), mm6[i])) / Z6[i]
}, 0)
hi6b <- vapply(seq_len(80), function(i) sum(dpois(3:d6$y[i], mm6[i])) / Z6[i], 0)
r6 <- osa(f6)
ok6 <- r6 >= qnorm(lo6b) - 1e-6 & r6 <= qnorm(pmin(hi6b, 1)) + 1e-6
cat(sprintf("%-44s window mass %.2g; outside step %d of 80, nonfinite %d, y range %d-%d\n",
            "poisson trunc(lb = 3), mu about 0.05", median(Z6), sum(!ok6),
            sum(!is.finite(r6)), min(d6$y), max(d6$y)))

# 6. weight 0 rows with a group-level term: the fit ignores them, the
#    one-step route conditions on them
set.seed(407)
d7 <- data.frame(g = factor(rep(1:15, each = 8)), t = rep(1:8, 15),
                 x = rnorm(120))
d7$y <- rpois(120, exp(0.3 + 0.4 * d7$x + rnorm(15, 0, 1)[d7$g]))
d7$w <- ifelse(d7$t <= 4, 0, 1)
d7$y[d7$t <= 4] <- d7$y[d7$t <= 4] + 6L
f7 <- frm(y | weights(w) ~ x + (1 | g), family = poisson(), data = d7)
f7d <- frm(y ~ x + (1 | g), family = poisson(), data = d7[d7$t > 4, ])
r7 <- osa(f7)
r7d <- osa(f7d)
keep <- d7$t > 4
# r7d is on positions 1..60; the uniforms differ, so compare the mean
cat(sprintf(paste("%-44s estimates equal: %s; rows t > 4: mean r weighted",
                  "fit %.3f, mean r fit without the rows %.3f\n"),
            "weight 0 on rows t <= 4 (y + 6)",
            isTRUE(all.equal(f7$estimates$beta, f7d$estimates$beta,
                             tolerance = 1e-5)),
            mean(r7[keep]), mean(r7d)))
cat("END review-osa-edge\n")
