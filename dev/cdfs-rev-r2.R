# Reviewer, fix round 2: (a) cens() inside a two-sided trunc() against a
# hand-rolled likelihood; (b) simulate() on a window with almost no
# mass, inside the window, without hanging, distributed as the
# truncated CDF; (c) the weibull fuzz spec fills every row; (d) brms's
# truncated-mean coverage.
.libPaths(c("C:/Users/adf44/source/r/wt-cdfs-rev-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(testthat); library(frmtmb)})
options(frmtmb.notices = FALSE)
cat("lib:", find.package("frmtmb"), "\n")
lse <- function(v) { m <- max(v); m + log(sum(exp(v - m))) }
cat("== (a) cens() inside trunc(lb = 2, ub = 9), poisson; and gaussian trunc(-1, 3)\n")
set.seed(4)
n <- 200
x <- rnorm(n)
y <- rpois(3 * n, exp(1.3 + 0.3 * rep(x, 3))); xx <- rep(x, 3)
k <- which(y >= 2 & y <= 9)[1:n]
d <- data.frame(x = xx[k], y = y[k], cc = sample(c(-1, 0, 0, 1, 2), n, TRUE))
d$y2 <- pmin(d$y + 2, 9)
d$cc[d$cc == 2 & d$y2 <= d$y] <- 0
d$cc[d$cc == 1 & d$y >= 9] <- 0
d$cc[d$cc == -1 & d$y <= 2] <- 0
ob <- frm(y | cens(cc, y2) + trunc(lb = 2, ub = 9) ~ x, family = poisson(), data = d,
          dry_run = "objective")$obj
p <- ob$par + c(0.05, -0.03)
mu <- exp(p[1] + p[2] * d$x)
win <- function(a, b, m) lse(dpois(a:b, m, log = TRUE))   # P(a <= Y <= b)
ll <- vapply(seq_len(nrow(d)), function(i) {
  m <- mu[i]; yi <- d$y[i]
  num <- switch(as.character(d$cc[i]), "0" = dpois(yi, m, log = TRUE),
                "1" = win(yi + 1, 9, m), "-1" = win(2, yi, m),
                "2" = win(yi + 1, d$y2[i], m))
  num - win(2, 9, m)
}, 0)
cat(sprintf("poisson: -fn %.10f  hand (window composition) %.10f  codes %s\n", -ob$fn(p), sum(ll),
            paste(table(d$cc), collapse = "/")))
g <- data.frame(x = x, y = 1 + 0.5 * x + rnorm(n))
g <- g[g$y > -1 & g$y < 3, ]
g$cc <- sample(c(-1, 0, 1, 2), nrow(g), TRUE)
g$y2 <- pmin(g$y + 0.5, 2.99)
g$cc[g$cc == 2 & g$y2 <= g$y] <- 0
ob <- frm(y | cens(cc, y2) + trunc(lb = -1, ub = 3) ~ x, family = gaussian(), data = g,
          dry_run = "objective")$obj
p <- ob$par
m <- p[1] + p[2] * g$x; s <- exp(p[3])
P <- function(a, b) pnorm(b, m, s) - pnorm(a, m, s)
num <- ifelse(g$cc == 0, dnorm(g$y, m, s, log = TRUE),
       ifelse(g$cc == 1, log(P(g$y, 3)), ifelse(g$cc == -1, log(P(-1, g$y)), log(P(g$y, g$y2)))))
cat(sprintf("gaussian: -fn %.10f  hand %.10f\n", -ob$fn(p), sum(num - log(P(-1, 3)))))

cat("== (b) simulate() from windows with almost no mass\n")
one <- function(lab, f, fam, dd, lo, hi, cdf) {
  fit <- suppressWarnings(frm(f, family = fam, data = dd))
  t0 <- proc.time()[["elapsed"]]
  s <- tryCatch(simulate(fit, nsim = 200, seed = 3), error = function(e) conditionMessage(e))
  el <- proc.time()[["elapsed"]] - t0
  if (is.character(s)) { cat(lab, "ERROR:", s, "\n"); return(invisible()) }
  v <- as.numeric(as.matrix(s))
  cat(sprintf("%-46s %.1f s | draws %d | outside [%g, %g]: %d | NA %d\n", lab, el, length(v), lo, hi,
              sum(v < lo | v > hi, na.rm = TRUE), sum(is.na(v))))
  invisible(list(fit = fit, v = v))
}
# poisson mu 2 in [40, 45]: P(window) about 4e-31; data must lie in it,
# and the fitted mean has to stay near 2, so mu is held through an
# offset-free model fitted to a few window rows plus a fixed intercept
# start is not possible; instead the fit runs on data in the window and
# simulate() is called from a model whose mean is then forced to 2
dp <- data.frame(y = c(40, 41, 40, 42, 45, 40, 43, 41))
fitp <- frm(y | trunc(lb = 40, ub = 45) ~ 1, family = poisson(), data = dp)
cat("poisson [40, 45] fit mu:", exp(fitp$opt$par), "\n")
ns <- asNamespace("frmtmb")
fam <- ns$as_frmtmb_family(poisson())
t0 <- proc.time()[["elapsed"]]
sv <- tryCatch(ns$sim_response(fam, list(mu = rep(2, 5000)),
                               list(trunc_lb = rep(40, 5000), trunc_ub = rep(45, 5000)), 5000),
               error = function(e) conditionMessage(e))
cat(sprintf("sim_response poisson mu 2 in [40, 45], 5000 rows: %.1f s\n", proc.time()[["elapsed"]] - t0))
if (is.character(sv)) cat("ERROR:", sv, "\n") else {
  tb <- table(factor(sv, levels = 40:45))
  pr <- dpois(40:45, 2) / sum(dpois(40:45, 2))
  cat("counts 40..45:", tb, " | expected:", round(5000 * pr, 1), " | outside:", sum(sv < 40 | sv > 45),
      " | chisq p", signif(suppressWarnings(chisq.test(tb, p = pr))$p.value, 3), "\n")
}
for (cs in list(c(12, 15), c(100, 110), c(0, 0), c(0, 1))) {
  sv <- tryCatch(ns$sim_response(fam, list(mu = rep(2, 2000)),
                                 list(trunc_lb = rep(cs[1], 2000), trunc_ub = rep(cs[2], 2000)), 2000),
                 error = function(e) conditionMessage(e))
  if (is.character(sv)) cat("[", cs, "] ERROR:", sv, "\n") else
    cat(sprintf("poisson mu 2 in [%g, %g]: outside %d, mean %.4f (truncated mean %.4f)\n", cs[1], cs[2],
                sum(sv < cs[1] | sv > cs[2]), mean(sv),
                sum((cs[1]:cs[2]) * dpois(cs[1]:cs[2], 2)) / sum(dpois(cs[1]:cs[2], 2))))
}
famg <- ns$as_frmtmb_family(gaussian())
sv <- ns$sim_response(famg, list(mu = rep(0, 4000), sigma = rep(1, 4000)),
                      list(trunc_lb = rep(9, 4000), trunc_ub = rep(9.5, 4000)), 4000)
Ft <- function(q) (pnorm(q, lower.tail = FALSE, log.p = FALSE) - pnorm(9, lower.tail = FALSE)) /
  (pnorm(9.5, lower.tail = FALSE) - pnorm(9, lower.tail = FALSE))
cat(sprintf("gaussian in [9, 9.5]: outside %d, KS p %.3g\n", sum(sv < 9 | sv > 9.5),
            suppressWarnings(ks.test(sv, function(q) -Ft(q))$p.value)))
famG <- ns$as_frmtmb_family(Gamma(link = "log"))
sv <- ns$sim_response(famG, list(mu = rep(1, 4000), shape = rep(3, 4000)), list(trunc_lb = rep(8, 4000)), 4000)
cat(sprintf("Gamma(shape 3, mean 1) above 8: outside %d, KS p %.3g\n", sum(sv < 8),
            ks.test(sv, function(q) 1 - pgamma(q, 3, 3, lower.tail = FALSE) / pgamma(8, 3, 3, lower.tail = FALSE))$p.value))

cat("== (c) the weibull fuzz spec 20446531\n")
testthat::source_test_helpers("C:/Users/adf44/source/r/frmtmb-wt-cdfs/tests/testthat", env = globalenv())
plan <- fuzz_plan(seed = 20260901L, size = 300L)
sp <- as.list(plan[plan$seed == 20446531, , drop = FALSE][1, ])
d <- fuzz_data(sp)
fit <- suppressWarnings(eval(parse(text = fuzz_call_text(sp))))
s <- simulate(fit, nsim = 20, seed = 1)
v <- as.matrix(s)
cat(sprintf("rows %d x %d, NA %d, below 0.5: %d, min %.4f\n", nrow(v), ncol(v), sum(is.na(v)),
            sum(v < 0.5, na.rm = TRUE), min(v, na.rm = TRUE)))

cat("== (d) brms 2.23.0's truncated means\n")
bn <- asNamespace("brms")
cat("posterior_epred_trunc_* in brms:",
    paste(sub("posterior_epred_trunc_", "", grep("^posterior_epred_trunc_", ls(bn), value = TRUE)), collapse = ", "), "\n")
