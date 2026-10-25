# The brms 2.23.0 log-density identities for fcor() and sar(), printed:
# brms's Stan log_prob at frmtmb's estimates minus frmtmb's logLik (the
# constant brms leaves out) and the largest gradient there. The same
# checks run in test-fcor.R and frmtmb.spatial's test-sar.R; this prints
# the numbers. Data: test-fcor.R's fcor_data() (seed 1, n = 30) and
# test-sar.R's sar_data() (seed 2, 7 x 7 lattice, symmetric W).
# Run: Rscript dev/sarfcor-lpcheck.R > dev/sarfcor-out/lpcheck.txt
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-sarfcor/dev/stan-cache")
.libPaths(c("C:/Users/adf44/source/r/wt-sarfcor-lib",
            "C:/Users/adf44/source/r/rellib-r12",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({
  library(testthat)
  library(frmtmb.spatial)
})
wt <- "C:/Users/adf44/source/r/frmtmb-wt-sarfcor"
env <- new.env(parent = asNamespace("frmtmb"))
sys.source(file.path(wt, "extensions/frmtmb.spatial/tests/testthat",
                     "helper-brms.R"), envir = env)
options(frmtmb.brms_lp_report = TRUE)
lp <- get("brms_lp_check", env)
cat("lib:", find.package("frmtmb"), "|", find.package("frmtmb.spatial"),
    "\n")

set.seed(1)
n <- 30
V <- 0.6^abs(outer(seq_len(n), seq_len(n), "-"))
d <- data.frame(x = rnorm(n), g = factor(rep(seq_len(n / 5), each = 5)))
d$y <- 1 + 0.5 * d$x + as.vector(t(chol(V)) %*% rnorm(n)) * 0.8
show <- function(lab, r, expect) {
  cat(sprintf("%-14s lp - logLik %.12g  expected %.12g  diff %.3g  max|grad| %.3g\n",
              lab, r$measured_const, expect, r$measured_const - expect,
              r$max_grad))
}
f <- frm(y ~ x + fcor(V), data = d, data2 = list(V = V))
show("fcor", lp(brms::bf(y ~ x + fcor(V)), gaussian(), d, f,
                data2 = list(V = V)), 0)
f <- frm(bf(y ~ x + fcor(V), sigma ~ x), data = d, data2 = list(V = V))
show("fcor het", lp(brms::bf(y ~ x + fcor(V), sigma ~ x), gaussian(), d, f,
                    data2 = list(V = V)), 0)
f <- frm(y ~ x + fcor(V), family = student(), data = d, data2 = list(V = V))
show("fcor student", lp(brms::bf(y ~ x + fcor(V)), brms::student(), d, f,
                        data2 = list(V = V)), 0)

set.seed(2)
k <- 7
gr <- expand.grid(r = seq_len(k), c = seq_len(k))
n <- nrow(gr)
A <- outer(seq_len(n), seq_len(n), function(i, j) {
  as.numeric(abs(gr$r[i] - gr$r[j]) + abs(gr$c[i] - gr$c[j]) == 1)
})
W <- A / rowSums(A)
ds <- data.frame(x = rnorm(n))
IA <- diag(n) - 0.5 * W
ds$ylag <- as.vector(solve(IA, 1 + 0.5 * ds$x + rnorm(n, 0, 0.7)))
ds$yerr <- 1 + 0.5 * ds$x + as.vector(solve(IA, rnorm(n, 0, 0.7)))
dg <- 1 / sqrt(rowSums(A))
Ws <- A * outer(dg, dg)
c2 <- 0.5 * n * log(2 * pi)
f <- frm(ylag ~ x + sar(Ws), data = ds, data2 = list(Ws = Ws))
show("sar lag", lp(brms::bf(ylag ~ x + sar(Ws)), gaussian(), ds, f,
                   const = c2, data2 = list(Ws = Ws)), c2)
f <- frm(yerr ~ x + sar(Ws, type = "error"), data = ds,
         data2 = list(Ws = Ws))
show("sar error", lp(brms::bf(yerr ~ x + sar(Ws, type = "error")),
                     gaussian(), ds, f, const = c2,
                     data2 = list(Ws = Ws)), c2)
cp <- 0.5 * n * log(pi)
f <- frm(ylag ~ x + sar(Ws), family = student(), data = ds,
         data2 = list(Ws = Ws))
show("sar student", lp(brms::bf(ylag ~ x + sar(Ws)), brms::student(), ds, f,
                       const = cp, data2 = list(Ws = Ws)), cp)
# the nu the two student rows are evaluated at: lgamma((nu + n) / 2) -
# lgamma(nu / 2) as brms's Stan code writes it cancels there
f1 <- frm(y ~ x + fcor(V), family = student(), data = d,
          data2 = list(V = V))
cat("nu at the fcor student estimate:",
    format(frm_linpred(f1, dpar = "nu", type = "response")[1],
           digits = 4), "\n")
cat("nu at the sar student estimate:",
    format(frm_linpred(f, dpar = "nu", type = "response")[1],
           digits = 4), "\n")
