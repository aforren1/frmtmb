# Short check, item 2: the one gaussian row where this build is far short and 0.78.0 is not: `logs` s 0.1 n 150 seed 12, and logs s 0.01 n 30. Is a cheap
# grid change enough? gp_marginal_opt() patched in process, as review-ncgp-r2-remedy.R. Then the arm on all of dev/ncgp-r2fresh.R's and my new curves if it closes them.
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
rev_libs("lane")
ns <- asNamespace("frmtmb")
orig <- ns$gp_marginal_opt
txt <- paste(deparse(orig), collapse = "\n")
cat(substr(txt, regexpr("lrs <-", txt), regexpr("lrs <-", txt) + 60), "\n")
i <- regexpr("o0 <- ", txt); cat(substr(txt, i - 120, i + 260), "\n")
mk <- function(lrs = NULL, nd = NULL) {
  t2 <- txt
  if (!is.null(lrs)) { t3 <- sub("lrs <- log(10^c(1, -1, -3, -6, -9))", paste0("lrs <- log(10^c(", lrs, "))"), t2, fixed = TRUE); stopifnot(t3 != t2); t2 <- t3 }
  if (!is.null(nd)) { t3 <- sub("seq_len(min(distinct, length(o0)))", paste0("seq_len(min(", nd, "L, length(o0)))"), t2, fixed = TRUE); stopifnot(t3 != t2); t2 <- t3 }
  f <- eval(parse(text = t2), envir = ns); environment(f) <- ns; f
}
set_opt <- function(f) { unlockBinding("gp_marginal_opt", ns); assign("gp_marginal_opt", f, envir = ns); lockBinding("gp_marginal_opt", ns) }
arms <- list(shipped = orig, lr7 = mk(lrs = "1, 0, -1, -2, -3, -6, -9"), dist6 = mk(nd = 6), both = mk(lrs = "1, 0, -1, -2, -3, -6, -9", nd = 6))
fl <- list(logs = function(x) log1p(x) * 2 + 0.2 * sin(8 * x), bump3 = function(x) exp(-(x - 2)^2 / 0.5) - exp(-(x - 5)^2 / 0.3) + exp(-(x - 8)^2 / 0.2), am = function(x) (1 + 0.5 * cos(x)) * sin(5 * x))
rows <- list(list("logs", 0.1, 150L, 12), list("logs", 0.01, 30L, 12), list("bump3", 1, 150L, 12), list("am", 1, 150L, 12))
for (rw in rows) {
  set.seed(71000 + rw[[4]] * 17 + rw[[3]]); n <- rw[[3]]
  x <- sort(round(runif(n, 0, 10), 4)); y <- fl[[rw[[1]]]](x) + rw[[2]] * rnorm(n)
  m0 <- cf_max(cf_setup(y, x), ngrid = 30L)
  out <- sprintf("RLG %s s=%g n %d cf %.4f sd %.3g rho %.3g sigma %.3g |", rw[[1]], rw[[2]], n, m0$ll, m0$sd, m0$rho, m0$sigma)
  for (a in names(arms)) {
    set_opt(arms[[a]]); t0 <- proc.time()[[3]]
    r <- rev_fit(y ~ gp(x), data.frame(y = y, x = x)); el <- proc.time()[[3]] - t0
    g <- rev_gp_read(r$fit)
    out <- paste0(out, sprintf(" %s %.3f sd %.3g rho %.3g %.2fs |", a, g$ll - m0$ll, g$sd, g$rho, el))
  }
  set_opt(orig); cat(out, "\n"); flush.console()
}
# the arm on two whole tables
fa <- list(
  bump3 = fl$bump3, chirp = function(x) sin(0.4 * x^2), tri = function(x) abs((x %% 3) - 1.5) * 2 + 0.3 * x, step2 = function(x) (x > 3) - 2 * (x > 7) + 0.2 * sin(6 * x),
  am = fl$am, logs = fl$logs, quad = function(x) (x - 5)^2 / 5 + cos(3 * x), off = function(x) 200 + 20 * cos(1.3 * x) + 0.5 * x)
tb <- c(rows = 0, s0 = 0, s1 = 0, t0 = 0, t1 = 0)
best_arm <- arms$dist6
for (seed in 11:12) for (fn in names(fa)) for (sg in c(1, 0.1, 1e-2, 1e-3)) for (n in c(30L, 70L, 150L)) {
  set.seed(71000 + seed * 17 + n); x <- sort(round(runif(n, 0, 10), 4)); y <- fa[[fn]](x) + sg * rnorm(n)
  set_opt(orig); t0 <- proc.time()[[3]]; r0 <- rev_fit(y ~ gp(x), data.frame(y = y, x = x)); e0 <- proc.time()[[3]] - t0
  set_opt(best_arm); t0 <- proc.time()[[3]]; r1 <- rev_fit(y ~ gp(x), data.frame(y = y, x = x)); e1 <- proc.time()[[3]] - t0; set_opt(orig)
  g0 <- rev_gp_read(r0$fit); g1 <- rev_gp_read(r1$fit)
  m <- cf_max(cf_setup(y, x), extra = list(c(log(g0$rho), 2 * (log(g0$sigma) - log(g0$sd))), c(log(g1$rho), 2 * (log(g1$sigma) - log(g1$sd)))), ngrid = 22L)
  tb <- tb + c(1, g0$ll - m$ll < -1e-3, g1$ll - m$ll < -1e-3, e0, e1)
  if (g1$ll - m$ll < -1e-3) cat(sprintf("RLG2 dist6 short %s s=%g n %d seed %d %.3g (shipped %.3g)\n", fn, sg, n, seed, g1$ll - m$ll, g0$ll - m$ll))
}
print(round(tb, 1))
cat("DONE\n")
