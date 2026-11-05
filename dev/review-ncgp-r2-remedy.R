# Fix round 1: what would close the fresh rows. gp_marginal_opt() patched in this process: arm A, rho on
# factors of 2 instead of 4 over the same range, and the best cell of each DISTINCT length scale instead
# of the best three cells overall, keep 4. Same curves and seeds as review-ncgp-r2-fresh.R. Args: seeds
args <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
seeds <- eval(parse(text = args[1]))
rev_libs("lane")
ns <- asNamespace("frmtmb")
src <- deparse(ns$gp_marginal_opt)
txt <- paste(src, collapse = "\n")
txt2 <- sub("for (r in d/4^(0:5)) {", "for (r in d/2^(0:10)) {", txt, fixed = TRUE)
stopifnot(!identical(txt, txt2))
old_ord <- "ord <- order(vals)[seq_len(min(keep, length(vals)))]"
stopifnot(grepl(old_ord, txt2, fixed = TRUE))
new_ord <- paste("{",
  "rk <- vapply(cells, function(cc) paste(signif(unlist(cc$rho), 6), collapse = '_'), '')",
  "o0 <- order(vals); o0 <- o0[!duplicated(rk[o0])]",
  "ord <- o0[seq_len(min(keep, length(o0)))]",
  "}", sep = "\n")
txt2 <- sub(old_ord, new_ord, txt2, fixed = TRUE)
patched <- eval(parse(text = txt2), envir = ns)
environment(patched) <- ns
orig_opt <- ns$gp_marginal_opt
orig_max <- ns$gp_marginal_max
set_opt <- function(f) { unlockBinding("gp_marginal_opt", ns); assign("gp_marginal_opt", f, envir = ns); lockBinding("gp_marginal_opt", ns) }
set_keep <- function(k) { f <- orig_max; formals(f)$keep <- k; environment(f) <- ns; unlockBinding("gp_marginal_max", ns); assign("gp_marginal_max", f, envir = ns); lockBinding("gp_marginal_max", ns) }
fns <- list(
  saw = function(x) (x %% 2) - 1 + 0.2 * x, abs3 = function(x) abs(sin(x)) * 3 - 0.5 * x, dbl = function(x) exp(-(x - 3)^2) + 2 * exp(-(x - 7)^2 / 0.1),
  gauss4 = function(x) sin(7 * x) * exp(-x / 5) + 0.3 * x, poly = function(x) (x - 2) * (x - 5) * (x - 8) / 10 + sin(4 * x), sq = function(x) sign(sin(1.5 * x)) + 0.1 * x,
  tan1 = function(x) atan(5 * (x - 5)) + cos(9 * x) / 2, big = function(x) 50 + 10 * sin(2.5 * x) + 3 * x)
tally <- list(ship = c(0, 0, 0, 0), arm4 = c(0, 0, 0, 0))
for (seed in seeds) for (fn in names(fns)) for (sg in c(1, 0.1, 1e-2, 1e-3)) for (n in c(30L, 70L, 150L)) {
  set.seed(55000 + seed * 11 + n)
  x <- sort(round(runif(n, 0, 10), 4))
  y <- fns[[fn]](x) + sg * rnorm(n)
  d <- data.frame(y = y, x = x)
  set_opt(orig_opt); set_keep(3L)
  t0 <- proc.time()[[3]]; r1 <- rev_fit(y ~ gp(x), d); t1 <- proc.time()[[3]] - t0
  set_opt(patched); set_keep(4L)
  t0 <- proc.time()[[3]]; r2 <- rev_fit(y ~ gp(x), d); t2 <- proc.time()[[3]] - t0
  set_opt(orig_opt); set_keep(3L)
  g1 <- rev_gp_read(r1$fit); g2 <- rev_gp_read(r2$fit)
  s <- cf_setup(y, x)
  m <- cf_max(s, extra = list(c(log(g1$rho), 2 * (log(g1$sigma) - log(g1$sd))), c(log(g2$rho), 2 * (log(g2$sigma) - log(g2$sd)))), ngrid = 24L)
  a1 <- g1$ll - m$ll; a2 <- g2$ll - m$ll
  sil1 <- a1 < -1e-3 && g1$conv == 0 && !length(r1$warn); sil2 <- a2 < -1e-3 && g2$conv == 0 && !length(r2$warn)
  tally$ship <- tally$ship + c(1, a1 < -1e-3, sil1, t1); tally$arm4 <- tally$arm4 + c(1, a2 < -1e-3, sil2, t2)
  cat(sprintf("RRM %s s=%g n %d seed %d ship %.2e conv %d w %d %.2fs | arm %.2e conv %d w %d %.2fs\n", fn, sg, n, seed, a1, g1$conv, length(r1$warn), t1, a2, g2$conv, length(r2$warn), t2))
  flush.console()
}
cat(sprintf("RRMSUM ship rows %d short %d silent %d sec %.0f | arm rows %d short %d silent %d sec %.0f\n", tally$ship[1], tally$ship[2], tally$ship[3], tally$ship[4], tally$arm4[1], tally$arm4[2], tally$arm4[3], tally$arm4[4]))
cat("DONE\n")
