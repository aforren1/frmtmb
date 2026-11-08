# Lane arstep: the cause on the review's 25.9 row. Data as
# dev/review-ncgp-r3-ng.R, case gar, curve step, n 160, seed 1
# (set.seed(74000 + 31 * 1 + 160)). Args: arm [case n seed]
args <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-arstep/dev/arstep-lib.R")
arm <- args[1]
cs <- if (length(args) > 1) args[2] else "gar"
n <- if (length(args) > 2) as.integer(args[3]) else 160L
seed <- if (length(args) > 3) as.integer(args[4]) else 1L
rev_libs(arm)
ns <- asNamespace("frmtmb")
f <- function(x) 2 * (x > 5) - 1
set.seed(74000 + seed * 31 + n)
x <- sort(round(runif(n, 0, 10), 4))
eta <- f(x)
m <- switch(cs,
  student = list(d = data.frame(x = x, y = eta + 0.02 * rt(n, 5)), fml = y ~ gp(x), fam = student()),
  gar = list(d = data.frame(x = x, t = seq_len(n), y = eta + as.numeric(arima.sim(list(ar = 0.4), n, sd = 0.03))), fml = y ~ gp(x) + ar(time = t), fam = gaussian()))
fitit <- function(start = NULL) { a <- list(m$fml, m$d, family = m$fam); if (!is.null(start)) a$start <- start; do.call(rev_fit, a) }
show <- function(lbl, r) {
  e <- r$fit$estimates
  cat(sprintf("DIAG %s ll %.4f conv %d w %d | theta %s | betad %s | beta %s | other %s\n", lbl,
              as.numeric(logLik(r$fit)), r$fit$opt$convergence, length(r$warn),
              paste(signif(e$theta, 5), collapse = ","), paste(signif(e$betad, 5), collapse = ","),
              paste(signif(e$beta, 5), collapse = ","),
              paste(vapply(setdiff(names(e), c("theta", "betad", "beta", "b")), function(k) paste0(k, "=", paste(signif(e[[k]], 5), collapse = ",")), ""), collapse = " ")))
}
t0 <- proc.time()[[3]]
r0 <- fitit()
cat(sprintf("DIAG default wall %.2f s\n", proc.time()[[3]] - t0))
show("default", r0)
fr <- r0$fit$frame
cat("DIAG outer par names:", paste(names(r0$fit$obj$par), collapse = ","), "\n")

# 1. The working start: what gp_start_working() hands over, then each distinct
# local maximum of the working marginal over the length scale.
w <- ns$gp_start_working(fr)
cat(sprintf("DIAG working start theta %s (sdgp %.4g, lscale %.4g)\n",
            paste(signif(w$theta, 6), collapse = ","), exp(w$theta[1]), exp(w$theta[2])))
bk <- Filter(function(b) identical(b$covstruct, "gp"), fr$re_blocks)[[1]]
lp <- fr$linpreds[[ns$gp_block_lp(bk, fr)]]
yw <- ns$gp_working_y(fr, lp)
rb <- list(bk); rb[[1]]$theta_idx <- 1:2
fr2 <- list(re_blocks = rb, par_template = list(theta = numeric(2), beta = 0))
s <- ns$gp_marginal_data(yw, matrix(1, n, 1L), lp$Z, rb, rep(0, n), rep(1, n), TRUE, 1L, fr2)
d <- sqrt(max(Reduce(`+`, bk$aux_D2)))
# profile over the noise ratio at each length scale on a fine grid
prof <- function(lrho) {
  o <- optimize(function(lr) {
    v <- ns$gp_marginal_cell(s, ns$gp_marginal_theta(fr2, 0, list(`1` = exp(lrho))), lr)$value
    if (is.finite(v)) v else 1e300
  }, c(-40, 5))
  c(lrho = lrho, nll = o$objective, lr = o$minimum)
}
grid <- t(vapply(seq(log(d / 2^12), log(d * 4), length.out = 121), prof, numeric(3)))
loc <- which(diff(sign(diff(grid[, "nll"]))) > 0) + 1L
cat("DIAG working-marginal profile over log lscale (local minima of nll):\n")
for (i in loc) cat(sprintf("DIAG   local max: lscale %.4g  nll %.4f  sigma/sd %.3g\n",
                           exp(grid[i, "lrho"]), grid[i, "nll"], exp(grid[i, "lr"] / 2)))
cat(sprintf("DIAG   profile ends: lscale %.4g nll %.4f; lscale %.4g nll %.4f\n",
            exp(grid[1, 1]), grid[1, 2], exp(grid[nrow(grid), 1]), grid[nrow(grid), 2]))
mo <- ns$gp_marginal_opt(s)
cat(sprintf("DIAG gp_marginal_opt: value %.4f sdgp %.4g lscale %.4g\n", mo$value, exp(mo$theta[1]), exp(mo$theta[2])))

# 2. The full model's profile over the length scale: every other outer
# coordinate optimized on the fit's own tape at each lscale.
obj <- r0$fit$obj
pn <- names(obj$par)
ti <- which(pn == "theta")
il <- ti[2]
base_par <- r0$fit$opt$par
fullprof <- function(lrho, from) {
  free <- setdiff(seq_along(base_par), il)
  g <- function(q) { p <- from; p[free] <- q; p[il] <- lrho; v <- obj$fn(p); if (is.finite(v)) v else 1e300 }
  gg <- function(q) { p <- from; p[free] <- q; p[il] <- lrho; obj$gr(p)[free] }
  o <- tryCatch(nlminb(from[free], g, gg, control = list(iter.max = 300, eval.max = 400)), error = function(e) NULL)
  if (is.null(o)) return(list(v = NA, p = from))
  p <- from; p[free] <- o$par; p[il] <- lrho
  list(v = o$objective, p = p)
}
cat("DIAG full-model profile over lscale (warm-chained up and down from the default stop):\n")
lr0 <- base_par[il]
up <- seq(lr0, log(d * 2), by = 0.25); dn <- seq(lr0, log(d / 2^9), by = -0.25)
res <- list()
for (path in list(up, dn)) {
  from <- base_par
  for (lr in path) {
    pr <- fullprof(lr, from)
    if (is.finite(pr$v)) from <- pr$p
    res[[length(res) + 1]] <- c(lrho = lr, nll = pr$v, pr$p)
  }
}
res <- do.call(rbind, res); res <- res[order(res[, "lrho"]), ]
res <- res[!duplicated(round(res[, "lrho"], 6)), ]
for (i in seq_len(nrow(res))) cat(sprintf("DIAG   lscale %.4g  ll %.4f  par %s\n", exp(res[i, 1]), -res[i, 2], paste(signif(res[i, -(1:2)], 4), collapse = ",")))

# 3. The best of the review's 15 starts, and where it ends.
best <- r0; bv <- as.numeric(logLik(r0$fit))
for (rho0 in c(0.1, 0.3, 1, 3, 8)) for (sd0 in c(0.3, 1.5, 6)) {
  th <- r0$fit$estimates$theta
  th[bk$theta_idx[1]] <- log(sd0); th[bk$theta_idx[-1]] <- log(rho0)
  r2 <- fitit(list(theta = th))
  if (inherits(r2$fit, "error")) next
  v <- as.numeric(logLik(r2$fit))
  cat(sprintf("DIAG start sd %g lscale %g -> ll %.4f conv %d lscale %.4g sdgp %.4g\n", sd0, rho0, v, r2$fit$opt$convergence,
              exp(r2$fit$estimates$theta[bk$theta_idx[2]]), exp(r2$fit$estimates$theta[bk$theta_idx[1]])))
  if (is.finite(v) && v > bv) { bv <- v; best <- r2 }
}
show("best", best)
# 4. The default fit started from the working marginal's other local maxima.
for (i in loc) {
  th <- r0$fit$estimates$theta
  th[bk$theta_idx[1]] <- w$theta[1]; th[bk$theta_idx[-1]] <- grid[i, "lrho"]
  r2 <- fitit(list(theta = th))
  if (!inherits(r2$fit, "error")) show(sprintf("from-working-localmax-lscale-%.4g", exp(grid[i, "lrho"])), r2)
}
cat("DONE\n")
