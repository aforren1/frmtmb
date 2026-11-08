# Item 3: the outer Hessian of an exact gp() fit from the marginal's gradient
# against the z tape's, in one process, interleaved, minimum of `reps` rounds.
# Data as dev/review-ncgp-r2-prof.R: seed 20261103, y = sin(x) + 0.3 e, one row
# per position. Args: arm, positions, reps.
a <- commandArgs(TRUE)
arm <- a[1]; npos <- as.integer(a[2]); reps <- as.integer(a[3] %||% 3)
source("C:/Users/adf44/source/r/frmtmb-wt-fitedge/dev/fitedge-lib.R")
fe_libs(arm)
ns <- asNamespace("frmtmb")
set.seed(20261103)
# the review's design (dev/review-ncgp-r2-prof.R): `npos` rows on at most 500
# positions, since an exact gp() refuses more than 500
pos <- sort(runif(min(npos, 500L), 0, 10))
xr <- rep(pos, length.out = npos)
d <- data.frame(x = xr, y = sin(xr) + 0.3 * rnorm(npos))
t0 <- proc.time()[[3]]
if (identical(a[4], "warm")) {
  # past gp_marginal_max_rows the fit has no warm start; start it at the
  # marginal's own maximum so the descent is short, since what is timed is the
  # Hessian at the end (frm() itself is timed at 100 and 500)
  fr <- frm(y ~ gp(x), data = d, dry_run = "frame")
  m <- ns$gp_marginal_opt(ns$gp_marginal_setup(fr, max_rows = Inf), 3L, 4L)
  cat(sprintf("HT warm start from the marginal %.1f s\n", proc.time()[[3]] - t0))
  r <- fe_fit(y ~ gp(x), data = d,
              start = list(beta = m$beta, betad = m$betad, theta = m$theta))
} else {
  r <- fe_fit(y ~ gp(x), data = d)
}
tfit <- proc.time()[[3]] - t0
fit <- r$fit
cat(sprintf("HT %s pos %d frm %.2f s ll %.6f conv %d nwarn %d\n", arm, npos,
            tfit, as.numeric(logLik(fit)), fit$opt$convergence, length(r$warn)))
for (w in r$warn) cat("  W:", substr(w, 1, 160), "\n")
cat("HT SE", signif(sqrt(diag(as.matrix(vcov(fit, full = TRUE)))), 6), "\n")
has_new <- exists("gp_marginal_outer_gr", envir = ns, inherits = FALSE)
orig <- if (has_new) get("gp_marginal_outer_gr", envir = ns)
route <- function(new, k = 1L) {
  if (has_new) {
    f <- if (new) orig else function(fit, at) NULL
    unlockBinding("gp_marginal_outer_gr", ns)
    assign("gp_marginal_outer_gr", f, envir = ns)
    lockBinding("gp_marginal_outer_gr", ns)
  }
  t <- proc.time()[[3]]
  for (j in seq_len(k)) {
    rm(list = intersect(".frmtmb_outer_h", ls(fit$obj$env, all.names = TRUE)),
       envir = fit$obj$env)
    h <- ns$fit_outer_hessian(fit)
  }
  list(h = h, t = (proc.time()[[3]] - t) / k)
}
# blocks past 1.2 s on the slower route, so the 10 ms tick is under 1 percent
h1 <- route(FALSE)
k <- max(1L, ceiling(1.2 / max(h1$t, 0.01)))
tn <- to <- tc <- Inf
# the first z tape call is the control's other arm when blocks are of one call,
# so a slow size pays for three z tape Hessians and not four
if (k == 1L) tc <- h1$t
for (i in seq_len(reps)) {
  hn <- route(TRUE, k); tn <- min(tn, hn$t)
  ho <- route(FALSE, k); to <- min(to, ho$t)
  # the control: the z tape against itself, which must read 1
  if (k > 1L || i > 1L) tc <- min(tc, route(FALSE, k)$t)
}
cat(sprintf("HT %s pos %d control z tape %.4f s against %.4f s, ratio %.3f (block of %d)\n",
            arm, npos, tc, to, tc / to, k))
invisible(route(TRUE))
at <- hn$h$at
if (has_new) {
  g <- orig(fit, at)
  cat("HT marginal route taken:", is.function(g), "\n")
}
cat(sprintf("HT %s pos %d hessian marginal %.3f s  z tape %.3f s  ratio %.2f  (min of %d)\n",
            arm, npos, tn, to, to / tn, reps))
Hn <- hn$h$H; Ho <- ho$h$H
cat(sprintf("HT max |Hn - Ho| / max |Ho| %.3e ; max E new %.3e old %.3e\n",
            max(abs(Hn - Ho)) / max(abs(Ho)), max(hn$h$E), max(ho$h$E)))
sn <- sqrt(diag(solve(Hn))); so <- sqrt(diag(solve(Ho)))
cat("HT SE from H new", signif(sn, 7), "\nHT SE from H old", signif(so, 7),
    "\nHT max rel SE diff", signif(max(abs(sn / so - 1)), 3), "\n")
cat("DONE\n")
