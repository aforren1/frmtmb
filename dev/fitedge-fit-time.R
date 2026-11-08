# Item 3: whole frm() at `positions`, the marginal's Hessian against the z
# tape's (gp_marginal_outer_gr() patched to NULL, 0.79.1's path), interleaved
# in one process, `reps` rounds, minimum of each. Data as
# dev/fitedge-hess-time.R (seed 20261103). Args: arm, positions, reps.
a <- commandArgs(TRUE)
arm <- a[1]; npos <- as.integer(a[2]); reps <- as.integer(a[3])
source("C:/Users/adf44/source/r/frmtmb-wt-fitedge/dev/fitedge-lib.R")
fe_libs(arm)
ns <- asNamespace("frmtmb")
orig <- get("gp_marginal_outer_gr", envir = ns)
use <- function(new) {
  unlockBinding("gp_marginal_outer_gr", ns)
  assign("gp_marginal_outer_gr", if (new) orig else function(fit, at) NULL,
         envir = ns)
  lockBinding("gp_marginal_outer_gr", ns)
}
set.seed(20261103)
# the review's design (dev/review-ncgp-r2-prof.R): `npos` rows on at most 500
# positions, since an exact gp() refuses more than 500
pos <- sort(runif(min(npos, 500L), 0, 10))
xr <- rep(pos, length.out = npos)
d <- data.frame(x = xr, y = sin(xr) + 0.3 * rnorm(npos))
one <- function(new) {
  use(new)
  t0 <- proc.time()[[3]]
  f <- suppressWarnings(frm(y ~ gp(x), data = d))
  s <- summary(f)
  list(t = proc.time()[[3]] - t0, ll = as.numeric(logLik(f)),
       se = sqrt(diag(as.matrix(vcov(f, full = TRUE)))))
}
tn <- to <- Inf
# a first fit in a process pays for loading; the order alternates by round
invisible(one(FALSE))
for (i in seq_len(reps)) {
  if (i %% 2L) {
    rn <- one(TRUE); ro <- one(FALSE)
  } else {
    ro <- one(FALSE); rn <- one(TRUE)
  }
  tn <- min(tn, rn$t)
  to <- min(to, ro$t)
}
cat(sprintf("FT %s pos %d frm()+summary(): marginal Hessian %.2f s, z tape %.2f s, ratio %.2f (min of %d)\n",
            arm, npos, tn, to, to / tn, reps))
cat(sprintf("FT logLik equal: %s; max rel SE diff %.2e\n", identical(rn$ll, ro$ll),
            max(abs(rn$se / ro$se - 1))))
cat("DONE\n")
