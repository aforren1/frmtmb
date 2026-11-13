# Lane gpstop: dev/fitedge-r2-gpfresh.R (lane fitedge, fix round) pointed at
# this round: arm base is rellib-r20, arm lane this lane's build, arm ref
# the lane library with the linear spelling and 12 theta starts (a theta
# start turns every warm start off, so the reference is build-free).
# Args: set (fresh | lane96), arm (base | lane | ref).
#   arm base  rellib-r20's default fit
#   arm lane  the lane's default fit
#   arm pick  the lane, with the scan run from the body's start as well, so the
#             objective at the start picks between the two (a cheap rule)
#   arm ref   lane library: the linear spelling (where one exists) and 12
#             theta starts on the body
# set fresh: dev/review-fitedge-designs.R, 8 designs x 3 curves x n 60, 120 x
#   seeds 1, 2. set lane96: the lane's own sweep (dev/fitedge-gpbody-sweep.R).
a <- commandArgs(TRUE)
set <- a[1]; arm <- a[2]
source("C:/Users/adf44/source/r/frmtmb-wt-gpstop/dev/gpstop-lib.R")
fe_libs(if (arm == "base") "base" else "lane")
RV_OUT <- file.path(FE_WT, "dev/gpstop-out")
if (arm == "pick") {
  f <- get("gp_start_scan", asNamespace("frmtmb"))
  b <- body(f)
  stopifnot(identical(b[[length(b)]],
                      quote(gp_start_scan_units(obj, template, frame, val,
                                                covered))))
  b[[length(b)]] <- quote(gp_start_scan_units(obj, template, frame, val))
  body(f) <- b
  assignInNamespace("gp_start_scan", f, "frmtmb")
  cat("PATCH pick: the scan also runs from the body's start\n")
}
if (set == "fresh") {
  source(file.path(FE_WT, "dev/review-fitedge-designs.R"))
  rows <- expand.grid(seed = 1:2, n = c(60L, 120L), ci = 1:3,
                      design = rv_designs, stringsAsFactors = FALSE)
  getd <- function(r) rv_mk(r$design, r$seed, r$n, r$ci)
  getsp <- function(r) rv_spell[[r$design]]
  lab <- function(r) c("sin", "wig", "exp")[r$ci]
} else {
  curves <- c("sin", "wig", "exp", "step")
  mk <- function(seed, n, ci, body) {
    set.seed(18000 + 100 * seed + n + ci)
    x <- sort(runif(n, 0, 10)); z <- rnorm(n); w <- runif(n, 1, 2)
    f <- switch(curves[ci], sin = sin(x), wig = 2 * sin(3 * x) + x,
                exp = exp(0.3 * x), step = sign(sin(1.5 * x)))
    y <- switch(body,
      A = f + 0.3 * z + 0.1 * rnorm(n),
      B = rpois(n, exp(1 + 0.5 * f / max(abs(f)) + 0.3 * z)),
      C = exp(0.5 + 0.5 * f / max(abs(f))) + 0.3 * z + 0.1 * rnorm(n),
      D = 1 + f * w + 0.1 * rnorm(n))
    data.frame(x, z, w, y)
  }
  spell <- list(
    A = list(nl = bf(y ~ a + b * z, a ~ gp(x), b ~ 1, nl = TRUE),
             lin = y ~ z + gp(x), fam = gaussian()),
    B = list(nl = bf(y ~ a + b * z, a ~ gp(x), b ~ 1, nl = TRUE),
             lin = y ~ z + gp(x), fam = poisson()),
    C = list(nl = bf(y ~ exp(a) + b * z, a ~ gp(x), b ~ 1, nl = TRUE),
             lin = NULL, fam = gaussian()),
    D = list(nl = bf(y ~ b + a * w, a ~ gp(x), b ~ 1, nl = TRUE),
             lin = NULL, fam = gaussian()))
  rows <- expand.grid(seed = 1:3, n = c(60L, 120L), ci = 1:4,
                      design = names(spell), stringsAsFactors = FALSE)
  getd <- function(r) mk(r$seed, r$n, r$ci, r$design)
  getsp <- function(r) spell[[r$design]]
  lab <- function(r) curves[r$ci]
}
# theta starts in each design's layout: (1 | g) first, a gp block is
# (log sdgp, log lscale), gp(by = ) holds one block per level
th_start <- function(design, sd0, rho) {
  switch(design, H = c(log(0.5), log(sd0), log(rho)),
         I = rep(c(log(sd0), log(rho)), 2), c(log(sd0), log(rho)))
}
out <- file.path(RV_OUT, sprintf("gpf-%s-%s.tsv", set, arm))
cat("design\tcurve\tn\tseed\tll\tconv\tnwarn\twclass\tsec\ttheta\n", file = out)
for (i in seq_len(nrow(rows))) {
  r <- rows[i, ]
  d <- getd(r)
  sp <- getsp(r)
  t0 <- proc.time()[[3]]
  th <- ""
  if (arm != "ref") {
    q <- fe_fit(sp$nl, data = d, family = sp$fam)
    ok <- !inherits(q$fit, "error")
    ll <- if (ok) as.numeric(logLik(q$fit)) else NA
    cv <- if (ok) q$fit$opt$convergence else NA
    wc <- if (ok) paste(unique(sub(":.*", "", q$warn)), collapse = ",") else
      paste("ERROR", substr(conditionMessage(q$fit), 1, 60))
    nw <- length(q$warn)
    if (ok) {
      th <- paste(round(q$fit$opt$par[names(q$fit$opt$par) == "theta"], 4),
                  collapse = ",")
    }
  } else {
    best <- -Inf
    if (!is.null(sp$lin)) {
      q <- fe_fit(sp$lin, data = d, family = sp$fam)
      if (!inherits(q$fit, "error")) best <- as.numeric(logLik(q$fit))
    }
    for (rho in c(0.15, 0.5, 1.5, 5)) for (sd0 in c(0.3, 1, 3)) {
      q <- fe_fit(sp$nl, data = d, family = sp$fam,
                  start = list(theta = th_start(r$design, sd0, rho)))
      if (!inherits(q$fit, "error")) {
        v <- as.numeric(logLik(q$fit))
        if (is.finite(v) && v > best) best <- v
      }
    }
    ll <- best; cv <- NA; nw <- NA; wc <- ""
  }
  cat(sprintf("%s\t%s\t%d\t%d\t%.6f\t%s\t%s\t%s\t%.2f\t%s\n", r$design, lab(r),
              r$n, r$seed, ll, cv, nw, wc, proc.time()[[3]] - t0, th),
      file = out, append = TRUE)
}
cat("DONE", out, "\n")
