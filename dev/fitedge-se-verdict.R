# Item 3: the SE verdicts of exact gp() fits with the outer Hessian from the
# marginal's gradient (new) against the z tape's (old), in one process: the old
# route is the same build with gp_marginal_outer_gr() returning NULL, which is
# 0.79.1's code path. Sets, as the ncgp review built them:
#   modes2  dev/review-ncgp-modes2.R, gaussian part: 7 curves x 3 noises x n 40,
#           120 x seeds 1:4, set.seed(31000 + 7 seed + n), 168 rows
#   fresh   dev/review-ncgp-r2-fresh.R: 8 curves x 4 noises x n 30, 70, 150 x
#           seeds 1:2, set.seed(55000 + 11 seed + n), 192 rows
#   other   dev/review-ncgp-r2-se.R's 4 fits (seed 60000 + s), the 10 spellings
#           of test-gp-start.R's "the marginal reads the tape's objective"
#           (seed 6101), and dev/review-ncgp-r2-edge2.R's data (seed 58001)
# Args: arm, set.
a <- commandArgs(TRUE)
arm <- a[1]; set <- a[2]
source("C:/Users/adf44/source/r/frmtmb-wt-fitedge/dev/fitedge-lib.R")
fe_libs(arm)
ns <- asNamespace("frmtmb")
orig <- get("gp_marginal_outer_gr", envir = ns)
# FE_FORCE=1 drops every gate after the structural one, so the marginal's
# Hessian is read wherever it can be built: what a gate has to choose from
if (identical(Sys.getenv("FE_FORCE"), "1")) {
  # the force runs predate gp_marginal_kappa_max; it is lifted here as well
  for (nm in intersect(c("gp_marginal_value_tol", "gp_marginal_kappa_max"),
                       ls(ns, all.names = TRUE))) {
    unlockBinding(nm, ns)
    assign(nm, Inf, envir = ns)
    lockBinding(nm, ns)
  }
  out_tag <- "force"
} else {
  out_tag <- arm
}
use <- function(new) {
  unlockBinding("gp_marginal_outer_gr", ns)
  assign("gp_marginal_outer_gr", if (new) orig else function(fit, at) NULL,
         envir = ns)
  lockBinding("gp_marginal_outer_gr", ns)
}
out <- file.path(FE_OUT, sprintf("severdict-%s-%s.tsv", set, out_tag))
cat("set\trow\tll_old\tll_new\tnan_old\tnan_new\twarn_same\twarn_prefix_same\tse_maxrel\tfit_se_maxrel\tvaldiff\troute\tsig_sd\tt_old\tt_new\n",
    file = out)
seen <- function(fit, w) {
  se <- tryCatch(sqrt(diag(as.matrix(vcov(fit, full = TRUE)))),
                 error = function(e) NA)
  list(ll = as.numeric(logLik(fit)), se = se,
       nan = paste(which(!is.finite(se)), collapse = ","), w = sort(w))
}
row <- function(lab, ...) {
  res <- list()
  for (new in c(FALSE, TRUE)) {
    use(new)
    t0 <- proc.time()[[3]]
    r <- fe_fit(...)
    tt <- proc.time()[[3]] - t0
    if (inherits(r$fit, "error")) {
      cat(lab, "ERROR", conditionMessage(r$fit), "\n")
      return(invisible(NULL))
    }
    # fitted standard errors at the first five rows, which read the joint
    # covariance the outer Hessian gives
    fse <- tryCatch(suppressWarnings(fitted(
      r$fit, newdata = list(...)$data[1:5, , drop = FALSE])[, 2L]),
      error = function(e) rep(NA_real_, 5L))
    res[[as.character(new)]] <- c(seen(r$fit, r$warn), t = tt,
                                  list(fit = r$fit, fse = fse))
  }
  use(TRUE)
  o <- res[["FALSE"]]; n <- res[["TRUE"]]
  fit <- n$fit
  s <- tryCatch(ns$gp_marginal_setup(fit$frame, max_rows = Inf),
                error = function(e) NULL)
  vd <- NA_real_
  if (!is.null(s)) {
    p <- list(beta = fit$estimates$beta[s$idx], theta = fit$estimates$theta)
    if (s$has_sigma) p$betad <- fit$estimates$betad
    v <- tryCatch(ns$gp_marginal_nll(s)(p), error = function(e) NA_real_)
    vd <- abs(v - fit$opt$objective) / max(1, abs(fit$opt$objective))
  }
  route <- is.function(orig(fit, fit$opt$par))
  ok <- is.finite(o$se) & is.finite(n$se)
  fo <- is.finite(o$fse) & is.finite(n$fse)
  # sigma over the first field's sd, the noise the marginal's covariance is
  # conditioned by
  bk <- Filter(function(b) identical(b$covstruct, "gp"), fit$frame$re_blocks)[[1]]
  sg <- if (length(fit$estimates$betad)) {
    exp(fit$estimates$betad[[1]] - fit$estimates$theta[[bk$theta_idx[1]]])
  } else NA_real_
  cat(sprintf("%s\t%s\t%.8f\t%.8f\t%s\t%s\t%s\t%s\t%.3e\t%.3e\t%.3e\t%s\t%.3e\t%.2f\t%.2f\n",
              set, lab, o$ll, n$ll, o$nan, n$nan, identical(o$w, n$w),
              identical(substr(o$w, 1, 60), substr(n$w, 1, 60)),
              if (any(ok)) max(abs(n$se[ok] / o$se[ok] - 1)) else NA,
              if (any(fo)) max(abs(n$fse[fo] / o$fse[fo] - 1)) else NA,
              vd, route, sg, o$t, n$t), file = out, append = TRUE)
}
if (set == "modes2") {
  fns <- list(mix = function(x) 2 * sin(3 * x) + x, mix2 = function(x) sin(4 * x) + 0.5 * x,
              mix3 = function(x) 3 * cos(2 * x) - 2 * x, mix4 = function(x) sin(5 * x) + x^2 / 10,
              two = function(x) sin(x) + 0.3 * sin(10 * x), sin3 = function(x) sin(3 * x),
              sin10 = function(x) sin(10 * x))
  for (seed in 1:4) for (fn in names(fns)) for (sg in c(0.3, 0.05, 1e-2)) for (n in c(40L, 120L)) {
    set.seed(31000 + seed * 7 + n)
    x <- sort(round(runif(n, 0, 10), 4))
    y <- fns[[fn]](x) + sg * rnorm(n)
    row(sprintf("%s_s%g_n%d_seed%d", fn, sg, n, seed), y ~ gp(x),
        data = data.frame(y = y, x = x))
  }
}
if (set == "fresh") {
  fns <- list(
    saw = function(x) (x %% 2) - 1 + 0.2 * x,
    abs3 = function(x) abs(sin(x)) * 3 - 0.5 * x,
    dbl = function(x) exp(-(x - 3)^2) + 2 * exp(-(x - 7)^2 / 0.1),
    gauss4 = function(x) sin(7 * x) * exp(-x / 5) + 0.3 * x,
    poly = function(x) (x - 2) * (x - 5) * (x - 8) / 10 + sin(4 * x),
    sq = function(x) sign(sin(1.5 * x)) + 0.1 * x,
    tan1 = function(x) atan(5 * (x - 5)) + cos(9 * x) / 2,
    big = function(x) 50 + 10 * sin(2.5 * x) + 3 * x)
  for (seed in 1:2) for (fn in names(fns)) for (sg in c(1, 0.1, 1e-2, 1e-3)) for (n in c(30L, 70L, 150L)) {
    set.seed(55000 + seed * 11 + n)
    x <- sort(round(runif(n, 0, 10), 4))
    y <- fns[[fn]](x) + sg * rnorm(n)
    row(sprintf("%s_s%g_n%d_seed%d", fn, sg, n, seed), y ~ gp(x),
        data = data.frame(y = y, x = x))
  }
}
if (set == "other") {
  for (cs in list(list(seed = 2, sg = 1, c0 = 1), list(seed = 2, sg = 1, c0 = 1000),
                  list(seed = 4, sg = 0.1, c0 = 1), list(seed = 4, sg = 0.1, c0 = 3))) {
    set.seed(60000 + cs$seed); n <- 50L
    x <- sort(round(runif(n, 0, 8), 4)); z <- rnorm(n)
    g <- factor(sample(letters[1:5], n, TRUE))
    y <- cs$c0 * (sin(2 * x) + 0.3 * z + cs$sg * rnorm(n))
    fml <- if (cs$seed == 2) y ~ z + gp(x) else y ~ gp(x, z)
    row(sprintf("r2se_seed%d_c%g", cs$seed, cs$c0), fml,
        data = data.frame(x, z, g, y))
  }
  set.seed(6101)
  n <- 60L
  d <- data.frame(x = runif(n, 0, 8), z = runif(n, 0, 8),
                  g = factor(sample(letters[1:5], n, TRUE)),
                  f = factor(sample(c("a", "b"), n, TRUE)),
                  w = runif(n, 0.5, 2), s = 0.05,
                  xr = rep(sort(runif(15L, 0, 8)), each = 4L))
  d$y <- sin(2 * d$x) + 0.5 * sin(d$z) + sin(2 * d$xr) + 0.1 * rnorm(n) +
    0.2 * rnorm(5L)[as.integer(d$g)]
  specs <- list(y ~ gp(x), y ~ z + gp(x), y ~ gp(x, by = f), y ~ gp(x, by = w),
                y ~ gp(x) + gp(z), y ~ z + gp(x) + (1 | g), y ~ gp(xr),
                y ~ gp(x, z, iso = FALSE), y ~ gp(x) + s(z),
                y | se(s) ~ gp(x), y | weights(w) ~ gp(x))
  for (i in seq_along(specs)) {
    row(sprintf("spec%d_%s", i, gsub("[ \t]", "", deparse1(specs[[i]]))),
        specs[[i]], data = d)
  }
  set.seed(58001); x <- sort(runif(60, 0, 10)); z <- rnorm(60)
  de <- data.frame(x, z, y = 2 * sin(3 * x) + x + 0.3 * z + 0.1 * rnorm(60))
  row("edge_nointercept", y ~ 0 + gp(x), data = de)
  row("edge_z", y ~ z + gp(x), data = de)
}
cat("DONE", out, "\n")
