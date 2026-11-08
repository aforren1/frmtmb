# Review of lane fitedge, item 3. Lane library only. Args: part.
#   near  fits whose condition bound lands near the gate (1e4 to 1e8): the
#         SE verdicts with the route on and off (two fits each), and on the
#         route-off fit the outer Hessian three ways: z tape, the gate as
#         shipped, the gate lifted. Also the true 2-norm condition of V.
#   fall  fits the gate must hand back to the z tape (or admit), each with
#         the route on and off: no error, the same SEs, NaN sets and warnings
#   time  frm() + vcov() at 100 and 500 positions, the route on and off,
#         interleaved in one process; and the outer Hessian alone
# Seeds: near set.seed(71000 + 13 seed + 7 i); fall set.seed(72000 + i).
a <- commandArgs(TRUE)
part <- a[1]
source("C:/Users/adf44/source/r/frmtmb-wt-fitedge/dev/fitedge-lib.R")
# RV_LIB names another library for the lane arm (the merged scratch tree)
FE_LIB <- Sys.getenv("RV_LIB", FE_LIB)
RV_TAG <- Sys.getenv("RV_TAG", "")
fe_libs("lane")
RV_OUT <- file.path(FE_WT, "dev/review-fitedge-out")
ns <- asNamespace("frmtmb")
orig <- get("gp_marginal_outer_gr", envir = ns)
put <- function(nm, val) {
  unlockBinding(nm, ns)
  assign(nm, val, envir = ns)
  lockBinding(nm, ns)
}
use <- function(new) put("gp_marginal_outer_gr",
                         if (new) orig else function(fit, at) NULL)
kmax0 <- get("gp_marginal_kappa_max", envir = ns)
clear_h <- function(fit) {
  rm(list = intersect(".frmtmb_outer_h", ls(fit$obj$env, all.names = TRUE)),
     envir = fit$obj$env)
}
ses <- function(fit) tryCatch(suppressWarnings(
  sqrt(diag(as.matrix(vcov(fit, full = TRUE))))), error = function(e) NA)
hse <- function(H) tryCatch(sqrt(diag(solve(H))), error = function(e) NA)
# the code's own bound, and the true condition of V at the fit's point
kappas <- function(fit) {
  s <- tryCatch(ns$gp_marginal_setup(fit$frame, max_rows = Inf),
                error = function(e) NULL)
  if (is.null(s)) return(c(NA, NA))
  at <- fit$opt$par
  pn <- names(at)
  sig2 <- if (s$has_sigma) exp(2 * at[[which(pn == "betad")]]) else 0
  noise <- (sig2 / s$w + s$se2) / s$nj
  V <- as.matrix(ns$gp_marginal_cov(s, unname(at[pn == "theta"])))
  ev <- eigen(V + diag(noise, s$u), symmetric = TRUE, only.values = TRUE)$values
  c((sum(diag(V)) + max(noise)) / min(noise), max(ev) / min(ev))
}
row_onoff <- function(set, lab, ...) {
  res <- list()
  for (new in c(FALSE, TRUE)) {
    use(new)
    r <- fe_fit(...)
    if (inherits(r$fit, "error")) {
      cat(set, lab, "ERROR", conditionMessage(r$fit), "\n")
      use(TRUE)
      return(NULL)
    }
    se <- ses(r$fit)
    res[[as.character(new)]] <- list(fit = r$fit, se = se, w = sort(substr(
      r$warn, 1, 60)), nan = paste(which(!is.finite(se)), collapse = ","),
      ll = as.numeric(logLik(r$fit)))
  }
  use(TRUE)
  o <- res[["FALSE"]]; n <- res[["TRUE"]]
  route <- tryCatch(is.function(orig(o$fit, o$fit$opt$par)),
                    error = function(e) paste("ERROR", conditionMessage(e)))
  ok <- is.finite(o$se) & is.finite(n$se)
  list(o = o, n = n, route = route,
       se_rel = if (any(ok)) max(abs(n$se[ok] / o$se[ok] - 1)) else NA,
       nan_same = identical(o$nan, n$nan), warn_same = identical(o$w, n$w),
       ll_same = identical(o$ll, n$ll))
}
if (part == "near") {
  out <- file.path(RV_OUT, "se-near.tsv")
  cat("lab\tkap_code\tcond_true\troute\tll_same\tnan_same\twarn_same\tse_rel_onoff\tse_rel_gated\tse_rel_forced\tnan_old\tnan_new\n",
      file = out)
  curves <- list(sin = function(x) sin(x), mix = function(x) 2 * sin(3 * x) + x)
  specs <- list(plain = y ~ gp(x), cov = y ~ z + gp(x),
                grp = y ~ gp(x) + (1 | g), rep = y ~ gp(x))
  i <- 0L
  for (seed in 1:2) for (cv in names(curves)) for (sp in names(specs))
  for (sg in c(0.03, 0.018, 0.012, 0.01, 0.008, 0.006, 0.004)) {
    i <- i + 1L
    set.seed(71000 + 13 * seed + 7 * i)
    u <- 100L
    xs <- sort(round(runif(u, 0, 10), 4))
    x <- if (sp == "rep") rep(xs, each = 2L) else xs
    n <- length(x)
    z <- rnorm(n); g <- factor(sample(1:6, n, TRUE))
    y <- curves[[cv]](x) + 0.3 * z * (sp == "cov") +
      0.5 * rnorm(6)[g] * (sp == "grp") + sg * rnorm(n)
    d <- data.frame(x, z, g, y)
    lab <- sprintf("%s_%s_s%g_seed%d", cv, sp, sg, seed)
    r <- row_onoff("near", lab, specs[[sp]], data = d)
    if (is.null(r)) next
    fit <- r$o$fit
    kk <- kappas(fit)
    # the outer Hessian three ways on one fit
    use(FALSE); clear_h(fit); Hz <- ns$fit_outer_hessian(fit)$H
    use(TRUE); clear_h(fit); Hg <- ns$fit_outer_hessian(fit)$H
    put("gp_marginal_kappa_max", Inf); clear_h(fit)
    Hf <- ns$fit_outer_hessian(fit)$H
    put("gp_marginal_kappa_max", kmax0); clear_h(fit)
    sz <- hse(Hz); sgt <- hse(Hg); sf <- hse(Hf)
    rel <- function(s) {
      ok <- is.finite(s) & is.finite(sz)
      if (any(ok)) max(abs(s[ok] / sz[ok] - 1)) else NA
    }
    cat(sprintf("%s\t%.3e\t%.3e\t%s\t%s\t%s\t%s\t%.3e\t%.3e\t%.3e\t%s\t%s\n",
                lab, kk[1], kk[2], r$route, r$ll_same, r$nan_same, r$warn_same,
                r$se_rel, rel(sgt), rel(sf), r$o$nan, r$n$nan),
        file = out, append = TRUE)
  }
  cat("DONE", out, "\n")
}
if (part == "fall") {
  out <- file.path(RV_OUT, paste0("se-fall", RV_TAG, ".tsv"))
  cat("lab\troute\tll_same\tnan_same\twarn_same\tse_rel\n", file = out)
  set.seed(72000)
  n <- 80L
  x <- sort(round(runif(n, 0, 8), 4)); z <- rnorm(n)
  t <- rep(1:10, 8); gg <- factor(rep(1:8, each = 10))
  f <- sin(x)
  d <- data.frame(x, z, t, gg, w = runif(n, 0.5, 2), s = 0.1,
                  y = f + 0.2 * rnorm(n),
                  k = rpois(n, exp(1 + 0.5 * f)),
                  b = rbinom(n, 1, plogis(f)),
                  ln = exp(0.5 + 0.3 * f + 0.2 * rnorm(n)),
                  yt = f + 0.2 * rt(n, 4))
  d$ym <- d$y; d$ym[c(5, 17, 40)] <- NA
  d$xo <- d$z + rnorm(n, 0, 0.1)
  cases <- list(
    poisson = list(k ~ gp(x), family = poisson()),
    bernoulli = list(b ~ gp(x), family = bernoulli()),
    lognormal = list(ln ~ gp(x), family = lognormal()),
    student = list(yt ~ gp(x), family = student()),
    latent_pois_ar = list(k ~ gp(x) + ar(t, gg, cov = TRUE), family = poisson()),
    gauss_ar = list(y ~ gp(x) + ar(t, gg)),
    gauss_arcov = list(y ~ gp(x) + ar(t, gg, cov = TRUE)),
    sigma_x = list(bf(y ~ gp(x), sigma ~ z)),
    smooth = list(y ~ gp(x) + s(z)),
    mi = list(bf(ym | mi() ~ gp(x))),
    me = list(y ~ me(xo, s) + gp(x)),
    weights = list(y | weights(w) ~ gp(x)),
    se = list(y | se(s, sigma = TRUE) ~ gp(x)),
    se_only = list(y | se(s) ~ gp(x)),
    prior = list(y ~ gp(x), prior = set_prior("normal(0, 1)", class = "sdgp")),
    reml = list(y ~ z + gp(x), REML = TRUE),
    approx = list(y ~ gp(x, k = 15)),
    by_fac = list(y ~ gp(x, by = gg)),
    two_gp = list(y ~ gp(x) + gp(z)),
    nl_body = list(bf(y ~ a + b * z, a ~ gp(x), b ~ 1, nl = TRUE)),
    sigma_const = list(bf(y ~ gp(x), sigma = 0.2)))
  for (nm in names(cases)) {
    a1 <- cases[[nm]]
    r <- tryCatch(do.call(row_onoff, c(list("fall", nm), a1[1],
                                       list(data = d), a1[-1])),
                  error = function(e) {
                    cat("fall", nm, "ERROR", conditionMessage(e), "\n")
                    NULL
                  })
    if (is.null(r)) next
    cat(sprintf("%s\t%s\t%s\t%s\t%s\t%.3e\n", nm, r$route, r$ll_same,
                r$nan_same, r$warn_same, r$se_rel), file = out, append = TRUE)
  }
  cat("DONE", out, "\n")
}
if (part == "time") {
  for (u in c(100L, 500L)) {
    set.seed(20261103)
    x <- sort(runif(u, 0, 10))
    d <- data.frame(x = x, y = sin(x) + 0.3 * rnorm(u))
    tt <- list(on = numeric(0), off = numeric(0))
    hh <- list(on = numeric(0), off = numeric(0))
    reps <- if (u == 100L) 6L else 2L
    for (i in seq_len(reps)) {
      for (arm in if (i %% 2L) c("on", "off") else c("off", "on")) {
        use(arm == "on")
        t0 <- proc.time()[[3]]
        fit <- suppressWarnings(frm(y ~ gp(x), data = d))
        v <- vcov(fit, full = TRUE)
        tt[[arm]] <- c(tt[[arm]], proc.time()[[3]] - t0)
        # the outer Hessian alone, cache cleared, on this fit
        clear_h(fit)
        t0 <- proc.time()[[3]]
        h <- ns$fit_outer_hessian(fit)
        hh[[arm]] <- c(hh[[arm]], proc.time()[[3]] - t0)
        cat(sprintf("TIME pos %d rep %d arm %s frm+vcov %.2f s hessian %.3f s ll %.6f\n",
                    u, i, arm, tail(tt[[arm]], 1), tail(hh[[arm]], 1),
                    as.numeric(logLik(fit))))
      }
    }
    use(TRUE)
    cat(sprintf("TIMESUM pos %d frm+vcov min on %.2f off %.2f ratio %.2f | hessian min on %.3f off %.3f ratio %.2f\n",
                u, min(tt$on), min(tt$off), min(tt$off) / min(tt$on),
                min(hh$on), min(hh$off), min(hh$off) / min(hh$on)))
  }
  cat("DONE time\n")
}
