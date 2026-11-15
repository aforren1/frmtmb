# Review of lane predraw2: importance-ratio correctness, RNG, cache, cost.
# Usage: Rscript dev/reviews/predraw2-review-a.R [part ...]
local({
.libPaths(c("C:/Users/adf44/source/r/wt-predraw2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("lib:", dirname(find.package("frmtmb")), "\n")
parts <- commandArgs(trailingOnly = TRUE)
if (!length(parts)) parts <- c("jac", "pois", "rng", "rds", "ord", "cost")
q <- function(expr) suppressWarnings(suppressMessages(expr))

mix_fit <- function() {
  set.seed(101)
  x <- stats::rnorm(300)
  invisible(stats::runif(300))
  d <- data.frame(x = x, y = stats::rpois(300, exp(1 + 0.3 * x)))
  q(frm(bf(y ~ x), family = mixture(poisson(), poisson()), data = d))
}

# ---- 1. the log Jacobian anchor, on a 1-D target with a known answer ----
# target: -log density f(t) = t^2/2 for t >= 0 and (t/a)^2/2 for t < 0,
# so P(t < 0) = a / (1 + a). The side law for the lower side, built as
# draw_edge_side() builds it, has knots z * a, and maps x -> a x there:
# the proposal is exactly the target's shape on each side with mass 1/2.
if ("jac" %in% parts) {
  cat("\n== jac: 1-D target, lower side narrowed\n")
  for (a in c(0.3, 0.1)) {
    env <- new.env()
    obj <- list(fn = function(p) {
      t <- p[1]
      if (t >= 0) t^2 / 2 else (t / a)^2 / 2
    }, env = env)
    object <- list(obj = obj, opt = list(par = 0))
    v <- a * (1:4)
    law <- list(lo = list(v), hi = list(1), s = 1, Lc = NULL,
                edge = integer(0), narrowed = 1L)
    set.seed(1)
    N <- 20000
    Z <- matrix(stats::rnorm(N), 1, N)
    L <- matrix(1, 1, 1)
    D <- 0 + ns$draw_edge_apply(law, L, Z)
    dw <- list(D = D, Z = Z, L = L, law = law, v0 = 0, held = integer(0),
               off = 0L)
    r <- ns$draw_logratio(object, dw, D, Z)
    set.seed(2)
    rs <- ns$draw_resample(r$lw, r$stratum)
    lane <- mean(D[1, rs$idx] < 0)
    # the proposal's density carries the map's slope itself, not its
    # slope over the slope at 0 on the side the draw fell on
    x <- Z[1, ]
    lj_abs <- ifelse(x < 0, log(ns$draw_edge_slope(abs(x), v)), 0)
    lw_abs <- r$lw + ifelse(x < 0, log(ns$draw_edge_slope(0, v)), 0)
    set.seed(2)
    rs2 <- ns$draw_resample(lw_abs, r$stratum)
    fixed <- mean(D[1, rs2$idx] < 0)
    cat(sprintf("a %.2f exact P(t<0) %.4f | proposal %.4f | lane reweighted %.4f (k %.2f) | absolute Jacobian %.4f (k %.2f) | median e over 20 %.4f\n",
                a, a / (1 + a), mean(D[1, ] < 0), lane, rs$k, fixed, rs2$k,
                stats::median(r$e[1:20])))
  }
}

# ---- 2. exact posterior: poisson intercept-only, flat prior on log rate --
if ("pois" %in% parts) {
  cat("\n== pois: reweighting forced on, exact Gamma(S, n) posterior\n")
  set.seed(9)
  d <- data.frame(y = c(rep(0L, 47), 1L, 1L, 1L))
  fit <- q(frm(bf(y ~ 1), family = poisson(), data = d))
  S <- sum(d$y); n <- nrow(d)
  real <- ns$draw_check
  for (force in c(FALSE, TRUE)) {
    if (force) {
      assignInNamespace("draw_check", function(object, dw) {
        list(fire = TRUE, excess = 9, n = 20L)
      }, ns = "frmtmb")
    }
    set.seed(1)
    dr <- q(ns$predict_par_drawer(fit, TRUE, 4000L))
    eta <- vapply(1:4000, function(s) attr(dr, "delta")(s), 0) +
      fit$opt$par[1]
    lam <- exp(eta)
    pr <- c(0.025, 0.1, 0.5, 0.9, 0.975)
    cat(sprintf("%s quantiles of rate: %s | exact %s | mean %.4f exact %.4f\n",
                if (force) "reweighted" else "normal    ",
                paste(format(stats::quantile(lam, pr), digits = 4), collapse = " "),
                paste(format(stats::qgamma(pr, S, n), digits = 4), collapse = " "),
                mean(lam), S / n))
  }
  assignInNamespace("draw_check", real, ns = "frmtmb")
}

# ---- 3. RNG: first call (check computed) against second (cached) --------
if ("rng" %in% parts) {
  cat("\n== rng\n")
  kinds <- list(default = NULL,
                lecuyer_bm_round = list("L'Ecuyer-CMRG", "Box-Muller",
                                        "Rounding"),
                mt_ak = list("Mersenne-Twister", "Ahrens-Dieter",
                             "Rejection"))
  for (nm in names(kinds)) {
    RNGkind("default", "default", "default")
    fit <- mix_fit()
    if (!is.null(kinds[[nm]])) {
      suppressWarnings(do.call(RNGkind, kinds[[nm]]))
    } else {
      RNGkind("default", "default", "default")
    }
    k0 <- RNGkind()
    res <- list()
    for (i in 1:2) {
      suppressWarnings(set.seed(1))
      p <- q(predict(fit, ndraws = 200, summary = FALSE))
      res[[i]] <- list(p = p, seed = .Random.seed, kind = RNGkind(),
                       after = stats::rnorm(3))
    }
    cat(sprintf("%-18s fire %s | draws identical %s | .Random.seed identical %s | kind kept %s %s | next rnorm identical %s (%s vs %s)\n",
                nm, format(fit$cache[["draw_check"]]$fire),
                identical(res[[1]]$p, res[[2]]$p),
                identical(res[[1]]$seed, res[[2]]$seed),
                identical(res[[1]]$kind, k0), identical(res[[2]]$kind, k0),
                identical(res[[1]]$after, res[[2]]$after),
                format(res[[1]]$after[1], digits = 8),
                format(res[[2]]$after[1], digits = 8)))
  }
  RNGkind("default", "default", "default")
  # no .Random.seed before the call
  fit <- mix_fit()
  if (exists(".Random.seed", globalenv())) rm(".Random.seed", envir = globalenv())
  invisible(q(predict(fit, ndraws = 20)))
  cat("after a call with no seed: .Random.seed exists",
      exists(".Random.seed", globalenv()), "kind", RNGkind(), "\n")
}

# ---- 4. cache through saveRDS / readRDS -----------------------------------
if ("rds" %in% parts) {
  cat("\n== rds\n")
  fit <- mix_fit()
  set.seed(1)
  p1 <- q(predict(fit, ndraws = 200, summary = FALSE))
  f <- tempfile(fileext = ".rds")
  saveRDS(fit, f)
  g <- readRDS(f)
  cat("verdict carried:", !is.null(g$cache[["draw_check"]]), "\n")
  said <- character(0)
  set.seed(1)
  p2 <- withCallingHandlers(predict(g, ndraws = 200, summary = FALSE),
    message = function(m) {
      said <<- c(said, class(m)[1])
      invokeRestart("muffleMessage")
    }, warning = function(w) {
      said <<- c(said, class(w)[1])
      invokeRestart("muffleWarning")
    })
  cat("read fit: draws identical to before save", identical(p1, p2),
      "| said:", said, "\n")
  # a fit read back with no verdict
  saveRDS(mix_fit(), f)
  h <- readRDS(f)
  set.seed(1)
  p3 <- q(predict(h, ndraws = 200, summary = FALSE))
  cat("fresh read fit: check fire", h$cache[["draw_check"]]$fire,
      "excess", h$cache[["draw_check"]]$excess,
      "| draws identical to the live fit's", identical(p1, p3), "\n")
}

# ---- 5. the 3277 ordinal mixture: narrowed sides and their slopes at 0 ----
if ("ord" %in% parts) {
  cat("\n== ord (3277 test fit)\n")
  set.seed(7)
  d <- data.frame(x = stats::rnorm(400))
  k <- stats::rbinom(400, 1, 0.3)
  d$yo <- cut(d$x + ifelse(k == 1, 2, 0) + stats::rlogis(400),
              c(-Inf, -1, 0.5, 2, Inf), labels = FALSE)
  fit <- q(frm(bf(yo ~ x), family = mixture(cumulative(), cumulative()),
               data = d))
  set.seed(13)
  dr <- q(ns$predict_par_drawer(fit, TRUE, 4000L))
  law <- fit$cache[["draw_edge_law"]]
  chk <- fit$cache[["draw_check"]]
  nms <- ns$outer_par_names(fit)
  cat("check fire", chk$fire, "excess", chk$excess, "n", chk$n, "\n")
  cat("edge:", nms[law$edge], "| narrowed:", nms[law$narrowed], "\n")
  one <- function(v) length(v) == 1L && v == 1
  for (j in which(!vapply(law$lo, one, TRUE) | !vapply(law$hi, one, TRUE))) {
    cat(sprintf("  %-40s lo %s slope0 %.3f | hi %s slope0 %.3f\n", nms[j],
                paste(format(law$lo[[j]], digits = 3), collapse = ","),
                ns$draw_edge_slope(0, law$lo[[j]]),
                paste(format(law$hi[[j]], digits = 3), collapse = ","),
                ns$draw_edge_slope(0, law$hi[[j]])))
  }
  # the reweighted draws under the lane's ratio and under the proposal's
  # own density (absolute log Jacobian), on the same proposal draws
  ds <- ns$fit_draw_space(fit)
  v0 <- ns$fit_outer_vector(fit, ds$map)
  set.seed(13)
  Z <- matrix(stats::rnorm(length(v0) * 4000), length(v0), 4000)
  L <- chol(ds$V)
  D <- v0 + ns$draw_edge_apply(law, L, Z)
  dw <- list(D = D, Z = Z, L = L, law = law, v0 = v0, held = integer(0),
             off = length(v0) - length(fit$opt$par))
  r <- ns$draw_logratio(fit, dw, D, Z)
  X <- crossprod(law$Lc %||% L, Z)
  add <- numeric(4000)
  for (j in setdiff(which(!vapply(law$lo, one, TRUE) |
                          !vapply(law$hi, one, TRUE)), law$edge)) {
    neg <- X[j, ] < 0
    add[neg] <- add[neg] + log(ns$draw_edge_slope(0, law$lo[[j]]))
    add[!neg] <- add[!neg] + log(ns$draw_edge_slope(0, law$hi[[j]]))
  }
  set.seed(3)
  a1 <- ns$draw_resample(r$lw, r$stratum)
  set.seed(3)
  a2 <- ns$draw_resample(r$lw + add, r$stratum)
  cat("log ratio added by the anchor: range", format(range(add), digits = 3),
      "| k lane", format(a1$k, digits = 3), "k absolute", format(a2$k, digits = 3), "\n")
  for (j in law$narrowed) {
    cat(sprintf("  %-40s share below estimate: proposal %.3f lane %.3f absolute %.3f\n",
                nms[j], mean(X[j, ] < 0), mean(X[j, a1$idx] < 0),
                mean(X[j, a2$idx] < 0)))
  }
}

# ---- 6. cost of the check on a large GLMM --------------------------------
if ("cost" %in% parts) {
  cat("\n== cost: poisson, 10000 rows, (1 | g) with 2000 levels\n")
  set.seed(21)
  n <- 10000; G <- 2000
  d <- data.frame(x = stats::rnorm(n), g = factor(sample.int(G, n, TRUE)))
  d$y <- stats::rpois(n, exp(0.2 + 0.3 * d$x + stats::rnorm(G, 0, 0.5)[d$g]))
  tf <- system.time(fit <- q(frm(bf(y ~ x + (1 | g)), family = poisson(),
                                 data = d)))
  cat("fit", tf[["elapsed"]], "s\n")
  t1 <- system.time(q(predict(fit, ndraws = 200)))[["elapsed"]]
  chk <- fit$cache[["draw_check"]]
  t2 <- system.time(q(predict(fit, ndraws = 200)))[["elapsed"]]
  rm("draw_check", envir = fit$cache)
  t3 <- system.time(q(predict(fit, ndraws = 200)))[["elapsed"]]
  rm("draw_check", envir = fit$cache)
  law <- fit$cache[["draw_edge_law"]]
  ds <- ns$fit_draw_space(fit)
  v0 <- ns$fit_outer_vector(fit, ds$map)
  dw <- list(L = chol(ds$V), law = if (isFALSE(law)) NULL else law, v0 = v0,
             held = integer(0), off = length(v0) - length(fit$opt$par))
  t4 <- system.time(ck <- ns$draw_check(fit, dw))[["elapsed"]]
  cat(sprintf("first predict %.2f s | second (all cached) %.2f s | law cached, check recomputed %.2f s | draw_check alone %.2f s | excess %.4f n %d\n",
              t1, t2, t3, t4, chk$excess, chk$n))
  t5 <- system.time(q(predict(fit, ndraws = 1000)))[["elapsed"]]
  cat(sprintf("predict 1000 draws (cached) %.2f s\n", t5))
}
})
