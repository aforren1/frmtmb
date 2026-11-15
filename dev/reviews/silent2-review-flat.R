# Review of lane silent2, item 5207: is a threshold flagged by the SE cut
# really not placed by the data? For each fit, per component, the lane's
# latent threshold SE (mixture_ord_thres_se()); for each flagged component
# (and any above 20), the worst threshold is moved by D latent units along
# the quadratic profile direction (V J' / J V J'), unprofiled, and then
# re-optimized with the threshold held (penalty 1e4 (tau - target)^2):
# the profile drop in log-likelihood. A placed threshold loses > 1.92 at
# D = 20 one way or the other; a flat one loses ~0 both ways.
#   Rscript dev/reviews/silent2-review-flat.R <set> <out.txt>
# set: gained (fits of the lane's 21 gained), new (cases the lane did not
# run: weights, disc ~ x, small n, probit, rare category), thres4.
args <- commandArgs(TRUE)
set <- args[1]
.libPaths(c("C:/Users/adf44/source/r/wt-silent2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
`%||%` <- function(a, b) if (is.null(a)) b else a
mk <- function(seed, n = 300, mix = FALSE) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n))
  lat <- if (mix) {
    cls <- rbinom(n, 1, 0.4)
    ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) + rlogis(n)
  } else 0.8 * d$x + rlogis(n)
  d$y <- 1L + (lat > -1 + 0.1 * d$z) + (lat > 0) + (lat > 1 - 0.1 * d$z)
  d
}
gen <- function(seed, n = 300, K = 4, b1 = 1.5, b2 = -0.8, c1 = 1, c2 = -1,
                p = 0.4, noise = "logit", probs = NULL) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n))
  cls <- rbinom(n, 1, p)
  e <- switch(noise, logit = rlogis(n), probit = rnorm(n),
              cloglog = log(rexp(n)))
  lat <- ifelse(cls == 1, b1 * d$x + c1, b2 * d$x + c2) + e
  cuts <- stats::quantile(lat, probs %||% (seq_len(K - 1) / K))
  d$y <- 1L + rowSums(outer(lat, cuts, ">"))
  d
}
omx_data <- function(seed, n = 400) {
  set.seed(seed)
  x <- rnorm(n); z <- rnorm(n)
  g <- factor(sample(c("a", "b"), n, TRUE))
  cls <- rbinom(n, 1, 0.4)
  lat <- ifelse(cls == 1, 1.5 * x + 1, -0.8 * x - 1) + rlogis(n)
  y <- as.integer(cut(lat, c(-Inf, -1.5, 0, 1.5, Inf)))
  data.frame(y, x, z, g)
}
cases <- list()
add <- function(name, data, formula, family) {
  cases[[name]] <<- list(data = data, formula = formula, family = family)
}
CC <- quote(mixture(cumulative(), cumulative()))
CS <- quote(mixture(cumulative(), sratio()))
if (set == "gained") {
  add("cum_cum_seed16", quote(mk(16, mix = TRUE)), bf(y ~ x + z), CC)
  add("cum_sratio_seed5", quote(mk(5, mix = TRUE)), bf(y ~ x + z), CS)
  add("cc_logit_K7_seed6", quote(gen(6, K = 7)), bf(y ~ x + z), CC)
  add("dpar_formula_seed2", quote(gen(2, K = 4)),
      bf(y ~ 1, mu1 ~ x, mu2 ~ x + z), CS)
  add("steep_n300_s3_K6_seed2",
      quote(gen(2, n = 300, K = 6, b1 = 3, b2 = -3, c1 = 2, c2 = -2)),
      bf(y ~ x + z), CC)
  add("cs_sr_acat_seed1", quote(gen(1, K = 4)), bf(y ~ cs(x) + z),
      quote(mixture(sratio(), acat())))
} else if (set == "thres4") {
  add("thres4_test", quote(omx_data(20261016)), bf(y | thres(4) ~ x), CS)
} else if (set == "new") {
  for (s in 1:3) {
    add(sprintf("w_const5_seed%d", s),
        substitute(transform(mk(S, mix = TRUE), w = 5), list(S = s)),
        bf(y | weights(w) ~ x + z), CC)
    add(sprintf("w_none_seed%d", s), substitute(mk(S, mix = TRUE), list(S = s)),
        bf(y ~ x + z), CC)
    add(sprintf("w_runif_seed%d", s),
        substitute({dd <- mk(S, mix = TRUE); set.seed(100 + S)
                    dd$w <- runif(nrow(dd), 0.2, 1.8); dd}, list(S = s)),
        bf(y | weights(w) ~ x + z), CC)
    add(sprintf("disc2_x_seed%d", s), substitute(gen(S, K = 5), list(S = s)),
        bf(y ~ x + z, disc2 ~ x), CC)
  }
  for (s in 1:4) for (n in c(80, 150)) {
    add(sprintf("small_n%d_cc_seed%d", n, s),
        substitute(gen(S, n = N, K = 4), list(S = s, N = n)), bf(y ~ x), CC)
    add(sprintf("small_n%d_probit_seed%d", n, s),
        substitute(gen(S, n = N, K = 4, noise = "probit"), list(S = s, N = n)),
        bf(y ~ x), quote(mixture(cumulative("probit"), cumulative("probit"))))
  }
  for (s in 1:3) {
    add(sprintf("rare_cat_seed%d", s),
        substitute(gen(S, n = 300, K = 5, probs = c(0.03, 0.3, 0.6, 0.97)),
                   list(S = s)), bf(y ~ x + z), CC)
    add(sprintf("separated_n2000_seed%d", s),
        substitute(gen(S, n = 2000, K = 4, b1 = 2, b2 = -2, c1 = 2, c2 = -2),
                   list(S = s)), bf(y ~ x + z), CC)
  }
}

thres_fun <- function(f, k, resp = "y") {
  lp <- f$frame[["linpreds"]][[ns$linpred_key(resp, paste0("mu", k))]]
  ob <- ns$ord_lp_block(f$frame, f$spec, lp)
  p <- f$opt$par
  pos <- which(names(p) == ob$comp)
  list(pos = pos, fam = ob$fam,
       tau = function(pp) as.numeric(ns$ord_threshold_values(ob$fam, pp[pos])),
       unident = ob$fam[["thres"]][["unident"]] %||% integer(0))
}
flat_check <- function(f, k, resp = "y") {
  tf <- thres_fun(f, k, resp)
  p <- f$opt$par
  V <- f$cache$sdr$cov.fixed %||% ns$se_fit_analysis(f)$V
  nll0 <- f$obj$fn(p)
  t0 <- tf$tau(p)
  J <- vapply(tf$pos, function(i) {
    h <- 1e-6 * max(1, abs(p[i])); up <- p; dn <- p
    up[i] <- up[i] + h; dn[i] <- dn[i] - h
    (tf$tau(up) - tf$tau(dn)) / (2 * h)
  }, numeric(length(t0)))
  J <- matrix(J, nrow = length(t0))
  Vb <- V[tf$pos, tf$pos, drop = FALSE]
  se <- sqrt(pmax(rowSums((J %*% Vb) * J), 0))
  se[tf$unident] <- NA
  j <- which.max(se)
  g <- rep(0, length(p)); g[tf$pos] <- J[j, ]
  Vg <- as.vector(V %*% g)
  dvec <- Vg / sum(g * Vg)
  out <- character()
  for (D in c(-20, -5, 5, 20)) {
    pt <- p + D * dvec
    un <- f$obj$fn(pt) - nll0
    target <- t0[j] + D
    obj <- function(pp) {
      v <- f$obj$fn(pp)
      if (!is.finite(v)) return(1e10)
      v + 1e4 * (tf$tau(pp)[j] - target)^2
    }
    o <- tryCatch(stats::nlminb(if (is.finite(un)) pt else p, obj,
                                control = list(iter.max = 400, eval.max = 800)),
                  error = function(e) NULL)
    prof <- if (is.null(o)) NA else f$obj$fn(o$par) - nll0
    reached <- if (is.null(o)) NA else tf$tau(o$par)[j] - t0[j]
    out <- c(out, sprintf("D=%+g unprof %.3g prof %.3g (moved %.3g)", D, un,
                          prof, reached))
  }
  f$obj$fn(p)
  sprintf("comp %d thres %d at %.3g se(natural) %.3g | %s", k, j, t0[j], se[j],
          paste(out, collapse = " ; "))
}

for (nm in names(cases)) {
  cs <- cases[[nm]]
  d <- eval(cs$data)
  w <- character()
  f <- tryCatch(withCallingHandlers(
    frm(cs$formula, family = eval(cs$family), data = d),
    warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") },
    message = function(m) invokeRestart("muffleMessage")),
    error = function(e) e)
  if (inherits(f, "error")) { cat("CASE", nm, "ERROR", conditionMessage(f), "\n"); next }
  se <- ns$mixture_ord_thres_se(f, "y")
  degen <- grepl("degenerate boundary", w)
  unpl <- grepl("a threshold the data do not place", w)
  other <- any(degen) && any(grepl("step function|run off|same number", w[degen]))
  cat(sprintf("CASE %s ll %.4f code %d latentSE %s degen %s by_se %s other_reason %s nwarn %d\n",
              nm, as.numeric(logLik(f)), f$opt$convergence,
              paste(signif(se, 4), collapse = "/"), any(degen), any(unpl), other,
              length(w)))
  for (k in seq_along(se)) {
    if (is.finite(se[k]) && se[k] > 20) {
      cat("   ", tryCatch(flat_check(f, k), error = function(e) paste("flat_check error", conditionMessage(e))), "\n")
    }
  }
}
cat("DONE\n")
