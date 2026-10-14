# Reviewer of lane illse, attack 1-3: a degree-7 raw polynomial with
# (1 | g) under gaussian. Exact covariance at the fit's own point: the
# closed-form marginal likelihood (compound symmetry per group) in 256
# bits, Hessian by 256-bit second differences, inverted in 256 bits.
# Every consumer against the exact answer.
#   REV_ARM=lane|base Rscript dev/illse-rev-lmm.R <seeds> [deg]
source("C:/Users/adf44/source/r/frmtmb-wt-illse/dev/illse-rev-mp.R")
a <- commandArgs(TRUE)
seeds <- if (length(a)) eval(parse(text = a[1])) else 101:103
deg <- if (length(a) > 1) as.integer(a[2]) else 7L
xn <- paste0("x", seq_len(deg))
fo <- stats::as.formula(paste("y ~", paste(xn, collapse = " + "),
                              "+ (1 | g)"))
fo0 <- stats::as.formula(paste("~", paste(xn, collapse = " + ")))
mk <- function(seed) {
  set.seed(seed)
  d <- data.frame(x = stats::runif(200, 1, 2),
                  g = factor(rep(1:20, length.out = 200)))
  d$y <- sin(3 * d$x) + stats::rnorm(20, 0, 0.5)[d$g] +
    stats::rnorm(200, 0, 0.3)
  for (k in seq_len(deg)) d[[xn[k]]] <- d$x^k
  d
}
nd_of <- function(d) {
  nd <- data.frame(x = c(1, 1.25, 1.5, 1.75, 2), g = d$g[1])
  for (k in seq_len(deg)) nd[[xn[k]]] <- nd$x^k
  nd
}
for (seed in seeds) {
  d <- mk(seed)
  r <- cap(frm(fo, data = d))
  if (inherits(r$v, "rev_err")) { cat("seed", seed, "ERROR", r$v, "\n"); next }
  fit <- r$v
  p <- fit$opt$par
  pn <- names(p)
  pb <- which(pn == "beta"); pu <- which(pn == "betad")
  pw <- which(pn == "theta")
  X <- stats::model.matrix(fo0, d)
  stopifnot(max(abs(X - as.matrix(fit$frame$linpreds[[1]]$X))) == 0)
  # sufficient statistics, exact in 256 bits
  Xm <- mM(X); ym <- m(d$y)
  XX <- crossprod(Xm); Xy <- crossprod(Xm, ym); yy <- sum(ym^2)
  gi <- split(seq_len(nrow(d)), d$g)
  S <- lapply(gi, function(ix) list(s = colSums(Xm[ix, , drop = FALSE]),
                                     sy = sum(ym[ix]), n = length(ix)))
  nll <- function(th) {
    b <- th[pb]; a2 <- exp(2 * th[pu]); c2 <- exp(2 * th[pw])
    bm <- mpfr2array(b, dim = c(length(pb), 1L))
    rss <- yy - 2 * sum(Xy * b) + sum(bm * (XX %*% bm))
    out <- rss / a2
    for (st in S) {
      sr <- st$sy - sum(st$s * b)
      out <- out + (st$n - 1) * log(a2) + log(a2 + st$n * c2) -
        c2 * sr^2 / (a2 * (a2 + st$n * c2))
    }
    out / 2
  }
  # the objective agrees with frmtmb's up to a constant
  c0 <- as.numeric(fit$obj$fn(p)) - num(nll(m(p)))
  pp <- p; pp[pb[deg + 1]] <- pp[pb[deg + 1]] + 1e-4; pp[pu] <- pp[pu] + 0.01
  c1 <- as.numeric(fit$obj$fn(pp)) - num(nll(m(pp)))
  fit$obj$fn(p)
  H <- mp_hessian(nll, p)
  Ve <- mp_inv(H)
  g_ex <- mp_grad(nll, p)
  ge <- num(g_ex)
  newton_gain <- num(sum(g_ex * (Ve %*% mpfr2array(g_ex, dim = c(length(p), 1L))))) / 2
  V <- cap(vcov(fit, full = TRUE))
  e <- cov_err(V$v, Ve)
  cat(sprintf(paste0("seed %d deg %d | code %d | fit warn %s | lost %d | ",
                     "obj const diff %.2e | Newton gain at fit %.2e | ",
                     "root %s | vcov warn %s\n"),
              seed, deg, fit$opt$convergence, wtag(r$w),
              length(ns$sdr_of(fit)$se_lost), abs(c0 - c1), newton_gain,
              !is.null(ns$sdr_of(fit)$cov_fixed_root), wtag(V$w)))
  cat(sprintf("  vcov(full) every entry: cor %.2e | rel %.2e | diag %.2e\n",
              e[["cor"]], e[["rel"]], e[["diag"]]))
  se_ex <- sqrt(num(diag(Ve)))
  nd <- nd_of(d)
  A <- stats::model.matrix(fo0, nd)
  v_pop <- mp_rowvar(A, Ve[pb, pb])
  # a group's prediction: x0'beta + bhat_g(theta), bhat the conditional
  # mode; variance = conditional + propagated (TMB's construction)
  g1 <- as.character(nd$g[1])
  st1 <- S[[g1]]
  bhat <- function(th) {
    b <- th[pb]; a2 <- exp(2 * th[pu]); c2 <- exp(2 * th[pw])
    c2 * (st1$sy - sum(st1$s * b)) / (a2 + st1$n * c2)
  }
  a2 <- exp(2 * m(p[pu])); c2 <- exp(2 * m(p[pw]))
  vcond <- 1 / (st1$n / a2 + 1 / c2)
  v_grp <- vapply(seq_len(nrow(A)), function(i) {
    fi <- function(th) sum(m(A[i, ]) * th[pb]) + bhat(th)
    gg <- mp_grad(fi, p)
    num(vcond + sum(gg * (Ve %*% mpfr2array(gg, dim = c(length(p), 1L)))))
  }, 0)
  out <- list()
  # summary()
  s <- cap(summary(fit))
  if (!inherits(s$v, "rev_err")) {
    fx <- s$v$fixed
    out$summary_fixed <- fx[, "Est.Error"] / se_ex[pb]
    sp <- s$v$spec_pars
    if (!is.null(sp)) out$summary_sigma <- sp["sigma", "Est.Error"] /
      (exp(p[pu]) * se_ex[pu])
    rg <- s$v$random$g
    if (!is.null(rg)) out$summary_sd <- rg[1, "Est.Error"] /
      (exp(p[pw]) * se_ex[pw])
  }
  out$summary_w <- wtag(s$w)
  # confint() Wald
  ci <- cap(confint(fit))
  if (!inherits(ci$v, "rev_err")) {
    hw <- (ci$v[, 2] - ci$v[, 1]) / 2
    out$confint_beta <- hw[seq_along(pb)] / (stats::qnorm(0.975) * se_ex[pb])
  }
  out$confint_w <- wtag(ci$w)
  # frm_linpred(se.fit = TRUE) and fitted()
  lp <- cap(frm_linpred(fit, newdata = nd, re_formula = NA, se.fit = TRUE))
  if (!inherits(lp$v, "rev_err")) out$linpred_pop <- lp$v$se.fit / sqrt(v_pop)
  out$linpred_w <- wtag(lp$w)
  fp <- cap(fitted(fit, newdata = nd, re_formula = NA))
  if (!inherits(fp$v, "rev_err")) out$fitted_pop <- fp$v[, "Est.Error"] / sqrt(v_pop)
  out$fitted_pop_w <- wtag(fp$w)
  fg <- cap(fitted(fit, newdata = nd))
  if (!inherits(fg$v, "rev_err")) out$fitted_grp <- fg$v[, "Est.Error"] / sqrt(v_grp)
  out$fitted_grp_w <- wtag(fg$w)
  # hypothesis(): the polynomial at x = 1.5
  hs <- paste0("Intercept + ", paste0(1.5^(seq_len(deg)), " * ", xn,
                                      collapse = " + "), " = 0")
  h <- cap(hypothesis(fit, hs, class = "b"))
  if (!inherits(h$v, "rev_err")) out$hypothesis <- h$v$hypothesis$Est.Error / sqrt(v_pop[3])
  out$hypothesis_w <- wtag(h$w)
  # emmeans
  em <- cap(summary(emmeans::emmeans(fit, fo0,
                                     at = as.list(nd[3, xn]))))
  if (!inherits(em$v, "rev_err")) out$emmeans <- em$v$SE / sqrt(v_pop[3])
  out$emmeans_w <- wtag(em$w)
  # marginaleffects
  mf <- cap(marginaleffects::predictions(fit, newdata = nd, re.form = NA))
  if (!inherits(mf$v, "rev_err")) out$marginaleffects <- mf$v$std.error / sqrt(v_pop)
  out$marginaleffects_w <- wtag(mf$w)
  # parameter draws: predict()'s drawer (internal), 20000 draws, Monte
  # Carlo error of an SD near 0.5 percent
  set.seed(7)
  dr <- cap(ns$predict_par_drawer(fit, TRUE, 20000L))
  if (!inherits(dr$v, "rev_err") && is.function(attr(dr$v, "delta"))) {
    dl <- vapply(seq_len(20000L), attr(dr$v, "delta"), numeric(length(p)))
    out$draws_sd <- apply(A %*% dl[pb, , drop = FALSE], 1L, stats::sd) /
      sqrt(v_pop)
  } else out$draws_err <- if (inherits(dr$v, "rev_err")) dr$v else "no delta"
  out$draws_w <- wtag(dr$w)
  # frm_multiple(): three identical copies, so the pooled SE is the
  # within-imputation one, which is the single fit's
  mm <- cap(frm_multiple(fo, data = list(d, d, d)))
  if (!inherits(mm$v, "rev_err")) {
    hm <- cap(hypothesis(mm$v, hs, class = "b"))
    if (!inherits(hm$v, "rev_err")) {
      out$multiple_hyp <- hm$v$hypothesis$Est.Error / sqrt(v_pop[3])
    }
    out$multiple_hyp_w <- wtag(hm$w)
    sm <- cap(summary(mm$v))
    if (!inherits(sm$v, "rev_err")) {
      out$multiple_fixed <- sm$v$fixed[, "Est.Error"] / se_ex[pb]
    }
    # D1 Wald test of the two top powers
    tst <- xn[deg - 1:0]
    fcol <- "F"
    an <- cap(anova(mm$v, method = "D1", constraint = tst))
    q <- p[pb][deg + 1 - 1:0]
    q <- p[pb[match(tst, colnames(X))]]
    Ub <- Ve[pb[match(tst, colnames(X))], pb[match(tst, colnames(X))]]
    Fe <- num(sum(mpfr2array(m(q), dim = c(2L, 1L)) *
                    (mp_inv(Ub) %*% mpfr2array(m(q), dim = c(2L, 1L))))) / 2
    if (!inherits(an$v, "rev_err")) {
      tab <- as.data.frame(an$v)
      fcol <- intersect(c("F", "F.value", "statistic", "Fstat", "D1"),
                        names(tab))
      out$D1_F <- if (length(fcol)) suppressWarnings(as.numeric(tab[[fcol[1]]])) /
        Fe else NA
      if (!length(fcol)) print(tab)
    } else out$D1_F <- an$v
    out$D1_w <- wtag(an$w)
    # D1 over every power: Q' U^-1 Q with U the whole block
    pk <- pb[match(xn, colnames(X))]
    an2 <- cap(anova(mm$v, method = "D1", constraint = xn))
    qq <- mpfr2array(m(p[pk]), dim = c(deg, 1L))
    Fe2 <- num(sum(qq * (mp_inv(Ve[pk, pk]) %*% qq))) / deg
    if (!inherits(an2$v, "rev_err")) {
      tab2 <- as.data.frame(an2$v)
      out$D1all_F <- if (fcol[1] %in% names(tab2)) suppressWarnings(as.numeric(tab2[[fcol[1]]])) / Fe2 else NA
    } else out$D1all_F <- an2$v
    out$D1all_w <- wtag(an2$w)
  } else out$multiple <- mm$v
  # vcov_cluster(): bread whole; exact sandwich with the same scores
  cl <- ns$resolve_cluster(fit, ~ g)
  Sc <- ns$cluster_scores_at(fit, cl)
  vc <- cap(vcov_cluster(fit, ~ g, full = TRUE))
  if (!inherits(vc$v, "rev_err")) {
    adj <- ns$cr_adjust("CR0", nrow(Sc), nrow(d), length(p))
    SS <- crossprod(mM(Sc))
    Vce <- adj * (Ve %*% SS %*% Ve)
    ec <- cov_err(vc$v, Vce)
    out$cluster_cor <- ec[["cor"]]
    out$cluster_diag <- ec[["diag"]]
    vcl_pop <- mp_rowvar(A, Vce[pb, pb])
    hc <- cap(hypothesis(fit, hs, class = "b", vcov = vc$v))
    if (!inherits(hc$v, "rev_err")) {
      out$cluster_hyp <- hc$v$hypothesis$Est.Error / sqrt(vcl_pop[3])
    } else out$cluster_hyp <- hc$v
    out$cluster_hyp_w <- wtag(hc$w)
  } else out$cluster <- vc$v
  out$cluster_w <- wtag(vc$w)
  # the same call on a fresh fit, on which vcov() has not yet warned
  fit2 <- suppressWarnings(frm(fo, data = d))
  vc2 <- cap(vcov_cluster(fit2, ~ g, full = TRUE))
  out$cluster_fresh_w <- wtag(vc2$w)
  hc2 <- cap(hypothesis(fit2, hs, class = "b", vcov = vc2$v))
  out$cluster_fresh_hyp_w <- wtag(hc2$w)
  # influence: cooks.distance by groups, solve(V) of the whole matrix
  inf <- cap(stats::influence(fit, groups = "g"))
  if (!inherits(inf$v, "rev_err")) {
    cd <- cap(stats::cooks.distance(inf$v))
    Vi <- suppressWarnings(ns$interop_vcov(fit))
    posv <- which(pn %in% c("beta", "betad"))
    stopifnot(ncol(Vi) == length(posv))
    Pe <- mp_inv(Ve[posv, posv])
    dd <- sweep(inf$v$fixed, 2, inf$v$fixed_full)
    cde <- vapply(seq_len(nrow(dd)), function(i) {
      di <- mpfr2array(m(dd[i, ]), dim = c(ncol(dd), 1L))
      num(sum(di * (Pe %*% di))) / ncol(dd)
    }, 0)
    if (!inherits(cd$v, "rev_err")) out$cooks <- cd$v / cde else out$cooks <- cd$v
    out$cooks_w <- wtag(c(inf$w, cd$w))
    dfb <- cap(stats::dfbetas(inf$v))
    if (!inherits(dfb$v, "rev_err")) {
      dfe <- sweep(-dd, 2, sqrt(num(diag(Ve[posv, posv]))), `/`)
      out$dfbetas <- as.vector(dfb$v / dfe)
    }
  } else out$influence <- inf$v
  for (nm in names(out)) {
    v <- out[[nm]]
    cat(sprintf("  %-18s %s\n", nm,
                if (is.character(v)) substr(v[1], 1, 100) else
                  if (length(v) == 1L) sprintf("%.9g", v) else ratio_str(v)))
  }
}
