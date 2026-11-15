# Lane silent2, item 5207: what separates a flat ordinal-mixture threshold
# from a clean fit, measured on the controls ordstart used: the 40 default
# fits (dev/ordstart-starts.R data, n = 300, cum_cum and cum_sratio seeds
# 1 to 20) and the review's 145-fit sweep (dev/reviews/
# ordstart-review-sweep.R sets sweep, steep and wide). Each fit is the
# default frm() call, as the arm's frmtmb makes it (jittered starts on).
# Per fit and component: the degenerate rule's numbers, the smallest gap
# between ordered thresholds, the largest standard error of a threshold
# and of a coefficient (fixef(), natural units), and per outer parameter
# the log-likelihood change of a shift of +-10 (the smallest over the
# parameters of the larger of the two changes: a parameter that moves 10
# both ways at no cost is on a flat direction).
#   Rscript dev/silent2-ordmeasure.R base|lane set out.tsv
args <- commandArgs(trailingOnly = TRUE)
arm <- args[1]
set <- args[2]
out <- args[3]
libs <- c(if (arm == "lane") "C:/Users/adf44/source/r/wt-silent2-lib",
          "C:/Users/adf44/source/r/rellib-r21",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(libs)
suppressPackageStartupMessages(library(frmtmb))
cat("arm", arm, "set", set, "frmtmb", format(packageVersion("frmtmb")),
    "from", find.package("frmtmb"), "\n")

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
                p = 0.4, noise = "logit", nx = 0, single = FALSE) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n))
  cls <- rbinom(n, 1, p)
  if (single) cls[] <- 1L
  e <- switch(noise, logit = rlogis(n), probit = rnorm(n),
              cloglog = log(rexp(n)))
  lat <- ifelse(cls == 1, b1 * d$x + c1, b2 * d$x + c2) + e
  if (nx > 0) {
    X <- matrix(rnorm(n * nx), n, nx)
    colnames(X) <- paste0("w", seq_len(nx))
    lat <- lat + drop(X %*% rep(0.2, nx))
    d <- cbind(d, X)
  }
  cuts <- stats::quantile(lat, seq_len(K - 1) / K)
  d$y <- 1L + rowSums(outer(lat, cuts, ">"))
  d
}

cases <- list()
add <- function(name, data, formula, family, prior = NULL) {
  cases[[name]] <<- list(data = data, formula = formula, family = family,
                         prior = prior)
}
if (set == "forty") {
  for (fm in c("cum_cum", "cum_sratio")) for (s in 1:20) {
    add(sprintf("%s_seed%d", fm, s),
        substitute(mk(S, mix = TRUE), list(S = s)), bf(y ~ x + z),
        if (fm == "cum_cum") quote(mixture(cumulative(), cumulative())) else
          quote(mixture(cumulative(), sratio())))
  }
} else if (set == "steep") {
  for (n in c(300, 3000, 20000)) for (s in c(3, 5, 8)) for (K in c(3, 6)) {
    for (sd in if (n == 20000) 1 else 1:2) {
      add(sprintf("steep_n%d_s%g_K%d_seed%d", n, s, K, sd),
          substitute(gen(SEED, n = N, K = KK, b1 = S, b2 = -S, c1 = 2,
                         c2 = -2), list(SEED = sd, N = n, KK = K, S = s)),
          bf(y ~ x + z), quote(mixture(cumulative(), cumulative())))
    }
  }
  for (n in c(300, 3000)) for (s in c(3, 8)) for (sd in 1:2) {
    add(sprintf("single_n%d_s%g_seed%d", n, s, sd),
        substitute(gen(SEED, n = N, K = 4, b1 = S, c1 = 0, single = TRUE),
                   list(SEED = sd, N = n, S = s)),
        bf(y ~ x + z), quote(mixture(cumulative(), cumulative())))
  }
} else if (set == "sweep") {
  defs <- list(
    cc_probit = list(K = 5, noise = "probit", f = bf(y ~ x + z),
                     fam = quote(mixture(cumulative("probit"),
                                         cumulative("probit")))),
    sr_cr_cloglog = list(K = 5, noise = "cloglog", f = bf(y ~ x + z),
                         fam = quote(mixture(sratio("cloglog"), cratio()))),
    acat_acat_K3 = list(K = 3, f = bf(y ~ x + z),
                        fam = quote(mixture(acat(), acat()))),
    cum3_K5 = list(K = 5, f = bf(y ~ x + z),
                   fam = quote(mixture(cumulative(), cumulative(),
                                       cumulative()))),
    theta_formula = list(K = 4, f = bf(y ~ x + z, theta2 ~ z),
                         fam = quote(mixture(cumulative(), cumulative()))),
    dpar_formula = list(K = 4, f = bf(y ~ 1, mu1 ~ x, mu2 ~ x + z),
                        fam = quote(mixture(cumulative(), sratio()))),
    cs_sr_acat = list(K = 4, f = bf(y ~ cs(x) + z),
                      fam = quote(mixture(sratio(), acat()))),
    disc_free = list(K = 5, f = bf(y ~ x + z, disc2 ~ 1),
                     fam = quote(mixture(cumulative(), cumulative()))),
    order_mu = list(K = 4, f = bf(y ~ x + z),
                    fam = quote(mixture(cumulative(), sratio(),
                                        order = "mu"))),
    prior_b = list(K = 4, f = bf(y ~ x + z),
                   fam = quote(mixture(cumulative(), cumulative())),
                   prior = quote(c(set_prior("normal(0, 5)", class = "b",
                                             dpar = "mu1"),
                                   set_prior("normal(0, 5)", class = "b",
                                             dpar = "mu2")))),
    cr_cr_probit_K6 = list(K = 6, noise = "probit", f = bf(y ~ x + z),
                           fam = quote(mixture(cratio("probit"),
                                               cratio("probit")))),
    acat_cumcll_K7 = list(K = 7, f = bf(y ~ x + z),
                          fam = quote(mixture(acat(),
                                              cumulative("cloglog")))),
    cc_logit_K7 = list(K = 7, f = bf(y ~ x + z),
                       fam = quote(mixture(cumulative(), cumulative()))))
  for (nm in names(defs)) for (sd in 1:8) {
    df <- defs[[nm]]
    add(sprintf("%s_seed%d", nm, sd),
        substitute(gen(SEED, n = 300, K = KK, noise = NZ),
                   list(SEED = sd, KK = df$K, NZ = df$noise %||% "logit")),
        df$f, df$fam, prior = df$prior)
  }
} else if (set == "wide") {
  for (sd in 1:3) {
    add(sprintf("wide20_K7_n1000_seed%d", sd),
        substitute(gen(SEED, n = 1000, K = 7, nx = 20), list(SEED = sd)),
        bf(as.formula(paste("y ~ x + z +", paste0("w", 1:20,
                                                  collapse = " + ")))),
        quote(mixture(cumulative(), cumulative())))
  }
}

shift_flat <- function(f) {
  p <- f$opt$par
  nll0 <- f$obj$fn(p)
  d <- vapply(seq_along(p), function(i) {
    a <- p; a[i] <- a[i] + 10
    b <- p; b[i] <- b[i] - 10
    v <- abs(c(f$obj$fn(a), f$obj$fn(b)) - nll0)
    if (any(!is.finite(v))) Inf else max(v)
  }, 0)
  f$obj$fn(p)
  i <- which.min(d)
  list(val = d[i], par = paste0(names(p)[i], "[", i, "]"))
}

# per component, the largest change in the mixture's log-likelihood when
# the latent distances to ONE threshold are doubled (through the cs slot,
# as mixture_ord_degeneracy() does for cs() terms): near 0 for a threshold
# the data do not place; NA where every such doubling is not finite
thres_sharp <- function(f, resp) {
  ns <- asNamespace("frmtmb")
  rspec <- f$spec$responses[[resp]]
  fam <- rspec$family
  mx <- fam[["mix"]][["ord"]]
  dp <- ns$with_cs_offsets(f, rspec, ns$eval_dpars(f))[[resp]]
  av <- f$frame[["aterm_values"]][[resp]]
  ex <- ns$fit_extras(f, resp)
  y <- f$frame[["y"]][[resp]]
  n_all <- length(y)
  wts <- rep(as.numeric(av[["weights"]] %||% 1), length.out = n_all)
  wts <- wts / mean(wts)
  ll_at <- function(d) {
    lp_k <- fam[["mix"]][["log_pi"]](d)
    L <- vapply(seq_len(mx$K), function(k) {
      rep(as.numeric(lp_k[[k]]), length.out = n_all) +
        as.numeric(fam[["mix"]][["comp_lpdf"]](y, d, av, k, ex))
    }, numeric(n_all))
    m <- apply(L, 1L, max)
    sum(wts * (m + log(rowSums(exp(L - m)))))
  }
  ll0 <- ll_at(dp)
  vapply(seq_len(mx$K), function(k) {
    lp <- f$frame[["linpreds"]][[ns$linpred_key(resp, paste0("mu", k))]]
    ob <- ns$ord_lp_block(f$frame, f$spec, lp)
    tau <- as.numeric(ns$ord_threshold_values(ob$fam, f$estimates[[ob$comp]]))
    eta <- as.numeric(dp[[paste0("mu", k)]])
    eta <- rep_len(eta, n_all)
    slot <- ns$cs_slot(paste0("mu", k))
    cs <- dp[[slot]] %||% matrix(0, n_all, length(tau))
    if (ncol(cs) != length(tau)) return(NA_real_)
    sh <- vapply(seq_along(tau), function(j) {
      d3 <- dp
      cs3 <- cs
      cs3[, j] <- 2 * cs[, j] + eta - tau[j]
      d3[[slot]] <- cs3
      suppressWarnings(ll_at(d3)) - ll0
    }, 0)
    sh <- sh[is.finite(sh)]
    if (length(sh)) max(sh) else NA_real_
  }, 0)
}

rows <- list()
for (nm in names(cases)) {
  cs <- cases[[nm]]
  d <- eval(cs$data)
  w <- character()
  a <- list(cs$formula, family = eval(cs$family), data = d)
  if (!is.null(cs$prior)) a$prior <- eval(cs$prior)
  t0 <- proc.time()[["elapsed"]]
  f <- tryCatch(withCallingHandlers(do.call(frm, a), warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  }), error = function(e) e)
  el <- proc.time()[["elapsed"]] - t0
  if (inherits(f, "error")) {
    cat(sprintf("CASE %s ERROR %s\n", nm, conditionMessage(f)))
    next
  }
  resp <- names(f$spec$responses)[1]
  dg <- frmtmb:::mixture_ord_degeneracy(f, resp)
  fx <- tryCatch(suppressWarnings(fixef(f)), error = function(e) NULL)
  K <- nrow(dg)
  gap <- se_t <- rep(NA_real_, K)
  if (!is.null(fx)) {
    for (k in seq_len(K)) {
      r <- grep(paste0("^mu", k, "_Intercept\\["), rownames(fx))
      if (length(r)) {
        se_t[k] <- max(fx[r, "Est.Error"])
        if (length(r) > 1L) gap[k] <- min(diff(fx[r, "Estimate"]))
      }
    }
  }
  se_b <- if (is.null(fx)) NA else {
    r <- grep("_Intercept\\[", rownames(fx), invert = TRUE)
    if (length(r)) max(fx[r, "Est.Error"]) else NA
  }
  sf <- shift_flat(f)
  ts <- tryCatch(thres_sharp(f, resp), error = function(e) {
    cat("thres_sharp error:", conditionMessage(e), "\n")
    rep(NA_real_, K)
  })
  row <- data.frame(case = nm, ll = as.numeric(logLik(f)),
                    code = f$opt$convergence, elapsed = el,
                    degen = any(grepl("ended at a degenerate boundary", w)),
                    stepw = any(grepl("jittered starts", w)),
                    nwarn = length(w), se_b = se_b,
                    flat10 = sf$val, flat10_par = sf$par,
                    maxpar = max(abs(f$opt$par)))
  for (k in 1:3) {
    for (v in c("sharpen", "near", "cover", "collapsed")) {
      row[[paste0(v, k)]] <- if (k <= K) dg[[v]][k] else NA
    }
    row[[paste0("gap", k)]] <- if (k <= K) gap[k] else NA
    row[[paste0("se_t", k)]] <- if (k <= K) se_t[k] else NA
    row[[paste0("tsharp", k)]] <- if (k <= K) ts[k] else NA
  }
  row$warn <- paste(substr(gsub("[\t\n]", " ", w), 1, 80), collapse = " | ")
  rows[[nm]] <- row
  cat(sprintf(paste("CASE %s ll %.4f code %d degen %s se_t %s se_b %.3g",
                    "flat10 %.3g (%s) gap %s cover %s near %s tsharp %s\n"),
              nm, row$ll, row$code, row$degen,
              paste(signif(se_t, 3), collapse = "/"), se_b, sf$val, sf$par,
              paste(signif(gap, 3), collapse = "/"),
              paste(signif(dg$cover, 3), collapse = "/"),
              paste(signif(dg$near, 3), collapse = "/"),
              paste(signif(ts, 3), collapse = "/")))
}
R <- do.call(rbind, rows)
utils::write.table(R, out, sep = "\t", row.names = FALSE, quote = TRUE)
cat("DONE", nrow(R), "\n")
