# Review of lane meme: the L1 identity on cases the lane did not run.
#   Rscript dev/review-meme-l1.R  -> dev/review-meme-out/l1.txt
# A: three ungrouped me() terms, one 3 x 3 correlation (the lane's L1
#    has 2 x 2 correlations only), on
#    A1 frm()'s MAP objective with lkj(2) on corme,
#    A2 the frm_sample() objective with the same priors,
#    A3 the frm_sample() objective under its defaults, against brms's
#       defaults (lkj(1) on corme, flat meanme and sdme),
#    A4 the frm_sample() objective under prior = "flat", against brms
#       with "" on every class (corme included).
# B: me(xg, gr = g) + me(z) + me(wg, gr = g) + me(x): two correlated
#    groups, interleaved, under frm_sample() defaults and under lkj(3)
#    on corme with lkj(2) on group g.
# -f(theta) = log_prob(U(theta)) + log|det dU/dtheta| + C, as in
# tests/testthat/helper-me-lp.R; the closed form of log|det J| here is
# general in the size of each correlated group.
.libPaths(c("C:/Users/adf44/source/r/wt-meme-lib",
            "C:/Users/adf44/source/r/rellib-r15",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-meme/dev/stan-cache",
           NOT_CRAN = "true")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
  library(testthat)
})
cat("lib:", find.package("frmtmb"), find.package("frmtmb.sample"), "\n")
tt <- "C:/Users/adf44/source/r/frmtmb-wt-meme/tests/testthat"
te <- new.env(parent = asNamespace("frmtmb"))
for (h in c("helper-brms.R", "helper-me-lp.R")) {
  sys.source(file.path(tt, h), envir = te)
}
sns <- asNamespace("frmtmb.sample")

rv_d <- function(seed = 77, n = 40) {
  set.seed(seed)
  tx <- rnorm(n, 1, 0.8)
  tz <- 0.5 * tx + rnorm(n, 0, 0.7)
  tw <- -0.4 * tx + rnorm(n, 0, 0.7)
  d <- data.frame(x = tx + rnorm(n, 0, 0.3), sx = 0.3,
                  z = tz + rnorm(n, 0, 0.3), sz = runif(n, 0.2, 0.4),
                  w = tw + rnorm(n, 0, 0.3), sw = 0.35,
                  g = factor(rep(1:8, each = 5)))
  d$xg <- rnorm(8)[d$g]
  d$wg <- 0.5 * d$xg + rnorm(8, 0, 0.6)[d$g]
  d$sg <- 0.25
  d$y <- 2 + 0.7 * tx - 0.4 * tz + 0.3 * tw + 0.4 * d$xg +
    rnorm(n, 0, 0.5)
  d
}

# log |det dU/dtheta| in closed form: per correlated group of size d
# over N latent rows, -N sum(log sdme) - N sum(log diag L) plus the
# chart-to-CPC part, -1/2 sum_{i>j} log(sum_{k>=j} a_ik^2), where a is
# the unit-diagonal lower triangle us_chol_L() normalizes by rows
lfun_for <- function(me) {
  function(p, cmp) {
    lsd <- p[cmp == "logsdme"]
    th <- p[cmp == "thetame"]
    out <- 0
    for (g in me$groups) {
      K <- g$K
      if (isTRUE(g$cor)) {
        dd <- length(K)
        A <- diag(dd)
        A[lower.tri(A)] <- th[g$th_idx]
        N <- length(me$terms[[K[1]]]$idx)
        out <- out - N * sum(lsd[K]) + N * 0.5 * sum(log(rowSums(A^2)))
        for (i in 2:dd) for (j in 1:(i - 1)) {
          out <- out - 0.5 * log(sum(A[i, j:i]^2))
        }
      } else {
        for (k in K) out <- out - length(me$terms[[k]]$idx) * lsd[k]
      }
    }
    out
  }
}

points_check <- function(label, fit, ours_f, p0, bform, d, bp, npts = 4,
                         seed = 23, h = 1e-6) {
  cat("====", label, "\n")
  bfam <- stats::gaussian()
  code <- te$with_brms_me(function() {
    brms::make_stancode(bform, data = d, family = bfam, prior = bp)
  })
  pl <- grep("lprior \\+=|target \\+=", strsplit(code, "\n")[[1]],
             value = TRUE)
  cat("brms statements:\n", paste0("  ", trimws(pl), "\n"), sep = "")
  sdat <- te$with_brms_me(function() {
    te$brms_standata(bform, data = d, family = bfam, prior = bp)
  })
  rtab <- te$with_brms_me(function() te$brms_ranef_table(bform, d, bfam, bp))
  mod <- te$brms_stan_model(code)
  sf <- suppressMessages(rstan::sampling(mod, data = sdat, chains = 0))
  cmp <- names(p0)
  at <- function(p) {
    f2 <- fit
    for (cp in unique(cmp)) f2$estimates[[cp]][] <- p[cmp == cp]
    f2
  }
  U <- function(p) {
    rstan::unconstrain_pars(sf, te$stan_pars_from_fit(at(p), sdat, code,
                                                      rtab))
  }
  if (length(U(p0)) != length(p0)) stop("parameter counts differ")
  Lfun <- lfun_for(fit$frame$me)
  step <- function(p, j) replace(numeric(length(p)), j, h * max(1, abs(p[j])))
  set.seed(seed)
  pts <- c(list(p0), lapply(seq_len(npts - 1L), function(k) {
    p0 + stats::rnorm(length(p0), 0, 0.1)
  }))
  rows <- lapply(pts, function(p) {
    u <- U(p)
    J <- vapply(seq_along(p), function(j) {
      e <- step(p, j)
      (U(p + e) - U(p - e)) / (2 * e[j])
    }, numeric(length(p)))
    l_num <- as.numeric(determinant(J, logarithm = TRUE)$modulus)
    l_ana <- Lfun(p, cmp)
    ours <- -ours_f(p)
    lpb <- rstan::log_prob(sf, u, adjust_transform = TRUE)
    g_ours <- -as.numeric(ours_f(p, order = 1))
    g_b <- as.numeric(rstan::grad_log_prob(sf, u, adjust_transform = TRUE))
    g_l <- vapply(seq_along(p), function(j) {
      e <- step(p, j)
      (Lfun(p + e, cmp) - Lfun(p - e, cmp)) / (2 * e[j])
    }, 0)
    g_pred <- as.numeric(crossprod(J, g_b)) + g_l
    c(ours = ours, brms = lpb, Lnum = l_num, Lana = l_ana,
      const = ours - lpb - l_num, grad_err = max(abs(g_ours - g_pred)),
      grad_scale = max(abs(g_ours)))
  })
  tab <- do.call(rbind, rows)
  print(tab, digits = 12)
  cat(sprintf(paste0("RESULT %s: spread(const)/scale %.3g  max grad rel ",
                     "%.3g  max |Lnum - Lana| %.3g  const/log(2) %.10f\n"),
              label,
              diff(range(tab[, "const"])) /
                max(abs(tab[, c("ours", "brms", "Lnum")])),
              max(tab[, "grad_err"] / pmax(1, tab[, "grad_scale"])),
              max(abs(tab[, "Lnum"] - tab[, "Lana"]) /
                    pmax(1, abs(tab[, "Lnum"]))),
              mean(tab[, "const"]) / log(2)))
  invisible(tab)
}

sample_f <- function(fit0, prior) {
  rs <- suppressMessages(sns$sample_resolve_priors(fit0, prior))
  obj <- sns$prior_augmented_obj(fit0, rs$ri$entries)
  stopifnot(identical(names(obj$env$last.par), names(fit0$obj$env$last.par)))
  cat("sampled entries:", vapply(rs$ri$entries, function(e) {
    paste0(e$comp, "[", paste(e$idx, collapse = ","), "]:", e$dist$kind,
           if (!is.null(e$dist$eta)) e$dist$eta,
           if (isTRUE(e$dist$flat)) "(flat)")
  }, ""), "\n", fill = 78)
  obj$env$f
}

d <- rv_d()
base_rows <- list(list(prior = "normal(0, 1)", class = "b"),
                  list(prior = "student_t(3, 0, 2.5)", class = "Intercept"),
                  list(prior = "student_t(3, 0, 2.5)", class = "sigma"),
                  list(prior = "normal(0, 2)", class = "meanme"),
                  list(prior = "exponential(1)", class = "sdme"))

## A
fA <- bf(y ~ me(x, sx) + me(z, sz) + me(w, sw))
bA <- brms::bf(y ~ me(x, sx) + me(z, sz) + me(w, sw))
prA <- te$me_lp_priors(c(base_rows, list(list(prior = "lkj(2)",
                                              class = "corme"))))
fitA <- frm(fA, family = gaussian(), data = d, prior = prA$frm,
            control = frmtmb_control(grad_tol = 1e-6, restarts = 3))
p0A <- fitA$obj$env$last.par.best
points_check("A1 frm() MAP, lkj(2) on a 3 x 3 corme", fitA, fitA$obj$env$f,
             p0A, bA, d, prA$brm)
fitA0 <- frm(fA, family = gaussian(), data = d, dry_run = "objective")
points_check("A2 frm_sample(), the same priors", fitA,
             sample_f(fitA, prA$frm), p0A, bA, d, prA$brm)
points_check("A3 frm_sample() defaults vs brms defaults", fitA,
             sample_f(fitA, NULL), p0A, bA, d, NULL)
bflat <- c(brms::set_prior("", class = "Intercept"),
           brms::set_prior("", class = "sigma"),
           brms::set_prior("", class = "corme"))
points_check("A4 frm_sample(prior = \"flat\") vs brms all \"\"", fitA,
             sample_f(fitA, "flat"), p0A, bA, d, bflat)

## B
fB <- bf(y ~ me(xg, sg, gr = g) + me(z, sz) + me(wg, sg, gr = g) + me(x, sx))
bB <- brms::bf(y ~ me(xg, sg, gr = g) + me(z, sz) + me(wg, sg, gr = g) +
                 me(x, sx))
prB <- te$me_lp_priors(c(base_rows, list(
  list(prior = "lkj(3)", class = "corme"),
  list(prior = "lkj(2)", class = "corme", group = "g"),
  list(prior = "normal(1, 0.5)", class = "meanme", coef = "mewg"),
  list(prior = "student_t(3, 0, 1)", class = "sdme", coef = "mez"))))
fitB <- frm(fB, family = gaussian(), data = d, prior = prB$frm,
            control = frmtmb_control(grad_tol = 1e-6, restarts = 3))
cat("B groups:", vapply(fitB$frame$me$groups, function(g) {
  paste0("[gr=", g$gr, " K=", paste(g$K, collapse = ","), " cor=", g$cor,
         "]")
}, ""), "\n")
p0B <- fitB$obj$env$last.par.best
points_check("B1 frm() MAP, two correlated groups interleaved", fitB,
             fitB$obj$env$f, p0B, bB, d, prB$brm)
fitB0 <- frm(fB, family = gaussian(), data = d, dry_run = "objective")
points_check("B2 frm_sample() defaults vs brms defaults", fitB,
             sample_f(fitB, NULL), p0B, bB, d, NULL)
cat("DONE\n")
