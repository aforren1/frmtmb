# Release reviewer, 0.73.0: each BREAKING bullet of the core and
# frmtmb.sample NEWS, run on a build (rellib-r10 for 0.72.0, rellib-r11
# for 0.73.0), so that each says what changed.
#   Rscript dev/rel073-rev-breaking.R <lib>   (log breaking-<tag>.txt)
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.sample)})
options(frmtmb.notices = FALSE)
cat("lib: frmtmb", format(packageVersion("frmtmb")), "sample",
    format(packageVersion("frmtmb.sample")), find.package("frmtmb"), "\n")
run <- function(label, expr) {
  w <- character()
  r <- tryCatch(withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
  }), error = function(e) structure(conditionMessage(e), class = "err"))
  if (inherits(r, "err")) {
    cat(label, ": ERROR:", gsub("\\s+", " ", substr(r, 1, 230)), "\n")
  } else {
    ll <- if (inherits(r, "frmtmb_fit")) {
      sprintf("logLik %.4f conv %d", as.numeric(logLik(r)), r$opt$convergence)
    } else ""
    cat(label, ": ok", ll, if (length(w)) paste0("| warnings ", length(w),
      ": ", paste(unique(substr(gsub("\\s+", " ", w), 1, 90)),
                  collapse = " || ")), "\n")
  }
  invisible(r)
}
set.seed(20261021)
n <- 120
d <- data.frame(x = rnorm(n))
d$k <- rpois(n, exp(1 + 0.3 * d$x))
cat("\n== core 1: discrete bounds\n")
run("poisson trunc(lb = 0)", frm(k | trunc(lb = 0) ~ x, family = poisson(),
                                 data = d))
d$k2 <- d$k; d$ci <- "interval"
run("poisson interval y2 == y", frm(k | cens(ci, k2) ~ x, family = poisson(),
                                    data = d))
d$kf <- d$k + 0.5; d$cr <- "right"
run("poisson fractional right bound", frm(kf | cens(cr) ~ x,
                                          family = poisson(), data = d))
cat("\n== core 2: trunc() checks\n")
run("poisson trunc(lb = 2), rows below", frm(k | trunc(lb = 2) ~ x,
                                            family = poisson(), data = d))
run("gaussian trunc(lb = 1, ub = 1)", frm(x | trunc(lb = 1, ub = 1) ~ 1,
                                          data = d[abs(d$x - 1) < 1e-9, ,
                                                   drop = FALSE]))
dh <- data.frame(x = rnorm(200))
dh$y <- rnbinom(200, mu = exp(0.5 + 0.3 * dh$x), size = 2)
dh$y[dh$y > 0 & runif(200) < 0.3] <- 0
run("hurdle_negbinomial zeros below trunc(lb = 1)",
    frm(y | trunc(lb = 1) ~ x, family = hurdle_negbinomial(), data = dh))
cat("\n== core 3: the cens() refusal names the family\n")
run("bernoulli cens", {
  d$b <- rbinom(n, 1, 0.5); d$cn <- "none"
  frm(b | cens(cn) ~ x, family = bernoulli(), data = d)
})
cat("\n== core 4: dirichlet draws at small concentration\n")
set.seed(4)
for (phi in c(0.01, 1e-3)) {
  al <- matrix(phi / 3, 20000, 3)
  S <- frmtmb:::simplex_rdirichlet(al, c("a", "b", "c"))
  cat(sprintf("phi %g: NaN rows %d of 20000 | column means %s\n", phi,
              sum(!is.finite(rowSums(S))),
              paste(format(colMeans(S, na.rm = TRUE), digits = 3),
                    collapse = " ")))
}
cat("\n== core 5: quadrature over more than one block\n")
set.seed(5)
nq <- 240
dq <- data.frame(x = rnorm(nq), ga = factor(rep(1:20, each = 12)),
                 gb = factor(rep(1:2, length.out = nq)),
                 h = factor(rep(1:12, 20)), f = factor(rep(c("u", "v"), nq / 2)))
dq$y <- rbinom(nq, 1, plogis(0.3 * dq$x + rnorm(20, 0, 0.8)[dq$ga]))
dq$yg <- rnorm(nq, 0.3 * dq$x + rnorm(20, 0, 0.8)[dq$ga])
r <- run("bernoulli (1 | ga/gb) quadrature",
         frm(y ~ x + (1 | ga/gb), family = bernoulli(), data = dq,
             quadrature = TRUE))
if (inherits(r, "frmtmb_fit")) {
  fl <- suppressWarnings(frm(y ~ x + (1 | ga/gb), family = bernoulli(),
                             data = dq))
  cat("   quadrature minus Laplace logLik:",
      format(as.numeric(logLik(r)) - as.numeric(logLik(fl)), digits = 3), "\n")
  run("   summary() of it", capture.output(summary(r)))
}
run("bernoulli (1 | ga) + (1 | h) quadrature",
    frm(y ~ x + (1 | ga) + (1 | h), family = bernoulli(), data = dq,
        quadrature = TRUE))
run("gaussian (1 | ga) in mu and sigma, quadrature",
    frm(bf(yg ~ x + (1 | ga), sigma ~ (1 | ga)), data = dq, quadrature = TRUE))
run("bernoulli (1 | gr(ga, by = f)) quadrature",
    frm(y ~ x + (1 | gr(ga, by = f)), family = bernoulli(), data = dq,
        quadrature = TRUE))
run("bernoulli (1 | ga) quadrature", frm(y ~ x + (1 | ga), family = bernoulli(),
                                         data = dq, quadrature = TRUE))
cat("\n== core 6: a weighted likelihood with no maximum\n")
set.seed(6)
nw <- 80
dw <- data.frame(x = rnorm(nw), obs = factor(seq_len(nw)), w = 1.5,
                 f = factor(rep(c("a", "b"), nw / 2)), z = rnorm(nw))
dw$y <- 1 + 0.5 * dw$x + rnorm(nw)
dw$gy <- rgamma(nw, 3, 3 / exp(0.3 + 0.2 * dw$x))
run("gaussian w 1.5 (1 | obs)", frm(y | weights(w) ~ x + (1 | obs), data = dw))
run("gaussian w 1.5 (1 | obs), sigma ~ 0 + x (bounded)",
    frm(bf(y | weights(w) ~ x + (1 | obs), sigma ~ 0 + z), data = dw))
run("Gamma w 1.5 (1 | obs), gamma(2, 0.1) on shape (exempt)",
    frm(gy | weights(w) ~ x + (1 | obs), family = Gamma(link = "log"),
        data = dw, prior = set_prior("gamma(2, 0.1)", class = "shape")))
run("gaussian w 1.5 (1 | obs), student_t on sigma (refused)",
    frm(y | weights(w) ~ x + (1 | obs), data = dw,
        prior = set_prior("student_t(3, 0, 2.5)", class = "sigma")))
run("gaussian w 1.5 (1 | obs), sigma ~ 0 + f, normal on b_sigma (exempt)",
    frm(bf(y | weights(w) ~ x + (1 | obs), sigma ~ 0 + f), data = dw,
        prior = set_prior("normal(0, 1)", class = "b", dpar = "sigma")))
dw$w1 <- 0.9
run("gaussian w 0.9 (1 | obs)", frm(y | weights(w1) ~ x + (1 | obs),
                                    data = dw))
cat("\n== core 7: asym_laplace with an observation-level effect\n")
set.seed(7)
da <- data.frame(x = rnorm(80), obs = factor(1:80))
da$y <- 1 + 0.5 * da$x + rnorm(80)
run("asym_laplace (1 | obs) Laplace", frm(y ~ x + (1 | obs),
                                          family = asym_laplace(), data = da))
run("asym_laplace (1 | obs) quadrature", frm(y ~ x + (1 | obs),
                                             family = asym_laplace(),
                                             data = da, quadrature = TRUE))
cat("\n== sample 3: the moved generics refuse a fit with frmtmb alone\n")
fit <- suppressWarnings(frm(k ~ x, family = poisson(), data = d))
for (g in c("kfold", "reloo", "loo_moment_match", "loo_subsample", "loo_R2",
            "loo_epred")) {
  run(paste0("frmtmb::", g, "(fit)"),
      get(g, envir = asNamespace("frmtmb"))(fit))
}
cat("\n== sample 4: frm_sample() from a formula refuses a weighted spike\n")
run("frm_sample(y | weights(w) ~ x + (1 | obs))",
    frm_sample(bf(y | weights(w) ~ x + (1 | obs)), data = dw, chains = 1,
               iter = 20, refresh = 0))
cat("\n== sample 1 and 2: r_eff on a subset, posterior_average weights\n")
mk <- function(fit, nd = 400L, sd = 0.03, seed = 5L) {
  est <- unlist(fit$estimates[names(fit$frame$par_template)])
  set.seed(seed)
  M <- matrix(rep(est, each = nd) + rnorm(nd * length(est), 0, sd), nd,
              dimnames = list(NULL, frmtmb::brms_par_labels(fit)))
  M <- cbind(frmtmb.sample:::draws_to_natural(M, fit), lp__ = 0)
  structure(list(stanfit = NULL, draws = M, fit = fit),
            class = "frmtmb_draws")
}
d1 <- mk(fit)
d2 <- mk(suppressWarnings(frm(k ~ 1, family = poisson(), data = d)), seed = 6L)
q <- function(e) suppressWarnings(suppressMessages(e))
l_all <- q(loo(d1)); l_sub <- q(loo(d1, ndraws = 100))
cat(sprintf("loo(d1) elpd %.6f | loo(d1, ndraws = 100) elpd %.6f\n",
            l_all$estimates[1, 1], l_sub$estimates[1, 1]))
for (wt in c("stacking", "pseudobma")) {
  pa <- q(posterior_average(d1, d2, weights = wt, ndraws = 400, seed = 1))
  cat(wt, "weights:", format(attr(pa, "weights"), digits = 8), "\n")
}
d1c <- q(add_criterion(d1, "loo", ndraws = 50))
pa <- q(posterior_average(d1c, d2, weights = "stacking", ndraws = 400,
                          seed = 1))
cat("stacking with a stored criterion on 50 draws:",
    format(attr(pa, "weights"), digits = 8), "\n")
cat("END\n")
