# Reviewer: a third-party structure written from ?frmtmb_register_covstruct
# and ?frmtmb-extension-api only. myar(0 + t | g): an AR(1) covariance
# over the d terms of a level, theta = (log sd, atanh rho). It is the
# same model as core's ar1(0 + t | g), so the logLik has a reference.
.libPaths(c("C:/Users/adf44/source/r/wt-spatial-lib",
            "C:/Users/adf44/source/r/rellib-r11",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", as.character(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "\n")
show <- function(label, expr) {
  r <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    cat("  [warning]", label, ":", conditionMessage(w), "\n")
    invokeRestart("muffleWarning")
  }), error = function(e) {
    cat("ERROR", label, ":", conditionMessage(e), "\n")
    NULL
  })
  if (!is.null(r)) cat("OK   ", label, "\n")
  invisible(r)
}
ar_S <- function(theta, d) {
  rho <- tanh(theta[2])
  D <- abs(outer(seq_len(d), seq_len(d), "-"))
  exp(2 * theta[1]) * rho^D
}
frmtmb_register_covstruct(
  "myar",
  lpdf = function(b, theta, blk) {
    d <- blk[["dim"]]
    rho <- tanh(theta[2])
    D <- abs(outer(seq_len(d), seq_len(d), "-"))
    S <- exp(2 * theta[1]) * rho^D
    S <- RTMB::matrix(S, d, d)
    dim(b) <- c(d, length(b) %/% d)
    sum(RTMB::dmvnorm(t(b), 0, S, log = TRUE))
  },
  vcov = function(theta, blk) {
    V <- ar_S(theta, blk[["dim"]])
    dimnames(V) <- list(blk[["cnms"]], blk[["cnms"]])
    V
  },
  npar = 2L, sd_idx = 1L)
cat("registered:", "myar" %in% names(frmtmb:::cs_names()) ||
      "myar" %in% frmtmb:::cs_names(), "\n")

set.seed(42)
G <- 40; d <- 3
dat <- expand.grid(t = factor(seq_len(d)), g = factor(seq_len(G)))
S <- ar_S(c(log(0.8), atanh(0.6)), d)
B <- t(chol(S)) %*% matrix(rnorm(d * G), d)
dat$x <- rnorm(nrow(dat))
dat$y <- 1 + 0.5 * dat$x + B[cbind(as.integer(dat$t), as.integer(dat$g))] +
  rnorm(nrow(dat), sd = 0.5)

f1 <- show("fit myar", frm(y ~ x + myar(0 + t | g), data = dat))
f0 <- show("fit core ar1", frm(y ~ x + ar1(0 + t | g), data = dat))
if (!is.null(f1) && !is.null(f0)) {
  cat("logLik myar", format(as.numeric(logLik(f1)), digits = 12),
      " core ar1", format(as.numeric(logLik(f0)), digits = 12), "\n")
  cat("fixef myar", format(fixef(f1)[, 1], digits = 8), "\n")
  cat("fixef ar1 ", format(fixef(f0)[, 1], digits = 8), "\n")
}
fit <- f1
show("print VarCorr", print(VarCorr(fit)))
show("summary", capture.output(summary(fit)))
show("ranef", ranef(fit))
nd <- data.frame(t = factor(1:3, levels = 1:3), g = factor("new"), x = 0)
p_new <- show("predict new level", predict(fit, newdata = nd,
                                         allow_new_levels = TRUE))
if (!is.null(p_new)) print(p_new)
p_newu <- show("predict new level sample_new_levels=uncertainty",
               predict(fit, newdata = nd, allow_new_levels = TRUE,
                       sample_new_levels = "uncertainty"))
if (!is.null(p_newu)) print(p_newu)
pp <- show("posterior_predict new level", posterior_predict(
  fit, newdata = nd, allow_new_levels = TRUE, ndraws = 2000, seed = 1))
if (!is.null(pp)) {
  cat("pp var by t:", format(apply(pp, 2, var), digits = 4),
      " expected:", format(diag(VarCorr(fit)[[1]][[1]]) + sigma(fit)^2,
                           digits = 4), "\n")
}
sim <- show("simulate", simulate(fit, nsim = 2, seed = 3))
if (!is.null(sim)) print(dim(sim))
sim2 <- show("simulate re_formula = NA", simulate(fit, nsim = 2, seed = 3,
                                                 re_formula = NA))
sim3 <- show("simulate newdata new level", simulate(
  fit, nsim = 2, seed = 3, newdata = nd, allow_new_levels = TRUE))
ci <- show("confint", confint(fit))
if (!is.null(ci)) print(ci)
civ <- show("confint_varcorr", confint_varcorr(fit))
if (!is.null(civ)) print(civ)
cip <- show("confint profile", confint(fit, method = "profile"))
if (!is.null(cip)) print(cip)
ce <- show("conditional_effects", conditional_effects(fit, "x"))
if (!is.null(ce)) print(head(ce[[1]], 3))
ce2 <- show("conditional_effects re_formula NULL",
            conditional_effects(fit, "x", re_formula = NULL))
show("frm_compat", print(frm_compat("myar")))
show("variables", print(variables(fit)))
show("ngrps", print(ngrps(fit)))
show("vcov_cluster", vcov_cluster(fit))
show("loo/brms names: parnames", print(head(as.data.frame(fit), 2)))
# combined with other features
show("myar with gr(g, dist = student)", frm(y ~ x + myar(0 + t | gr(g, dist = "student")),
                                            data = dat))
show("myar with gr(g, by = )", {
  dat$f <- factor(as.integer(dat$g) %% 2)
  frm(y ~ x + myar(0 + t | gr(g, by = f)), data = dat)
})
show("myar |ID| with another term", frm(
  bf(y ~ x + myar(0 + t | q | g)) + bf(x ~ 1 + (1 | q | g)) +
    set_rescor(FALSE), data = dat))
show("myar dim 1 quadrature", frm(y ~ x + myar(1 | g), data = dat,
                                  quadrature = TRUE))
show("myar importance", frm(y ~ x + myar(0 + t | g), data = dat,
                            importance = 50))
show("myar autoscale", frm(y ~ x + myar(0 + t | g), data = dat,
                           control = frmtmb_control(autoscale = TRUE)))
show("myar REML", frm(y ~ x + myar(0 + t | g), data = dat, REML = TRUE))
show("myar set_prior sd", frm(y ~ x + myar(0 + t | g), data = dat,
                              prior = set_prior("normal(0, 1)", class = "sd")))
show("myar set_prior cor", frm(y ~ x + myar(0 + t | g), data = dat,
                               prior = set_prior("lkj(2)", class = "cor")))
show("drop1", drop1(fit))
show("update", update(fit, . ~ . - x))
show("hypothesis", hypothesis(fit, "x > 0"))
saveRDS(fit, "C:/Users/adf44/source/r/frmtmb-wt-spatial/dev/review-spatial-log/myar-fit.rds")
