## Can restoring the call break the case `.frm_` existed for?
## Fixed-effect calls named like covariance structures (diag, ar1, toep,
## exp) beside a real structure exp(0 + t | g), then every downstream
## method. Seed 77.   Rscript dev/review-brmsnames-unsub.R base|lane
arm <- commandArgs(trailingOnly = TRUE)[1L]
base <- c("C:/Users/adf44/source/r/rellib-r12",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (identical(arm, "lane")) {
  c("C:/Users/adf44/source/r/wt-brmsnames-lib", base)
} else base)
q <- function(e) suppressWarnings(suppressMessages(e))
q(library(frmtmb))
cat("arm", arm, dirname(system.file(package = "frmtmb")),
    "unsub:", exists("unsub_specials", asNamespace("frmtmb")), "\n")
# source equality of the installed lane functions
if (arm == "lane") {
  ns <- asNamespace("frmtmb")
  src <- new.env()
  for (f in c("R/parse.R", "R/priors.R", "R/brms-shapes.R")) {
    sys.source(f, envir = src, keep.source = FALSE)
  }
  for (fn in c("unsub_specials", "parse_linpred", "prior_coef_brms",
               "brms_fixef_rows")) {
    cat("installed == source:", fn,
        identical(deparse(get(fn, ns)), deparse(get(fn, src))), "\n")
  }
}
set.seed(77)
n <- 300
d <- data.frame(x = runif(n, 0.2, 2), t = rep(1:5, length.out = n),
                g = factor(rep(1:20, each = 15)),
                g2 = factor(rep(1:30, length.out = 300)))
d$pos <- num_factor(d$t)
d$ft <- factor(d$t)
d$y <- 1 + 0.5 * d$x^2 + 0.3 * sqrt(d$x) + 0.2 * exp(d$x) +
  rnorm(20, 0, 0.4)[d$g] + 0.1 * d$t * rnorm(20)[d$g] + rnorm(n, 0, 0.3)
tc <- function(label, e) {
  r <- tryCatch(q(e), error = function(err) paste("ERROR:",
                                                  conditionMessage(err)))
  cat(sprintf("%-28s", label))
  if (is.character(r)) cat(r, "\n") else print(r)
  invisible(r)
}
# User functions named after structures live in the formula's own
# environment, as a user's helper would
mk <- function() {
  ar1 <- function(v) v^2
  toep <- function(v) sqrt(v)
  diag <- function(v) v^3
  y ~ ar1(x) + toep(x) + diag(x) + exp(x) + exp(pos + 0 | g) +
    ar1(0 + ft | g2)
}
form <- mk()
fit <- tc("fit", frm(form, data = d))
if (inherits(fit, "frmtmb_fit")) {
  tc("fixef", round(fixef(fit)[, 1:2], 4))
  tc("block class", vapply(fit$frame$re %||% list(), function(r)
    r$cs %||% r$class %||% NA_character_, ""))
  tc("VarCorr names", names(VarCorr(fit)))
  tc("variables", grep("^(b|sd|cor)_", variables(fit), value = TRUE))
  tc("logLik", logLik(fit))
  nd <- d[1:4, ]
  nd$x <- c(0.5, 1, 1.5, 2)
  tc("predict newdata", fitted(fit, newdata = nd))
  tc("predict re_formula NA", fitted(fit, newdata = nd, re_formula = NA))
  tc("fitted check vs X%*%b", {
    b <- fixef(fit)[, 1]
    X <- cbind(1, nd$x^2, sqrt(nd$x), nd$x^3, exp(nd$x))
    max(abs(drop(X %*% b) - fitted(fit, newdata = nd, re_formula = NA)[, 1]))
  })
  u <- tc("update()", update(fit, . ~ . - toep(x)))
  if (inherits(u, "frmtmb_fit")) tc("update fixef", rownames(fixef(u)))
  u2 <- tc("update newdata", update(fit, newdata = d[1:200, ]))
  if (inherits(u2, "frmtmb_fit")) tc("update nd fixef", rownames(fixef(u2)))
  tc("drop1", drop1(fit))
  tc("ce names", names(conditional_effects(fit)))
  if (requireNamespace("emmeans", quietly = TRUE)) {
    tc("emmeans x", emmeans::emmeans(fit, ~ x, at = list(x = c(0.5, 1.5))))
  }
  if (requireNamespace("marginaleffects", quietly = TRUE)) {
    tc("avg_slopes x", marginaleffects::avg_slopes(fit, variables = "x"))
  }
  tc("default_prior coef", {
    pr <- as.data.frame(default_prior(form, data = d))
    pr$coef[nzchar(pr$coef)]
  })
  tc("set_prior coef ar1x", {
    v <- as.data.frame(validate_prior(set_prior("normal(0, 1)", coef =
                                                   "ar1x"), form, data = d))
    v[v$source == "user", c("class", "coef", "prior")]
  })
  if (requireNamespace("frmtmb.sample", quietly = TRUE)) {
    ds <- tc("frm_sample", frmtmb.sample::frm_sample(fit, chains = 1,
             iter = 60, warmup = 30, refresh = 0, seed = 1))
    if (inherits(ds, "frmtmb_draws")) {
      tc("draws b_", grep("^b_", variables(ds), value = TRUE))
      tc("draws fixef", rownames(fixef(ds)))
      tc("draws epred nd", dim(posterior_epred(ds, newdata = nd)))
    }
  }
  tc("saveRDS/readRDS predict", {
    f <- tempfile(fileext = ".rds")
    saveRDS(fit, f)
    f2 <- readRDS(f)
    fitted(f2, newdata = nd, re_formula = NA)[, 1]
  })
}
# The structure name used only as a structure, with the base function in
# scope: diag(0 + t | g) next to diag-free fixed terms
fit2 <- tc("diag struct + exp fixed",
           frm(y ~ exp(x) + diag(1 + t | g) + toep(0 + ft | g2), data = d))
if (inherits(fit2, "frmtmb_fit")) {
  tc("fit2 fixef", rownames(fixef(fit2)))
  tc("fit2 VarCorr", rownames(VarCorr(fit2)$g$sd))
  tc("fit2 cor present", !is.null(VarCorr(fit2)$g$cor))
  tc("fit2 predict", fitted(fit2, newdata = d[1:3, ])[, 1])
}
# exp() inside a bar whose wrapper is a structure
fit3 <- tc("ar1 struct around exp",
           frm(y ~ x + ar1(0 + factor(t) | g), data = d))
if (inherits(fit3, "frmtmb_fit")) tc("fit3 VarCorr", names(VarCorr(fit3)))
fit4 <- tc("diag(1 + exp(x) | g)", frm(y ~ exp(x) + diag(1 + exp(x) | g),
                                      data = d))
if (inherits(fit4, "frmtmb_fit")) {
  tc("fit4 VarCorr", rownames(VarCorr(fit4)$g$sd))
  tc("fit4 cor (diag => NULL)", is.null(VarCorr(fit4)$g$cor))
  tc("fit4 re_formula", fitted(fit4, newdata = d[1:3, ],
                                re_formula = ~ (1 + exp(x) | g))[, 1])
  tc("fit4 ranef cols", colnames(ranef(fit4)$g))
}
fit5 <- tc("(1 + exp(x) | g) re_formula",
           frm(y ~ exp(x) + (1 + exp(x) | g), data = d))
if (inherits(fit5, "frmtmb_fit")) {
  full <- fitted(fit5, newdata = d[1:3, ])[, 1]
  tc("re_formula same as full", all.equal(full, predict(fit5,
       newdata = d[1:3, ], re_formula = ~ (1 + exp(x) | g))[, 1]))
}
