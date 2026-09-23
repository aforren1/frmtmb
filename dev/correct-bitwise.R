# Lane wt-correct: nothing the items do not target may move. Fits every
# untargeted design under ML, REML = TRUE and control(profile = TRUE),
# and saves vcov(), fixef(), logLik(), summary(), predict() at a fixed
# seed, and the frm_sample() default prior's value at the estimate.
#   Rscript dev/correct-bitwise.R base   -> dev/correct-log/bitwise-base.rds
#   Rscript dev/correct-bitwise.R lane   -> dev/correct-log/bitwise-lane.rds
#   Rscript dev/correct-bitwise.R compare
arm <- commandArgs(trailingOnly = TRUE)[1]
OUT <- "C:/Users/adf44/source/r/frmtmb-wt-correct/dev/correct-log"
if (identical(arm, "compare")) {
  a <- readRDS(file.path(OUT, "bitwise-base.rds"))
  b <- readRDS(file.path(OUT, "bitwise-lane.rds"))
  stopifnot(identical(names(a), names(b)))
  n_ok <- 0L
  for (k in names(a)) {
    for (q in names(a[[k]])) {
      # environments are pointers, and two processes never share one:
      # a closure (a family's lpdf, the call's formula) is compared by
      # its code and its values, not by where it lives
      same <- identical(a[[k]][[q]], b[[k]][[q]], ignore.environment = TRUE,
                        ignore.bytecode = TRUE, ignore.srcref = TRUE)
      cat(sprintf("%-34s %-10s %s\n", k, q,
                  if (same) "identical" else "DIFFERS"))
      if (!same && is.list(a[[k]][[q]])) {
        # which components, and whether they differ in value at all
        x <- a[[k]][[q]]
        y <- b[[k]][[q]]
        for (f in union(names(x), names(y))) {
          if (identical(x[[f]], y[[f]], ignore.environment = TRUE,
                        ignore.bytecode = TRUE, ignore.srcref = TRUE)) next
          ae <- all.equal(x[[f]], y[[f]])
          cat(sprintf("    component %-10s all.equal: %s\n", f,
                      if (isTRUE(ae)) "TRUE (an environment pointer)" else
                        paste(length(ae), "differences")))
        }
      }
      n_ok <- n_ok + same
    }
  }
  cat("identical:", n_ok, "of", sum(lengths(a)), "\n")
  quit(save = "no")
}
source("C:/Users/adf44/source/r/frmtmb-wt-correct/dev/correct-prelude.R")
if (identical(arm, "base")) .libPaths(.libPaths()[-1L])
suppressMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("arm", arm, "frmtmb", format(packageVersion("frmtmb")), "from",
    dirname(find.package("frmtmb")), "\n")

set.seed(20260922)
n <- 240
d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:12, 20)))
u <- rnorm(12, 0, 0.6)
v <- rnorm(12, 0, 0.3)
d$y <- 1 + 0.5 * d$x + u[d$g] + v[d$g] * d$x +
  rnorm(n, 0, exp(-0.2 + 0.3 * d$z))
d$cnt <- rpois(n, exp(0.2 + 0.4 * d$x + u[d$g]))
d$o <- factor(cut(0.8 * d$x + rlogis(n), c(-Inf, -0.5, 0.7, Inf),
                  labels = FALSE), ordered = TRUE)
d$ym <- ifelse(rbinom(n, 1, 0.4) == 1, rnorm(n, 3, 0.8), rnorm(n, -1, 0.8))

designs <- list(
  gaussian = bf(y ~ x) + gaussian(),
  poisson = bf(cnt ~ x) + poisson(),
  distributional = bf(y ~ x, sigma ~ z) + gaussian(),
  mixed = bf(y ~ x + (1 + x | g)) + gaussian(),
  mixed_poisson = bf(cnt ~ x + (1 | g)) + poisson(),
  dist_mixed = bf(y ~ x + (1 | g), sigma ~ z) + gaussian()
)
settings <- list(ML = list(), REML = list(REML = TRUE),
                 profile = list(control = frmtmb_control(profile = TRUE)))
res <- list()
# an accessor that errors is recorded by its message, so a refusal that
# moved is a difference too
safe <- function(expr) {
  tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e)))
}
grab <- function(fit) {
  out <- list(vcov = safe(vcov(fit)), fixef = safe(fixef(fit)),
              logLik = safe(logLik(fit)), summary = safe(summary(fit)))
  set.seed(7)
  out$predict <- safe(predict(fit))
  out
}
for (nm in names(designs)) {
  for (st in names(settings)) {
    fit <- tryCatch(suppressWarnings(do.call(frm, c(list(designs[[nm]],
                                                          data = d),
                                                     settings[[st]]))),
                    error = function(e) e)
    key <- paste(nm, st)
    res[[key]] <- if (inherits(fit, "error")) {
      list(error = conditionMessage(fit))
    } else {
      grab(fit)
    }
    cat(key, if (inherits(fit, "error")) conditionMessage(fit) else "ok",
        "\n")
  }
}
# untargeted designs whose code paths the lane touched: a mixture with no
# theta formula (mixture() was refactored), theta1 ~ x on two components
# (the reference stays the last one), mo() (brms_coef_table), an ordinal
# fit (residuals(type = "osa") and fitted() stay), a MAP fit whose sd
# prior has only location blocks to reach
extra <- list(
  mixture = function() frm(bf(ym ~ 1) + mixture(gaussian(), gaussian()),
                           data = d),
  mixture_theta1 = function() frm(bf(ym ~ 1, theta1 ~ x) +
                                    mixture(gaussian(), gaussian()),
                                  data = d),
  mo = function() {
    d$m <- factor(sample(1:4, n, TRUE), ordered = TRUE)
    frm(bf(y ~ mo(m) + x) + gaussian(), data = d)
  },
  ordinal = function() frm(bf(o ~ x) + cumulative(), data = d),
  map_sd = function() frm(bf(y ~ x + (1 | g)) + gaussian(), data = d,
                          prior = set_prior("normal(0, 1)", class = "sd"))
)
for (nm in names(extra)) {
  fit <- suppressWarnings(extra[[nm]]())
  r <- list(vcov = vcov(fit), fixef = fixef(fit), logLik = logLik(fit),
            summary = summary(fit), variables = variables(fit))
  if (nm == "ordinal") {
    r$osa <- residuals(fit, type = "osa")
    r$fitted <- fitted(fit)
  } else {
    r$residuals <- residuals(fit)
  }
  res[[nm]] <- r
  cat(nm, "ok\n")
}
# frm_sample()'s default priors on untargeted univariate designs: the
# negative log prior at the estimate, and the specs themselves
# a sigma block, where the per-prefix sd default must give the density
# the one class-wide default gave
designs$sigma_re <- bf(y ~ x + (1 | g), sigma ~ z + (1 | g)) + gaussian()
for (nm in c("gaussian", "mixed", "mixed_poisson", "dist_mixed",
             "sigma_re")) {
  uf <- frm(designs[[nm]], data = d, dry_run = "objective")
  fit <- suppressWarnings(frm(designs[[nm]], data = d))
  ri <- suppressMessages(frmtmb.sample:::sample_resolve_priors(uf, NULL))$ri
  res[[paste("sample defaults", nm)]] <- list(
    nlp = neg_log_prior_fn(ri$entries)(fit$estimates))
}
saveRDS(res, file.path(OUT, paste0("bitwise-", arm, ".rds")))
cat("saved", length(res), "designs\n")
