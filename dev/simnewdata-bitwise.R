# Lane wt-simnewdata: nothing the items do not target may move. Fits
# designs under ML, REML = TRUE and control(profile = TRUE), and saves
# vcov(), fixef(), logLik(), summary(), predict() and fitted() at a
# fixed seed, the caller's .Random.seed after each seeded call, and
# simulate() under the re_formula values whose meaning did NOT change
# (NULL, and NA on a fit with no population smooth), plus the partial
# re_formula paths of predict()/fitted()/frm_linpred(), which now run
# through the split re_keep_plan()/re_resolve().
#   Rscript dev/simnewdata-bitwise.R base   -> dev/simnewdata-log/bitwise-base.rds
#   Rscript dev/simnewdata-bitwise.R lane   -> dev/simnewdata-log/bitwise-lane.rds
#   Rscript dev/simnewdata-bitwise.R compare
arm <- commandArgs(trailingOnly = TRUE)[1]
OUT <- "C:/Users/adf44/source/r/frmtmb-wt-simnewdata/dev/simnewdata-log"
# what each item may move, by the name of the quantity; everything else
# must be identical
intended <- c("simulate NA smooth", "simulate ~1", "simulate partial",
              "pp_check ~1")
# a draw that now redraws group effects takes more of the caller's
# stream, so the stream after it moves with it
intended <- c(intended, paste(intended, "seed"))
if (identical(arm, "compare")) {
  a <- readRDS(file.path(OUT, "bitwise-base.rds"))
  b <- readRDS(file.path(OUT, "bitwise-lane.rds"))
  stopifnot(identical(names(a), names(b)))
  tally <- c(identical = 0L, moved_as_intended = 0L, DIFFERS = 0L)
  for (k in names(a)) {
    stopifnot(identical(names(a[[k]]), names(b[[k]])))
    for (q in names(a[[k]])) {
      same <- identical(a[[k]][[q]], b[[k]][[q]], ignore.environment = TRUE,
                        ignore.bytecode = TRUE, ignore.srcref = TRUE)
      v <- if (same) "identical" else if (q %in% intended) {
        "moved_as_intended"
      } else "DIFFERS"
      tally[[v]] <- tally[[v]] + 1L
      cat(sprintf("%-26s %-22s %s\n", k, q, v))
    }
  }
  print(tally)
  quit(save = "no")
}
source("dev/simnewdata-prelude.R")
if (identical(arm, "base")) .libPaths(.libPaths()[-1L])
suppressMessages(library(frmtmb))
cat("arm", arm, "frmtmb", format(packageVersion("frmtmb")), "from",
    dirname(find.package("frmtmb")), "\n")

set.seed(20260923)
n <- 240
d <- data.frame(x = rnorm(n), z = rnorm(n), w = runif(n),
                g = factor(rep(1:12, 20)), h = factor(rep(1:8, 30)))
u <- rnorm(12, 0, 0.6)
v <- rnorm(12, 0, 0.3)
d$y <- 1 + 0.5 * d$x + u[d$g] + v[d$g] * d$x + rnorm(8, 0, 0.4)[d$h] +
  sin(4 * d$w) + rnorm(n, 0, exp(-0.2 + 0.3 * d$z))
d$cnt <- rpois(n, exp(0.2 + 0.4 * d$x + u[d$g]))

designs <- list(
  gaussian = bf(y ~ x) + gaussian(),
  poisson = bf(cnt ~ x) + poisson(),
  distributional = bf(y ~ x, sigma ~ z) + gaussian(),
  mixed = bf(y ~ x + (1 + x | g)) + gaussian(),
  two_terms = bf(y ~ x + (1 | g) + (1 | h)) + gaussian(),
  mixed_poisson = bf(cnt ~ x + (1 | g)) + poisson(),
  dist_mixed = bf(y ~ x + (1 | g), sigma ~ z) + gaussian(),
  smooth_mixed = bf(y ~ x + s(w) + (1 | g)) + gaussian()
)
settings <- list(ML = list(), REML = list(REML = TRUE),
                 profile = list(control = frmtmb_control(profile = TRUE)))
safe <- function(expr) {
  tryCatch(expr, error = function(e) paste("ERROR:", conditionMessage(e)))
}
seeded <- function(out, nm, expr) {
  set.seed(7)
  out[[nm]] <- safe(eval.parent(substitute(expr)))
  out[[paste(nm, "seed")]] <- get(".Random.seed", envir = globalenv())
  out
}
grab <- function(fit, nm) {
  has_re <- length(fit$frame$re_blocks) > 0L
  has_sm <- grepl("smooth", nm)
  two <- identical(nm, "two_terms")
  out <- list(vcov = safe(vcov(fit)), fixef = safe(fixef(fit)),
              logLik = safe(logLik(fit)), summary = safe(summary(fit)))
  out <- seeded(out, "predict", predict(fit, ndraws = 200))
  out <- seeded(out, "fitted", fitted(fit))
  out <- seeded(out, "predict NA", predict(fit, ndraws = 200,
                                           re_formula = NA))
  out <- seeded(out, "fitted ~1", fitted(fit, re_formula = ~1))
  out <- seeded(out, "simulate NULL", simulate(fit, nsim = 3))
  out <- seeded(out, "simulate NULL seed", simulate(fit, nsim = 3,
                                                    seed = 11))
  out <- seeded(out, if (has_sm) "simulate NA smooth" else "simulate NA",
                simulate(fit, nsim = 3, re_formula = NA))
  if (has_re) {
    out <- seeded(out, "simulate ~1", simulate(fit, nsim = 3,
                                               re_formula = ~1))
  }
  if (two) {
    out <- seeded(out, "predict partial",
                  predict(fit, ndraws = 200, re_formula = ~ (1 | g)))
    out <- seeded(out, "fitted partial",
                  fitted(fit, re_formula = ~ (1 | g)))
    out <- seeded(out, "linpred partial nd",
                  frm_linpred(fit, newdata = d[1:20, ],
                              re_formula = ~ (1 | h), se.fit = TRUE))
    out <- seeded(out, "simulate partial",
                  simulate(fit, nsim = 3, re_formula = ~ (1 | g)))
    out <- seeded(out, "nosuch",
                  predict(fit, re_formula = ~ (1 | nosuch)))
  }
  if (has_re && !has_sm) {
    out <- seeded(out, "pp_check NA",
                  pp_check(fit, ndraws = 3)$data)
    out <- seeded(out, "pp_check ~1",
                  pp_check(fit, ndraws = 3, re_formula = ~1)$data)
  }
  out
}
res <- list()
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
      suppressWarnings(grab(fit, nm))
    }
    cat(key, if (inherits(fit, "error")) conditionMessage(fit) else "ok",
        "\n")
  }
}
saveRDS(res, file.path(OUT, paste0("bitwise-", arm, ".rds")))
cat("saved", length(res), "designs\n")
