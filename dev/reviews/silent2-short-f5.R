# Short check of lane silent2's fix round, F5: aterm_scalar().
#   Rscript dev/reviews/silent2-short-f5.R base|lane
arm <- commandArgs(TRUE)[1]
libs <- c("C:/Users/adf44/source/r/rellib-r21",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (identical(arm, "lane")) libs <- c("C:/Users/adf44/source/r/wt-silent2-lib", libs)
.libPaths(libs)
suppressPackageStartupMessages(library(frmtmb))
cat("arm", arm, "frmtmb", find.package("frmtmb"), "\n")
q <- function(...) suppressWarnings(suppressMessages(frm(...)))
run <- function(e) {
  tryCatch(withCallingHandlers(e, message = function(m) {
    invokeRestart("muffleMessage")
  }), error = function(err) {
    structure(substr(gsub("\n", " ", conditionMessage(err)), 1, 140),
              class = "refused")
  })
}
say <- function(lab, r, ref = NULL) {
  out <- if (inherits(r, "refused")) paste("REFUSED:", r) else if (!is.null(ref)) {
    sprintf("ran, max|diff| vs reference %.3g", max(abs(r - ref)))
  } else "ran"
  cat(sprintf("%-52s %s\n", lab, out))
}
scaled <- frmtmb_family(
  "scaled_normal", dpars = c("mu", "sigma"),
  links = list(mu = "identity", sigma = "log"),
  lpdf = function(y, dpars, aterms) {
    m <- dpars[["mu"]] * aterms[["vreal1"]]
    -0.5 * ((y - m) / dpars[["sigma"]])^2 - log(dpars[["sigma"]]) -
      0.5 * log(2 * pi)
  },
  init_dpars = list(mu = function(y, aterms) mean(y / aterms[["vreal1"]]),
                    sigma = function(y, aterms) stats::sd(y)),
  post = list(mean_fn = function(dpars, aterms) dpars[["mu"]] * aterms[["vreal1"]]),
  sim = function(dpars, aterms, n) {
    stats::rnorm(n, dpars[["mu"]] * aterms[["vreal1"]], dpars[["sigma"]])
  })
set.seed(3)
n <- 120
da <- data.frame(x = rnorm(n))
da$y <- 2 * (1 + 0.5 * da$x) + rnorm(n, 0, 0.3)
est <- function(r) if (inherits(r, "refused")) r else unname(r[, "Estimate"])

cat("## vreal(rep(2, 120))\n")
fk <- q(bf(y | vreal(rep(2, 120)) ~ x), family = scaled, data = da)
ins <- unname(fitted(fk)[, "Estimate"])
say("fitted(newdata = data), 120 rows", est(run(fitted(fk, newdata = da))), ins)
say("fitted(newdata = data[1:10, ]), 10 rows",
    est(run(fitted(fk, newdata = da[1:10, ]))), ins[1:10])
say("predict(newdata = data[1:10, ])",
    run(predict(fk, newdata = da[1:10, ])[, "Estimate"]))
say("posterior_predict(newdata = data[1:10, ])",
    run(posterior_predict(fk, newdata = da[1:10, ], ndraws = 5)))
say("log_lik(newdata = data[1:10, ])",
    run(log_lik(fk, newdata = da[1:10, ], ndraws = 5)))
ce120 <- run(conditional_effects(fk, effects = "x", resolution = 120))
say("CE resolution = 120", if (inherits(ce120, "refused")) ce120 else {
  g <- ce120[[1]]
  g$estimate__ - unname(fitted(fk, newdata = data.frame(x = g$x))[, "Estimate"])
}, 0)
say("CE default", run(conditional_effects(fk, effects = "x")))
# the same model with the constant written as a scalar
fs <- q(bf(y | vreal(2) ~ x), family = scaled, data = da)
say("vreal(2) fit, coef vs vreal(rep(2,120))",
    unname(fixef(fs)[, "Estimate"]), unname(fixef(fk)[, "Estimate"]))

cat("## one distinct element and NA\n")
fna <- run(q(bf(y | vreal(c(rep(2, 110), rep(NA, 10))) ~ x), family = scaled,
             data = da))
say("fit vreal(c(rep(2,110), rep(NA,10)))",
    if (inherits(fna, "refused")) fna else "ok")
if (!inherits(fna, "refused")) {
  cat("  nobs", nobs(fna), "\n")
  say("  fitted(newdata = data)", est(run(fitted(fna, newdata = da))))
  say("  fitted(newdata = data[1:10, ])",
      est(run(fitted(fna, newdata = da[1:10, ]))))
  say("  CE default", run(conditional_effects(fna, effects = "x")))
}

cat("## factor constant\n")
ff <- run(q(bf(y | vreal(factor(rep("a", 120))) ~ x), family = scaled,
            data = da))
say("fit vreal(factor(rep(\"a\", 120)))",
    if (inherits(ff, "refused")) ff else "ok")
if (!inherits(ff, "refused")) {
  insf <- unname(fitted(ff)[, "Estimate"])
  say("  fitted(newdata = data)", est(run(fitted(ff, newdata = da))), insf)
  say("  fitted(newdata = data[1:10, ])",
      est(run(fitted(ff, newdata = da[1:10, ]))), insf[1:10])
  say("  CE default", run(conditional_effects(ff, effects = "x")))
}

cat("## trials(rep(10L, 60)) on 7 rows\n")
set.seed(5)
db <- data.frame(x = rnorm(60))
db$s <- rbinom(60, 10, plogis(db$x))
fb1 <- q(s | trials(rep(10L, 60)) ~ x, family = binomial(), data = db)
fb0 <- q(s | trials(10) ~ x, family = binomial(), data = db)
say("trials(rep(10L,60)) fitted(newdata 7) vs trials(10)",
    est(run(fitted(fb1, newdata = db[1:7, ]))),
    unname(fitted(fb0, newdata = db[1:7, ])[, "Estimate"]))

cat("## a varying variable-free value still refuses\n")
fe <- q(bf(y | vreal(rep_len(c(1, 3), 120)) ~ x), family = scaled, data = da)
say("rep_len(c(1,3),120) fitted(newdata = data)",
    est(run(fitted(fe, newdata = da))))
say("rep_len(c(1,3),120) CE resolution = 120",
    run(conditional_effects(fe, effects = "x", resolution = 120)))
cat("DONE\n")
