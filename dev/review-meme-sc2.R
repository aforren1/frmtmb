# Review of lane meme, short check of fix round 2: the nonlinear
# intercept spelling on frm_sample().
#   Rscript dev/review-meme-sc2.R <lane|base> > dev/review-meme-out/sc2/<arm>.txt
arm <- commandArgs(TRUE)[1]
libs <- c("C:/Users/adf44/source/r/rellib-r15",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (identical(arm, "lane")) libs <- c("C:/Users/adf44/source/r/wt-meme-lib", libs)
.libPaths(libs)
Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("lib:", find.package("frmtmb"), find.package("frmtmb.sample"), "\n")
q <- function(expr) suppressMessages(suppressWarnings(expr))
try_ <- function(lab, expr) {
  r <- tryCatch(q(expr), error = function(e) paste("ERROR:", conditionMessage(e)))
  cat("--", lab, "\n")
  if (is.character(r) && length(r) == 1) cat("  ", r, "\n") else print(r)
  invisible(r)
}
set.seed(2)
d <- data.frame(x = runif(60, 0, 3), w = rnorm(60),
                g = factor(rep(1:6, 10)))
d$y <- 2 * exp(-0.7 * d$x) + rnorm(60, 0, 0.1)
f1 <- bf(y ~ a * exp(-b * x), a ~ 1, b ~ 1, nl = TRUE)
fw <- bf(y ~ a * exp(-b * x), a ~ 1 + w, b ~ 1, nl = TRUE)
fg <- bf(y ~ a * exp(-b * x), a ~ 1 + (1 | g), b ~ 1, nl = TRUE)
ents <- function(f, pr) {
  fit <- q(frm(f, data = d, dry_run = "objective"))
  rs <- q(frmtmb.sample:::sample_resolve_priors(fit, pr))
  pt <- fit$frame$par_template
  list(entries = vapply(rs$ri$entries, function(e) {
    nm <- frmtmb:::par_template_names(pt[[e$comp]], e$comp)[e$idx]
    paste0(paste(nm, collapse = ","), ":", e$dist$kind, "(",
           paste(unlist(e$dist[-1]), collapse = ","), ")")
  }, ""), lower = rs$ri$lower, upper = rs$ri$upper)
}
draws <- function(f, pr) {
  ds <- q(frm_sample(f, data = d, family = gaussian(),
                     prior = pr + set_prior("exponential(1)", class = "sigma"),
                     sample_prior = "yes", chains = 2, iter = 600, seed = 3,
                     refresh = 0))
  pd <- prior_draws(ds)
  m <- as.matrix(ds)
  list(prior_names = names(pd), prior_min = round(vapply(pd, min, 0), 4),
       post_min = round(apply(m[, grep("^b_[ab]_", colnames(m)), drop = FALSE],
                              2, min), 4),
       post_mean = round(colMeans(m[, grep("^b_[ab]_", colnames(m)),
                                    drop = FALSE]), 4))
}
pI <- set_prior("normal(1, 1)", class = "Intercept", nlpar = "a", lb = 0) +
  set_prior("normal(1, 1)", class = "Intercept", nlpar = "b", lb = 0.4)
cat("#### 1. S1 rerun, a ~ 1\n")
try_("draws, Intercept nlpar spelling with lb 0 / 0.4", draws(f1, pI))
cat("\n#### 2. a ~ 1 + w\n")
try_("bound refused?", ents(fw, set_prior("normal(1, 1)", class = "Intercept",
                                          nlpar = "a", lb = 0)))
try_("no bound: a coefficient row?",
     ents(fw, set_prior("normal(5, 0.5)", class = "Intercept", nlpar = "a")))
try_("a ~ 1 + (1 | g), lb 0: class row?",
     ents(fg, set_prior("normal(1, 1)", class = "Intercept", nlpar = "a",
                        lb = 0)))
try_("a ~ 1 + (1 | g), draws", draws(fg, pI))
cat("\n#### 3. both spellings\n")
pA <- set_prior("normal(5, 0.1)", class = "Intercept", nlpar = "a")
pB <- set_prior("normal(0, 0.1)", class = "b", nlpar = "a")
try_("Intercept-spelling then class row", ents(f1, pA + pB))
try_("class row then Intercept-spelling", ents(f1, pB + pA))
try_("frm() MAP, Intercept-spelling then class row", {
  fit <- frm(f1, data = d, prior = pA + pB)
  fit$estimates$beta
})
try_("frm() MAP, class row then Intercept-spelling", {
  fit <- frm(f1, data = d, prior = pB + pA)
  fit$estimates$beta
})
try_("same spelling twice (duplicate)",
     {frm_sample(f1, data = d, family = gaussian(), prior = pA + pA,
                 chains = 1, iter = 50, refresh = 0); "sampled"})
try_("frm_sample, both spellings, posterior mean of b_a_Intercept", {
  ds <- frm_sample(f1, data = d, family = gaussian(), prior = pA + pB +
                     set_prior("normal(1, 1)", class = "b", nlpar = "b"),
                   chains = 1, iter = 600, seed = 3, refresh = 0)
  colMeans(as.matrix(ds)[, grep("^b_a", colnames(as.matrix(ds))),
                         drop = FALSE])
})
try_("fit's own prior (Intercept spelling) + call's class row", {
  fit <- frm(f1, data = d, prior = pA)
  fs <- frm_sample(fit, prior = pB, chains = 1, iter = 400, seed = 3,
                   refresh = 0)
  colMeans(as.matrix(fs)[, grep("^b_a", colnames(as.matrix(fs))),
                         drop = FALSE])
})
cat("\n#### 3b. guard absent: class b, nlpar a\n")
try_("class b nlpar a lb 0", ents(f1, set_prior("normal(1, 1)", class = "b",
                                               nlpar = "a", lb = 0)))
try_("draws, class b nlpar spelling",
     draws(f1, set_prior("normal(1, 1)", class = "b", nlpar = "a", lb = 0) +
             set_prior("normal(1, 1)", class = "b", nlpar = "b", lb = 0.4)))
cat("\n#### 4. frm()\n")
for (f in list(f1, fw)) {
  try_("frm logLik and estimates, Intercept nlpar lb 0", {
    fit <- frm(f, data = d, prior = set_prior("normal(1, 1)",
                                              class = "Intercept",
                                              nlpar = "a", lb = 0))
    c(logLik = sprintf("%.15g", as.numeric(logLik(fit))),
      sprintf("%.15g", fit$obj$env$last.par.best))
  })
}
cat("DONE\n")
