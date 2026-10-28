# Reviewer: autocor defaults under sample_prior. brms 2.23.0 gives an
# ar (cov = TRUE) or cosy parameter with a default prior no prior
# statement: no prior_ column under "yes", a refusal under "only".
#   Rscript dev/review-priordraws-acdefault.R [lane|rv]
arm <- commandArgs(TRUE)[1]
if (is.na(arm)) arm <- "lane"
.libPaths(c(if (arm == "rv") "C:/Users/adf44/source/r/wt-priordraws-rv-lib",
            "C:/Users/adf44/source/r/wt-priordraws-lib",
            "C:/Users/adf44/source/r/rellib-r14",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-priordraws/dev/review-priordraws-stancache")
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("lib:", find.package("frmtmb"), find.package("frmtmb.sample"), "\n")
set.seed(1)
d <- data.frame(x = rnorm(60), h = factor(rep(1:5, each = 12)))
d$y <- 1 + 0.5 * d$x + rnorm(60)
d$t <- rep(1:12, 5)
pr <- set_prior("normal(0, 1)", class = "b")
for (f in list(y ~ x + ar(time = t, gr = h, cov = TRUE),
               y ~ x + cosy(time = t, gr = h),
               y ~ x + ma(time = t, gr = h, cov = TRUE))) {
  cat("\n=====", deparse(f), "=====\n")
  vp <- tryCatch(validate_prior(pr, bf(f), data = d, family = gaussian()),
                 error = function(e) conditionMessage(e))
  if (is.character(vp)) cat("validate_prior:", vp, "\n") else
    print(as.data.frame(vp)[as.data.frame(vp)$class %in% c("ar", "ma", "cosy"), ])
  ds <- tryCatch(suppressWarnings(suppressMessages(frm_sample(
    bf(f), data = d, family = gaussian(), prior = pr, sample_prior = "yes",
    chains = 1, iter = 400, seed = 1, refresh = 0))),
    error = function(e) e)
  cat("yes:", if (inherits(ds, "error")) conditionMessage(ds) else
    paste(names(prior_draws(ds)), collapse = " "), "\n")
  if (!inherits(ds, "error")) {
    pa <- prior_draws(ds)
    for (nm in intersect(c("ar", "ma", "cosy"), names(pa))) {
      cat(" ", nm, "range", range(pa[[nm]]), "\n")
    }
  }
  do <- tryCatch(suppressWarnings(suppressMessages(frm_sample(
    bf(f), data = d, family = gaussian(), prior = pr, sample_prior = "only",
    chains = 1, iter = 400, seed = 1, refresh = 0))),
    error = function(e) e)
  cat("only:", if (inherits(do, "error")) conditionMessage(do) else
    paste("sampled", ndraws(do)), "\n")
}
# the mixture case with brms's own families
cat("\n===== mixture brms stancode only =====\n")
mp <- brms::set_prior("normal(0, 1)", class = "b", dpar = "mu1") +
  brms::set_prior("normal(0, 1)", class = "b", dpar = "mu2") +
  brms::set_prior("normal(0, 1)", class = "Intercept", dpar = "mu1") +
  brms::set_prior("normal(1, 1)", class = "Intercept", dpar = "mu2")
cat(tryCatch({brms::stancode(brms::bf(y ~ x), data = d,
                             family = brms::mixture("gaussian", "gaussian"),
                             prior = mp, sample_prior = "only"); "OK"},
             error = function(e) conditionMessage(e)), "\n")
cat("\n===== brms stancode prior_ statements, sample_prior = yes =====\n")
bp <- brms::set_prior("normal(0, 1)", class = "b")
for (f in list(y ~ x + ar(time = t, gr = h, cov = TRUE),
               y ~ x + cosy(time = t, gr = h),
               y ~ x + ma(time = t, gr = h, cov = TRUE),
               y ~ x + ar(time = t, gr = h))) {
  sc <- brms::stancode(brms::bf(f), data = d, family = brms::brmsfamily("gaussian"),
                       prior = bp, sample_prior = "yes")
  ln <- grep("prior_", strsplit(sc, "\n")[[1]], value = TRUE)
  cat(deparse(f), "\n", paste(" ", trimws(ln), collapse = "\n"), "\n")
  dl <- grep("^ *(vector|real).*(ar|ma|cosy);", strsplit(sc, "\n")[[1]], value = TRUE)
  cat("  decl:", trimws(dl), "\n")
}
