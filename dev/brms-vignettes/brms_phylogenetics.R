# Translation of brms's "Estimating Phylogenetic Multilevel Models with
# brms" (brms 2.23.0, doc/brms_phylogenetics.Rmd) onto the frmtmb
# surface.
#
# Read the comments, not just the code: every `<EDGE>` marker on a bv()
# call is a place where the brms line did not carry over verbatim.
#
# The script covers both inference paths. PATH 1 is `frm()` and the
# frequentist post-processing surface; PATH 2 is `frm_sample()` and the
# draws surface, which is the closer analogue of what the vignette
# literally prints. Every label carries an "ML: " or "SAMPLE: " prefix.
#
# Sizes. PATH 1 keeps the vignette's own sizes: 200 species for the
# simple, meta-analytic and count models, 1000 rows for the
# repeated-measurement model, every fit under 5 s. PATH 2 uses the first
# 40 species (and the 200 repeated-measurement rows that belong to
# them), because a 200-level phylogenetic covariance block does not
# finish a chain in the time budget. The measurement is in a BEHAVIOR
# row at the head of PATH 2.

source(file.path(Sys.getenv(
  "BV_DIR",
  unset = "C:/Users/adf44/source/r/frmtmb-wt-brmsvig/dev/brms-vignettes"
), "_harness.R"))
bv_load()
bv_init("brms_phylogenetics")

options(mc.cores = 1)
set.seed(1234)

## ---- The data ------------------------------------------------------
#
# The vignette reads a nexus tree and three tables from
# paul-buerkner.github.io. A download failure must not turn the whole
# audit into one error, so the fetch falls back to a simulated tree and
# simulated tables of the same shape. `bv_sim` records which of the two
# the numbers below come from.

bv_url <- "https://paul-buerkner.github.io/data/"
bv_dir <- file.path(tempdir(), "bv-phylo")
dir.create(bv_dir, showWarnings = FALSE, recursive = TRUE)

bv_fetch <- function(f) {
  dest <- file.path(bv_dir, f)
  if (file.exists(dest)) return(TRUE)
  ok <- tryCatch({
    utils::download.file(paste0(bv_url, f), dest, quiet = TRUE, mode = "wb")
    file.exists(dest) && file.size(dest) > 0
  }, error = function(e) FALSE, warning = function(w) FALSE)
  isTRUE(ok)
}

bv_sim <- !all(vapply(
  c("phylo.nex", "data_simple.txt", "data_repeat.txt",
    "data_effect.txt", "data_pois.txt"),
  bv_fetch, logical(1)))

if (bv_sim) {
  # SUBSTITUTED DATA. The download failed, so a 60-species coalescent
  # tree and four tables of the vignette's shape stand in. The edge
  # labels below still hold; only the printed numbers change.
  message("NOTE: download failed; simulated tree and data substituted")
  phylo <- ape::rcoal(60, tip.label = paste0("sp_", 1:60))
  A <- ape::vcv.phylo(phylo)
  A <- A / mean(diag(A))
  sp <- rownames(A)
  re <- as.numeric(t(chol(A)) %*% rnorm(60)) * 14
  names(re) <- sp
  cof <- rnorm(60, 10, 4)
  data_simple <- data.frame(
    phen = 38 + 5 * cof + re + rnorm(60, 0, 9),
    cofactor = cof, phylo = sp, stringsAsFactors = FALSE)
  data_repeat <- data.frame(
    phylo = rep(sp, each = 5), species = rep(sp, each = 5),
    cofactor = rnorm(300, 10, 4), stringsAsFactors = FALSE)
  data_repeat$phen <- 36 + 5 * data_repeat$cofactor +
    re[data_repeat$phylo] + rep(rnorm(60, 0, 5), each = 5) +
    rnorm(300, 0, 8)
  data_fisher <- data.frame(
    Zr = 0.16 + re / 200 + rnorm(60, 0, 0.07),
    N = sample(10:50, 60, TRUE), phylo = sp, stringsAsFactors = FALSE)
  data_pois <- data.frame(
    cofactor = cof, phylo = sp, stringsAsFactors = FALSE)
  data_pois$phen_pois <- stats::rpois(
    60, exp(-2.1 + 0.25 * cof + re / 70 + rnorm(60, 0, 0.2)))
} else {
  phylo <- ape::read.nexus(file.path(bv_dir, "phylo.nex"))
  data_simple <- utils::read.table(file.path(bv_dir, "data_simple.txt"),
                                   header = TRUE)
  data_repeat <- utils::read.table(file.path(bv_dir, "data_repeat.txt"),
                                   header = TRUE)
  data_fisher <- utils::read.table(file.path(bv_dir, "data_effect.txt"),
                                   header = TRUE)
  data_pois <- utils::read.table(file.path(bv_dir, "data_pois.txt"),
                                 header = TRUE)
  A <- ape::vcv.phylo(phylo)
}

bv("data", "head(data_simple)", print(utils::head(data_simple)),
   NA_character_, "")

bv("data", "A <- ape::vcv.phylo(phylo)", {
  cat("dim(A):", dim(A), " simulated:", bv_sim, "\n")
  dim(A)
}, NA_character_, "")

## ============ PATH 1: ML / Laplace (frm) ============

## ---- A simple phylogenetic model ------------------------------------
#
# brms carries four priors here and says in the text that they are not
# needed for convergence, only for sampling speed. That makes them
# MCMC-only, so they are dropped. The prior block is measured on its own
# further down, under "The prior surface".
#
# `(1 | gr(phylo, cov = A))` with `data2 = list(A = A)` is accepted
# verbatim. `family = gaussian()` is now a default as well (the v0.34
# audit recorded the missing default as FN-1; it is fixed), but the
# house style states the family, so it is stated.

model_simple <- bv("model", "ML: model_simple", {
  frm(bf(phen ~ cofactor + (1 | gr(phylo, cov = A))),
      data = data_simple, family = gaussian(),
      data2 = list(A = A))
}, NA_character_, "")

# brms prints sigma on the RESPONSE scale (9.24). frmtmb prints the
# sigma linear predictor on its LINK scale, so the same quantity reads
# 2.22 = log(9.2). A reader comparing the two summaries meets this
# first.
bv("post", "ML: summary(model_simple)", summary(model_simple), "BEHAVIOR",
   "sigma is printed on the log link scale (2.22) where brms reports the response-scale residual SD (9.24)")

# brms's plot() draws trace and density plots per parameter and `N`
# picks how many go on a page. frmtmb has no chain to trace, so the
# panels are residual diagnostics and `N` is absorbed by dots.
bv("post", "ML: plot(model_simple, N = 2, ask = FALSE)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  plot(model_simple, ask = FALSE)
}, "BEHAVIOR", "N = is silently absorbed by ...; the panels are residual and QQ diagnostics, not traces")

bv("post", "ML: plot(conditional_effects(model_simple), points = TRUE)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  plot(conditional_effects(model_simple), points = TRUE)
}, NA_character_, "")

# The phylogenetic signal. The brms string ports verbatim, including
# `class = NULL` and the `sd_phylo__Intercept` / `sigma` names, and
# lands on the vignette's lambda = 0.7.
hyp1 <- bv("post", "ML: hypothesis(model_simple, hyp, class = NULL)", {
  hypothesis(model_simple,
             "sd_phylo__Intercept^2 / (sd_phylo__Intercept^2 + sigma^2) = 0",
             class = NULL)
}, NA_character_, "")

# brms's plot(hyp) draws the posterior density of the derived quantity.
# There is no posterior here, so frmtmb draws the Wald normal curve
# implied by the delta-method standard error.
bv("post", "ML: plot(hyp)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  plot(hyp1)
}, "BEHAVIOR", "a Wald normal curve replaces brms's posterior density; no Evid.Ratio or Post.Prob column")

# Not the vignette's line, but the other half of brms's hypothesis
# surface. The v0.34 audit recorded the one-sided form as refused on a
# point fit (FN-5); it is accepted now, on the ML path as well as on the
# draws path, and the print says which rows are one-sided.
bv("post", "ML: hypothesis(model_simple, 'cofactor > 0')", {
  hypothesis(model_simple, "cofactor > 0")
}, "BEHAVIOR", "brms's one-sided form runs and reports a one-sided p and a half-open interval where brms prints an evidence ratio and a posterior probability")

## ---- A phylogenetic model with repeated measurements ----------------

data_repeat$spec_mean_cf <-
  with(data_repeat, sapply(split(cofactor, phylo), mean)[phylo])

bv("data", "head(data_repeat)", print(utils::head(data_repeat)),
   NA_character_, "")

# Two grouping factors over the same species, one of them covariance
# structured. brms adds sample_prior = TRUE, which is MCMC-only.
model_repeat1 <- bv("model", "ML: model_repeat1", {
  frm(bf(phen ~ spec_mean_cf + (1 | gr(phylo, cov = A)) + (1 | species)),
      data = data_repeat, family = gaussian(),
      data2 = list(A = A))
}, NA_character_, "")

bv("post", "ML: summary(model_repeat1)", summary(model_repeat1), "BEHAVIOR",
   "sigma on the log link scale again")

hyp2 <- bv("post", "ML: hypothesis(model_repeat1, hyp, class = NULL)", {
  hypothesis(model_repeat1, paste(
    "sd_phylo__Intercept^2 /",
    "(sd_phylo__Intercept^2 + sd_species__Intercept^2 + sigma^2) = 0"),
    class = NULL)
}, NA_character_, "")

bv("post", "ML: plot(hyp) [repeat1]", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  plot(hyp2)
}, "BEHAVIOR", "Wald curve, as above")

data_repeat$within_spec_cf <-
  data_repeat$cofactor - data_repeat$spec_mean_cf

# The vignette's update() line survives with only the MCMC arguments
# removed: `newdata` is accepted as an alias for `data`, and the
# one-sided `~ . + within_spec_cf` is expanded against the stored
# formula. The v0.34 audit recorded both as failures (FN-9); both are
# fixed.
model_repeat2 <- bv("model", "ML: model_repeat2 (update)", {
  update(model_repeat1, formula = ~ . + within_spec_cf,
         newdata = data_repeat)
}, NA_character_, "")

bv("post", "ML: summary(model_repeat2)", summary(model_repeat2), "BEHAVIOR",
   "sigma on the log link scale again")

bv("post", "ML: hypothesis(model_repeat2, hyp, class = NULL)", {
  hypothesis(model_repeat2, paste(
    "sd_phylo__Intercept^2 /",
    "(sd_phylo__Intercept^2 + sd_species__Intercept^2 + sigma^2) = 0"),
    class = NULL)
}, NA_character_, "")

## ---- A phylogenetic meta-analysis ------------------------------------

data_fisher$obs <- seq_len(nrow(data_fisher))
bv("data", "head(data_fisher)", print(utils::head(data_fisher)),
   NA_character_, "")

# `se()` with a known sampling SD, plus an observation-level grouping
# factor for the residual variance. Ports verbatim; adapt_delta and the
# two priors are MCMC-only.
model_fisher <- bv("model", "ML: model_fisher", {
  frm(bf(Zr | se(sqrt(1 / (N - 3))) ~ 1 + (1 | gr(phylo, cov = A)) +
           (1 | obs)),
      data = data_fisher, family = gaussian(),
      data2 = list(A = A))
}, NA_character_, "")

# The meta-analytic mean lands on the vignette's 0.16. The two variance
# components do not: brms reports sd(phylo) = 0.05 and sd(obs) = 0.07,
# frmtmb 0.01 and 0.07. The phylogenetic component sits at the boundary,
# where brms's half-t prior holds it away from zero and maximum
# likelihood does not.
bv("post", "ML: summary(model_fisher)", summary(model_fisher), "BEHAVIOR",
   "the boundary variance component sd(phylo) collapses toward zero without brms's half-t prior; the meta-analytic mean 0.16 agrees")

bv("post", "ML: plot(model_fisher)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  plot(model_fisher, ask = FALSE)
}, "BEHAVIOR", "diagnostics, not traces")

## ---- A phylogenetic count-data model ---------------------------------

data_pois$obs <- seq_len(nrow(data_pois))
bv("data", "head(data_pois)", print(utils::head(data_pois)),
   NA_character_, "")

model_pois <- bv("model", "ML: model_pois", {
  frm(bf(phen_pois ~ cofactor + (1 | gr(phylo, cov = A)) + (1 | obs)),
      data = data_pois, family = poisson("log"),
      data2 = list(A = A))
}, NA_character_, "")

bv("post", "ML: summary(model_pois)", summary(model_pois), NA_character_, "")

bv("post", "ML: plot(conditional_effects(model_pois), points = TRUE)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  plot(conditional_effects(model_pois), points = TRUE)
}, NA_character_, "")

model_normal <- bv("model", "ML: model_normal", {
  frm(bf(phen_pois ~ cofactor + (1 | gr(phylo, cov = A))),
      data = data_pois, family = gaussian(),
      data2 = list(A = A))
}, NA_character_, "")

bv("post", "ML: summary(model_normal)", summary(model_normal), "BEHAVIOR",
   "sigma on the log link scale")

bv("post", "ML: pp_check(model_pois)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  print(pp_check(model_pois))
}, NA_character_, "")

bv("post", "ML: pp_check(model_normal)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  print(pp_check(model_normal))
}, NA_character_, "")

# loo() is a frmtmb_draws method. A frmtmb_fit has no draws to build a
# pointwise log-likelihood matrix from.
bv("post", "ML: loo(model_pois, model_normal)", loo(model_pois, model_normal),
   "MISSING", "loo() dispatches only on frmtmb_draws; the frequentist route is AIC()/BIC()")

bv("post", "ML: AIC(model_pois, model_normal) [the substitute]", {
  print(AIC(model_pois, model_normal))
}, "SPELLING", "the loo comparison becomes AIC; no elpd, no Pareto k, no per-observation diagnostic")

## ---- The prior surface -----------------------------------------------
#
# The four priors of `model_simple` are dropped above because the
# vignette says they only help the sampler. They are still measured
# here, because the frmtmb prior surface diverges from brms in spelling
# and in meaning, and a sibling lane is working on it.
#
# In brms a prior is a prior. In frmtmb `set_prior()` is a penalty on
# the likelihood, so the fit below is penalized maximum likelihood and
# not the same estimator. That is by design, not a defect, but it means
# a reader cannot port the prior block and expect brms's numbers.

bv("post", "ML: prior(normal(0, 10), 'b') [brms spelling]", {
  prior(normal(0, 10), "b")
}, "MISSING", "brms's NSE prior() is not exported; the frmtmb spelling is set_prior('normal(0, 10)', class = 'b')")

bv("post", "ML: set_prior(..., class = 'sigma')", {
  set_prior("student_t(3, 0, 20)", class = "sigma")
}, "REFUSAL", "class 'sigma' is refused; the message lists b/Intercept/sd/cor/theta but does not name the replacement, class = 'Intercept' with dpar = 'sigma'")

bv_priors <- bv("post", "ML: the four-prior block, frmtmb spelling", {
  set_prior("normal(0, 10)", class = "b") +
    set_prior("normal(0, 50)", class = "Intercept") +
    set_prior("student_t(3, 0, 20)", class = "sd") +
    set_prior("student_t(3, 0, 20)", class = "Intercept", dpar = "sigma")
}, "SPELLING", "prior(dist, class) becomes set_prior('dist', class = ...) and the sigma prior moves to dpar = 'sigma'")

bv("post", "ML: frm(..., priors = ) [the penalized fit]", {
  summary(frm(bf(phen ~ cofactor + (1 | gr(phylo, cov = A))),
              data = data_simple, family = gaussian(),
              data2 = list(A = A), priors = bv_priors))
}, "BEHAVIOR", "priors = penalizes the likelihood rather than defining a posterior, and the summary does not say the fit is penalized")

bv("post", "ML: prior_summary(model_simple)", prior_summary(model_simple),
   "BEHAVIOR", "brms lists the default priors it chose; frmtmb answers 'No priors were set (plain maximum likelihood)', which is true but is not the same report")

bv("post", "ML: get_prior(...)", {
  get_prior(bf(phen ~ cofactor + (1 | gr(phylo, cov = A))),
            data = data_simple, family = gaussian(), data2 = list(A = A))
}, "BEHAVIOR", "the brms columns are there but every row reads (flat) and a 'theta' class appears that brms has no name for")

## ============ PATH 2: sampling (frm_sample) ============
#
# The vignette's printed output is posterior output, so this path is the
# closer analogue. Two things are different from PATH 1.
#
# Size. A 200-level phylogenetic covariance block does not finish a
# chain inside the budget: at 60 species one 300-iteration chain takes
# 18 s and diverges, at 200 species it does not finish in 120 s. PATH 2
# therefore uses the first 40 species and the 200 repeated-measurement
# rows that belong to them. The vignette's sizes are 200 species and
# 1000 rows.
#
# Route. `model_simple`, `model_repeat1` and `model_fisher` carry priors
# in brms, so they take the FORMULA route, where frm_sample() supplies
# brms 2.23's own default priors. `model_pois` and `model_normal` carry
# no prior in brms, so they take the FIT route from a 40-species frm()
# fit, which is what a reader who already has a fit would write.

bv_sp <- utils::head(unique(as.character(data_simple$phylo)), 40)
A_s <- A[bv_sp, bv_sp]
s_simple <- data_simple[as.character(data_simple$phylo) %in% bv_sp, ]
s_repeat <- data_repeat[as.character(data_repeat$phylo) %in% bv_sp, ]
s_fisher <- data_fisher[as.character(data_fisher$phylo) %in% bv_sp, ]
s_fisher$obs <- seq_len(nrow(s_fisher))
s_pois <- data_pois[as.character(data_pois$phylo) %in% bv_sp, ]
s_pois$obs <- seq_len(nrow(s_pois))

bv("data", "SAMPLE: the shrunken data", {
  cat("species:", length(bv_sp), " simple:", nrow(s_simple),
      " repeat:", nrow(s_repeat), " fisher:", nrow(s_fisher),
      " pois:", nrow(s_pois), "\n")
  length(bv_sp)
}, "BEHAVIOR", "the sample path needs 40 species where the vignette uses 200; a 200-level gr(cov = ) block does not finish a chain in 120 s")

bv_draws <- function(...) {
  frm_sample(..., chains = 1, iter = 400, warmup = 200, seed = 1,
             cores = 1, refresh = 0)
}

# frm_sample() prints the brms default priors it supplied. That is the
# closest the two packages come on the prior surface, and it is only
# reachable from the formula route.
sample_simple <- bv("model", "SAMPLE: model_simple", {
  bv_draws(bf(phen ~ cofactor + (1 | gr(phylo, cov = A))),
           data = s_simple, family = gaussian(), data2 = list(A = A_s))
}, "SPELLING", "frm_sample() replaces brm()'s sampler arguments and supplies brms 2.23 default priors in place of the vignette's prior block")

# brms groups the summary into Multilevel Hyperparameters, Regression
# Coefficients and Further Distributional Parameters, and names the
# variance component sd(Intercept). frmtmb prints one flat matrix and
# names the same component theta_1, on its internal scale.
bv("post", "SAMPLE: summary(model_simple)", summary(sample_simple),
   "BEHAVIOR", "one flat matrix instead of brms's three named blocks; the phylo SD appears as the internal theta_1 and sigma as sigma_Intercept on the log scale")

# The draws surface has the trace and density plots that PATH 1 cannot
# draw, under the bayesplot name rather than under plot().
bv("post", "SAMPLE: plot(model_simple, N = 2) -> mcmc_plot()", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  print(mcmc_plot(sample_simple))
}, "SPELLING", "brms's plot() becomes mcmc_plot(); N = has no counterpart")

bv("post", "SAMPLE: plot(conditional_effects(model_simple), points = TRUE)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  plot(conditional_effects(sample_simple), points = TRUE)
}, NA_character_, "")

# The phylogenetic signal, on the posterior this time. The estimate is a
# posterior mean of the derived quantity, which is what brms reports,
# but the columns are still frequentist: no Evid.Ratio, no Post.Prob, no
# Star.
hyp1s <- bv("post", "SAMPLE: hypothesis(model_simple, hyp, class = NULL)", {
  hypothesis(sample_simple,
             "sd_phylo__Intercept^2 / (sd_phylo__Intercept^2 + sigma^2) = 0",
             class = NULL)
}, "BEHAVIOR", "method = posterior, so the estimate is a posterior mean, but the columns are z and p where brms prints Evid.Ratio and Post.Prob")

bv("post", "SAMPLE: plot(hyp)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  plot(hyp1s)
}, NA_character_, "")

# The one-sided form on the draws path. Both paths accept it now (the
# v0.34 audit's FN-5 is fixed), and the only difference is the method
# line: posterior here, wald there.
bv("post", "SAMPLE: hypothesis(model_simple, 'cofactor > 0')", {
  hypothesis(sample_simple, "cofactor > 0")
}, "BEHAVIOR", "the one-sided form runs on the posterior, but the report is still a one-sided p rather than brms's evidence ratio and posterior probability")

sample_repeat1 <- bv("model", "SAMPLE: model_repeat1", {
  bv_draws(bf(phen ~ spec_mean_cf + (1 | gr(phylo, cov = A)) +
                (1 | species)),
           data = s_repeat, family = gaussian(), data2 = list(A = A_s))
}, "SPELLING", "frm_sample() with the default priors; sample_prior = TRUE has no counterpart")

bv("post", "SAMPLE: summary(model_repeat1)", summary(sample_repeat1),
   "BEHAVIOR", "flat matrix, internal theta names, as above")

bv("post", "SAMPLE: hypothesis(model_repeat1, hyp, class = NULL)", {
  hypothesis(sample_repeat1, paste(
    "sd_phylo__Intercept^2 /",
    "(sd_phylo__Intercept^2 + sd_species__Intercept^2 + sigma^2) = 0"),
    class = NULL)
}, "BEHAVIOR", "posterior method, frequentist columns, as above")

# update() has no draws method. The stored call is replayed without the
# sampler arguments, so it dies inside frm()'s own argument list.
bv("model", "SAMPLE: model_repeat2 (update)", {
  update(sample_repeat1, formula = ~ . + within_spec_cf,
         newdata = s_repeat)
}, "MISSING", "the draws surface has no update() method, so newdata = lands as an unused argument of frm(); the ML path accepts the same line")

sample_repeat2 <- bv("post", "SAMPLE: model_repeat2 workaround (re-sample)", {
  bv_draws(bf(phen ~ spec_mean_cf + within_spec_cf +
                (1 | gr(phylo, cov = A)) + (1 | species)),
           data = s_repeat, family = gaussian(), data2 = list(A = A_s))
}, "SPELLING", "the whole model is re-sampled from the extended formula because update() has no draws method")

bv("post", "SAMPLE: hypothesis(model_repeat2, hyp, class = NULL)", {
  hypothesis(sample_repeat2, paste(
    "sd_phylo__Intercept^2 /",
    "(sd_phylo__Intercept^2 + sd_species__Intercept^2 + sigma^2) = 0"),
    class = NULL)
}, "BEHAVIOR", "posterior method, frequentist columns, as above")

sample_fisher <- bv("model", "SAMPLE: model_fisher", {
  bv_draws(bf(Zr | se(sqrt(1 / (N - 3))) ~ 1 + (1 | gr(phylo, cov = A)) +
                (1 | obs)),
           data = s_fisher, family = gaussian(), data2 = list(A = A_s))
}, "SPELLING", "frm_sample() with the default priors in place of the vignette's two")

bv("post", "SAMPLE: summary(model_fisher)", summary(sample_fisher),
   "BEHAVIOR", "flat matrix, internal theta names")

bv("post", "SAMPLE: plot(model_fisher) -> mcmc_plot()", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  print(mcmc_plot(sample_fisher))
}, "SPELLING", "plot() becomes mcmc_plot()")

# The FIT route. brms gives these two models no prior, so a fit, which
# carries flat priors, is the honest translation. It does not work here,
# and the reason is worth the row: both fits are at a variance boundary,
# mode initialization starts the chain there, and the chain never
# leaves. Every draw is the same point, the reported sd is 0 and Rhat is
# 1.5 to 2.1. frm_sample() warns about exactly this and names both
# remedies, init = "random" and priors =.
sample_pois_fit <- bv("model", "SAMPLE: model_pois [fit route]", {
  bv_draws(frm(bf(phen_pois ~ cofactor + (1 | gr(phylo, cov = A)) +
                    (1 | obs)),
               data = s_pois, family = poisson("log"),
               data2 = list(A = A_s)))
}, "BEHAVIOR", "the fit route starts at a boundary ML mode with flat variance priors and the chain sticks there; the warning names init = 'random' and priors = as the remedies")

bv("post", "SAMPLE: summary(model_pois) [fit route]", summary(sample_pois_fit),
   "BEHAVIOR", "a stuck chain: sd = 0 on every parameter and Rhat above 1.5, so none of the vignette's numbers are reachable this way")

# The formula route on the same model. The default priors hold the
# variance parameters off the boundary and the chain moves, so this is
# the object the rest of the section uses.
sample_pois <- bv("model", "SAMPLE: model_pois [formula route]", {
  bv_draws(bf(phen_pois ~ cofactor + (1 | gr(phylo, cov = A)) + (1 | obs)),
           data = s_pois, family = poisson("log"), data2 = list(A = A_s))
}, "SPELLING", "the formula route replaces the fit route because brms's implicit priors, not flat ones, are what keep this chain off the boundary")

bv("post", "SAMPLE: summary(model_pois)", summary(sample_pois),
   "BEHAVIOR", "flat matrix and internal theta names; at 400 iterations Rhat stays above 1.1, so the printed numbers are not the vignette's")

bv("post", "SAMPLE: plot(conditional_effects(model_pois), points = TRUE)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  plot(conditional_effects(sample_pois), points = TRUE)
}, NA_character_, "")

sample_normal <- bv("model", "SAMPLE: model_normal", {
  bv_draws(bf(phen_pois ~ cofactor + (1 | gr(phylo, cov = A))),
           data = s_pois, family = gaussian(), data2 = list(A = A_s))
}, "SPELLING", "the formula route again, for the same reason")

bv("post", "SAMPLE: summary(model_normal)", summary(sample_normal),
   "BEHAVIOR", "flat matrix, internal theta names")

bv("post", "SAMPLE: pp_check(model_pois)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  print(pp_check(sample_pois))
}, NA_character_, "")

bv("post", "SAMPLE: pp_check(model_normal)", {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  print(pp_check(sample_normal))
}, NA_character_, "")

# The vignette's own loo() line, with two models in one call. loo() on
# draws takes one model: the second lands in `ndraws`.
bv("post", "SAMPLE: loo(model_pois, model_normal)", {
  loo(sample_pois, sample_normal)
}, "MISSING", "the two-model form errors with \"'list' object cannot be coerced to type 'integer'\", an internal message that neither refuses nor names loo_compare()")

# The vignette's conclusion, that the Poisson model fits better, does
# survive the rewrite.
bv("post", "SAMPLE: loo_compare(loo(a), loo(b)) [the substitute]", {
  loo_compare(loo(sample_pois), loo(sample_normal))
}, "SPELLING", "the comparison is spelled loo_compare(loo(x), loo(y)); the elpd difference still favors the Poisson model, but a 200-draw chain flags most Pareto k values as bad")

bv("post", "SAMPLE: bayes_R2(model_pois)", bayes_R2(sample_pois),
   NA_character_, "")

bv_done()
