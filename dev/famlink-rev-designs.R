## Reviewer check for lane wt-famlink, priorities 3 and 4: false alarms
## and silent changes beyond `y ~ x`.
##
## One list of constructions (family spelling x design) is run in one
## arm: "lane", "base" (frmtmb fits) or "brms" (brms::make_standata,
## which runs brms's family, link, mixture and response validation
## without compiling). Seed 20260916, n = 150, 15 groups.
##
## Usage: Rscript dev/famlink-rev-designs.R <lane|base|brms>
## Writes dev/famlink-rev-designs-<arm>.rds: per construction, status
## (ok / error), message, whether a bernoulli message fired and how many
## times, and for frmtmb arms logLik, estimates, SEs and fitted values.
ARM <- commandArgs(trailingOnly = TRUE)[1]
if (identical(ARM, "brms")) {
  .libPaths(c("C:/Users/adf44/source/r/pinlib",
              "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
  suppressMessages(library(brms))
} else {
  source("dev/famlink-rev-common.R")
}
set.seed(20260916)
n <- 150
d <- data.frame(x = rnorm(n), g = factor(rep(1:15, each = 10)))
re <- rnorm(15, 0, 0.3)[d$g]
d$yg <- 1 + 0.5 * d$x + re + rnorm(n)
d$ypos <- exp(0.5 + 0.3 * d$x + re + rnorm(n, 0, 0.4))
d$ycnt <- rpois(n, exp(1 + 0.3 * d$x + re))
d$y01 <- rbinom(n, 1, plogis(0.2 + 0.7 * d$x + re))
d$tr <- sample(3:10, n, TRUE)
d$ybin <- rbinom(n, d$tr, plogis(0.2 + 0.5 * d$x))
d$one <- 1L
d$yprop <- plogis(qlogis(0.4) + 0.5 * d$x + rnorm(n, 0, 0.5))
d$yord <- cut(d$yg, c(-Inf, 0, 1, 2, Inf), labels = FALSE)
d$yord_f <- factor(d$yord, ordered = TRUE)
d$y2cat <- ifelse(d$yg > 1, 2L, 1L)
d$y2cat_f <- factor(ifelse(d$yg > 1, "hi", "lo"), ordered = TRUE)
d$ycat <- factor(c("a", "b", "c")[d$yord %% 3 + 1])
d$cen <- ifelse(d$ypos > 3, "right", "none")
d$ycen <- pmin(d$ypos, 3)
d$yzi <- ifelse(runif(n) < 0.3, 0L, d$ycnt)
d$yangle <- atan2(sin(0.4 * d$x + rnorm(n, 0, 0.5)), cos(0.4 * d$x + rnorm(n, 0, 0.5)))
d$y1 <- d$yg; d$y2 <- d$ycnt
d$ymix <- c(rnorm(n / 2, -1), rnorm(n / 2, 2))
d$ymixc <- c(rpois(n / 2, 2), rpois(n / 2, 9))

F <- if (identical(ARM, "brms")) brms::bf else frmtmb::bf
mvb <- if (identical(ARM, "brms")) brms::mvbf else frmtmb::mvbf
cases <- list(
  # family as string, and brms's spellings
  str_gaussian = list(y ~ x, "gaussian", "yg"),
  str_normal = list(y ~ x, "normal", "yg"),
  str_gamma_lower = list(y ~ x, "gamma", "ypos"),
  str_Gamma = list(y ~ x, "Gamma", "ypos"),
  str_zi_poisson = list(y ~ x, "zi_poisson", "yzi"),
  str_hu_poisson = list(y ~ x, "hu_poisson", "yzi"),
  str_cumulative = list(y ~ x, "cumulative", "yord_f"),
  str_Bernoulli_case = list(y ~ x, "Bernoulli", "y01"),
  # c(name, link)
  pair_poisson_sqrt = list(y ~ x, c("poisson", "sqrt"), "ycnt"),
  pair_bernoulli_probit = list(y ~ x, c("bernoulli", "probit"), "y01"),
  pair_gamma_inverse = list(y ~ x, c("gamma", "inverse"), "ypos"),
  pair_cumulative_probit = list(y ~ x, c("cumulative", "probit"), "yord_f"),
  # stats family objects and their default spellings
  stats_Gamma_default = list(y ~ x, quote(stats::Gamma()), "ypos"),
  stats_Gamma_identity = list(y ~ x, quote(stats::Gamma("identity")), "ypos"),
  stats_invgauss_default = list(y ~ x, quote(stats::inverse.gaussian()), "ypos"),
  stats_invgauss_log = list(y ~ x, quote(stats::inverse.gaussian("log")), "ypos"),
  stats_gaussian_log = list(y ~ x, quote(stats::gaussian("log")), "ypos"),
  stats_gaussian_inverse = list(y ~ x, quote(stats::gaussian("inverse")), "ypos"),
  stats_poisson_sqrt = list(y ~ x, quote(stats::poisson("sqrt")), "ycnt"),
  stats_poisson_identity = list(y ~ x, quote(stats::poisson("identity")), "ycnt"),
  stats_binomial_probit = list(y | trials(tr) ~ x, quote(stats::binomial("probit")), "ybin"),
  stats_binomial_cauchit = list(y | trials(tr) ~ x, quote(stats::binomial("cauchit")), "ybin"),
  stats_binomial_cloglog = list(y | trials(tr) ~ x, quote(stats::binomial("cloglog")), "ybin"),
  stats_binomial_fn = list(y | trials(tr) ~ x, quote(stats::binomial), "ybin"),
  stats_quasipoisson = list(y ~ x, quote(stats::quasipoisson()), "ycnt"),
  # family as a function
  fn_poisson = list(y ~ x, quote(poisson), "ycnt"),
  fn_bernoulli = list(y ~ x, quote(bernoulli), "y01"),
  fn_negbinomial = list(y ~ x, quote(negbinomial), "ycnt"),
  fn_cumulative = list(y ~ x, quote(cumulative), "yord_f"),
  # brms family objects (built in the frmtmb arms too: brms is loaded
  # namespace-only there, so the object is brms's)
  brms_bernoulli_probit = list(y ~ x, quote(brms::bernoulli("probit")), "y01"),
  brms_student = list(y ~ x, quote(brms::student()), "yg"),
  brms_student_sigma_ident = list(y ~ x, quote(brms::student(link_sigma = "identity")), "ypos"),
  brms_negbinomial_sqrt = list(y ~ x, quote(brms::negbinomial("sqrt")), "ycnt"),
  brms_Beta = list(y ~ x, quote(brms::Beta()), "yprop"),
  brms_weibull_identity = list(y ~ x, quote(brms::weibull("identity")), "ypos"),
  brms_lognormal = list(y ~ x, quote(brms::lognormal()), "ypos"),
  brms_shifted_lognormal = list(y ~ x, quote(brms::shifted_lognormal()), "ypos"),
  brms_exgaussian = list(y ~ x, quote(brms::exgaussian()), "yg"),
  brms_skew_normal = list(y ~ x, quote(brms::skew_normal()), "yg"),
  brms_zip = list(y ~ x, quote(brms::zero_inflated_poisson()), "yzi"),
  brms_hurdle_poisson = list(y ~ x, quote(brms::hurdle_poisson()), "yzi"),
  brms_hurdle_gamma = list(y ~ x, quote(brms::hurdle_gamma()), "ycen"),
  brms_cumulative_probit = list(y ~ x, quote(brms::cumulative("probit")), "yord_f"),
  brms_sratio_cloglog = list(y ~ x, quote(brms::sratio("cloglog")), "yord_f"),
  brms_von_mises = list(y ~ x, quote(brms::von_mises()), "yangle"),
  brms_asym_laplace = list(y ~ x, quote(brms::asym_laplace()), "yg"),
  brms_invgauss = list(y ~ x, quote(brms::inverse.gaussian()), "ypos"),
  brms_geometric = list(y ~ x, quote(brms::geometric()), "ycnt"),
  brms_exponential = list(y ~ x, quote(brms::exponential()), "ypos"),
  brms_brmsfamily_gamma = list(y ~ x, quote(brms::brmsfamily("gamma", "identity")), "ypos"),
  brms_brmsfamily_zi_nb = list(y ~ x, quote(brms::brmsfamily("zero_inflated_negbinomial")), "yzi"),
  # distributional, random effects, addition terms, with non-default links
  dist_student_identity = list(quote(F(y ~ x, sigma ~ x)), quote(student(identity)), "yg"),
  dist_negbin_sqrt = list(quote(F(y ~ x, shape ~ 1)), quote(negbinomial(sqrt)), "ycnt"),
  dist_zip_identity = list(quote(F(y ~ x, zi ~ x)), quote(zero_inflated_poisson("identity")), "yzi"),
  re_bernoulli_cloglog = list(y ~ x + (1 | g), quote(bernoulli("cloglog")), "y01"),
  re_poisson_sqrt = list(y ~ x + (1 | g), quote(poisson("sqrt")), "ycnt"),
  re_gaussian_log = list(y ~ x + (1 | g), quote(gaussian("log")), "ypos"),
  cens_weibull_identity = list(y | cens(cen) ~ x, quote(weibull("identity")), "ycen"),
  cens_lognormal = list(y | cens(cen) ~ x, quote(lognormal()), "ycen"),
  trunc_gaussian_log = list(y | trunc(lb = 0) ~ x, quote(gaussian("log")), "ypos"),
  trials_binomial_probit = list(y | trials(tr) ~ x, quote(binomial("probit")), "ybin"),
  trials_betabin_cloglog = list(y | trials(tr) ~ x, quote(beta_binomial("cloglog")), "ybin"),
  trials_zib = list(y | trials(tr) ~ x, quote(zero_inflated_binomial()), "ybin"),
  # the bernoulli message
  msg_binomial_trials1 = list(y | trials(one) ~ x, quote(binomial()), "y01"),
  msg_binomial_notrials = list(y ~ x, quote(binomial()), "y01"),
  msg_betabin_trials1 = list(y | trials(one) ~ x, quote(beta_binomial()), "y01"),
  msg_cumulative_2cat = list(y ~ x, quote(cumulative()), "y2cat"),
  msg_cumulative_2cat_f = list(y ~ x, quote(cumulative()), "y2cat_f"),
  msg_categorical_2 = list(y ~ x, quote(categorical()), "y2cat"),
  msg_bernoulli = list(y ~ x, quote(bernoulli()), "y01"),
  msg_categorical_3 = list(y ~ x, quote(categorical()), "ycat"),
  msg_binomial_trials_var = list(y | trials(tr) ~ x, quote(binomial()), "ybin"),
  # ordinal response shapes
  ord_unordered_factor = list(y ~ x, quote(cumulative()), "ycat"),
  ord_integer = list(y ~ x, quote(cumulative()), "yord"),
  cat_unordered = list(y ~ x, quote(categorical()), "ycat"),
  # mixtures
  mix_gauss_gauss = list(y ~ 1, quote(mixture(gaussian, gaussian)), "ymix"),
  mix_gauss_student = list(y ~ 1, quote(mixture(gaussian, student)), "ymix"),
  mix_poisson_poisson = list(y ~ 1, quote(mixture(poisson, poisson)), "ymixc"),
  mix_poisson_negbin = list(y ~ 1, quote(mixture(poisson, negbinomial)), "ymixc"),
  mix_poisson_zip = list(y ~ 1, quote(mixture(poisson, zero_inflated_poisson)), "ymixc"),
  mix_poisson_hurdle = list(y ~ 1, quote(mixture(poisson, hurdle_poisson)), "ymixc"),
  mix_binom_betabin = list(y | trials(tr) ~ 1, quote(mixture(binomial, beta_binomial)), "ybin"),
  mix_poisson_binomial = list(y | trials(tr) ~ 1, quote(mixture(poisson, binomial)), "ybin"),
  mix_bern_bern = list(y ~ 1, quote(mixture(bernoulli, bernoulli)), "y01"),
  mix_lognormal_gamma = list(y ~ 1, quote(mixture(lognormal, Gamma)), "ypos"),
  mix_cum_cum = list(y ~ 1, quote(mixture(cumulative, cumulative)), "yord_f"),
  mix_zibeta_beta = list(y ~ 1, quote(mixture(Beta, zero_inflated_beta)), "yprop"),
  mix_exg_shiftln = list(y ~ 1, quote(mixture(exgaussian, shifted_lognormal)), "ypos"),
  mix_vonmises = list(y ~ 1, quote(mixture(von_mises, von_mises)), "yangle"),
  mix_gauss_hurdle_ln = list(y ~ 1, quote(mixture(gaussian, hurdle_lognormal)), "ypos"),
  # mvbf with per-response families
  mv_gauss_pois_sqrt = list(quote(mvb(F(y1 ~ x) + gaussian(), F(y2 ~ x) + poisson("sqrt"))), NULL, NULL),
  mv_student_negbin = list(quote(mvb(F(y1 ~ x + (1 | g)) + student(identity), F(y2 ~ x) + negbinomial(identity))), NULL, NULL),
  # added in the second pass
  brms_brmsfamily_invgauss = list(y ~ x, quote(brms::brmsfamily("inverse.gaussian")), "ypos"),
  brms_mixture_gauss = list(y ~ 1, quote(brms::mixture(gaussian, gaussian)), "ymix"),
  mix_gauss_ziasym = list(y ~ 1, quote(mixture(gaussian, zero_inflated_asym_laplace)), "ymix"),
  mix_hurdle_pois2 = list(y ~ 1, quote(mixture(hurdle_poisson, hurdle_poisson)), "ymixc"),
  mix_zip_zinb = list(y ~ 1, quote(mixture(zero_inflated_poisson, zero_inflated_negbinomial)), "ymixc"),
  mix_zib_binom = list(y | trials(tr) ~ 1, quote(mixture(zero_inflated_binomial, binomial)), "ybin"),
  mix_gauss_exgauss_str = list(y ~ 1, quote(mixture("gaussian", "exgaussian")), "ymix"),
  mix_weibull_exponential = list(y ~ 1, quote(mixture(weibull, exponential)), "ypos"),
  mix_student_skew = list(y ~ 1, quote(mixture(student, skew_normal)), "ymix"),
  mix_asym_asym = list(y ~ 1, quote(mixture(asym_laplace, asym_laplace)), "ymix")
)
if (!identical(ARM, "brms")) {
  # names frmtmb uses that brms spells differently
  cases$mix_lognormal_gamma[[2]] <- quote(mixture(lognormal, Gamma))
}
res <- list()
for (nm in names(cases)) {
  cs <- cases[[nm]]
  fo <- cs[[1]]
  if (is.call(fo) && !inherits(fo, "formula")) {
    fo <- tryCatch(eval(fo), error = function(e) e)
  }
  fam <- if (inherits(fo, "error")) fo else tryCatch(if (is.call(cs[[2]]) || is.name(cs[[2]])) eval(cs[[2]]) else cs[[2]],
                  error = function(e) e)
  dd <- d
  if (!is.null(cs[[3]])) dd$y <- d[[cs[[3]]]]
  n_msg <- 0L; msgs <- character()
  t0 <- Sys.time()
  out <- if (inherits(fam, "error")) fam else tryCatch(withCallingHandlers({
    if (identical(ARM, "brms")) {
      if (is.null(cs[[2]])) brms::make_standata(fo, data = dd) else
        brms::make_standata(fo, data = dd, family = fam)
    } else {
      if (is.null(cs[[2]])) frm(fo, data = dd) else frm(fo, data = dd, family = fam)
    }
  }, message = function(m) {
    if (grepl("bernoulli", conditionMessage(m))) {
      n_msg <<- n_msg + 1L; msgs <<- c(msgs, conditionMessage(m))
    }
    invokeRestart("muffleMessage")
  }, warning = function(w) invokeRestart("muffleWarning")),
  error = function(e) e)
  r <- list(status = if (inherits(out, "error")) "error" else "ok",
            message = if (inherits(out, "error")) conditionMessage(out) else "",
            n_bern_msg = n_msg)
  if (!inherits(out, "error") && !identical(ARM, "brms")) {
    r$logLik <- as.numeric(logLik(out))
    r$est <- unlist(out$estimates)
    r$se <- tryCatch(sqrt(diag(vcov(out))), error = function(e) NULL)
    r$fitted <- tryCatch(unlist(fitted(out)), error = function(e) NULL)
    r$family_link <- tryCatch({
      f <- family(out)
      if (inherits(f, "frmtmb_family")) f[["link"]] else lapply(f, `[[`, "link")
    }, error = function(e) paste("ERROR", conditionMessage(e)))
    # does update() or a refit message again?
    n_msg <- 0L
    withCallingHandlers(tryCatch(invisible(update(out)), error = function(e) NULL),
      message = function(m) {
        if (grepl("bernoulli", conditionMessage(m))) n_msg <<- n_msg + 1L
        invokeRestart("muffleMessage")
      }, warning = function(w) invokeRestart("muffleWarning"))
    r$n_bern_msg_update <- n_msg
  }
  res[[nm]] <- r
  cat(sprintf("%-28s %-5s msg %d %s\n", nm, r$status, r$n_bern_msg,
              substr(r$message, 1, 110)))
}
saveRDS(res, sprintf("dev/famlink-rev-designs-%s.rds", ARM))
