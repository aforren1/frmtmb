# Reviewer, lane wt-conditions (recheck 1): designs that frmtmb accepts,
# run on one library arm, to find a false alarm of the round's new
# refusals: check_frame_variables() in assemble_frame(),
# check_newdata_frame() before each newdata model.frame(), and
# prior_dist_params. Each case records "OK" or the error class and
# message; the base and lane records are then compared by case.
#   Rscript dev/conditions-rev-falsealarm.R base|lane   (seed 20260917)
arm <- commandArgs(TRUE)[1L]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/conditions-lib"
            else "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
options(frmtmb.notices = FALSE, warn = 1)
set.seed(20260917)
n <- 60
d <- data.frame(y = rnorm(n), x = rnorm(n), z = rnorm(n),
                g = factor(rep(1:10, each = 6)),
                g1 = factor(sample(1:8, n, TRUE)),
                g2 = factor(sample(1:8, n, TRUE)),
                f = factor(rep(c("a", "b", "c"), 20)),
                w = runif(n, 0.5, 1.5), expo = runif(n, 1, 3),
                w1 = 0.5, w2 = 0.5, t = rep(1:6, 10),
                loc = factor(rep(1:6, 10)))
d$k <- rbinom(n, 10, 0.4); d$ntr <- 10L
d$cnt <- rpois(n, 2 * d$expo)
d$cc <- sample(c(0, 1), n, TRUE, prob = c(0.8, 0.2))
d$ylo <- d$y - 0.5; d$yhi <- d$y + 0.5
d$ym <- ifelse(seq_len(n) %% 9 == 0, NA, d$y)
d$xm <- ifelse(seq_len(n) %% 8 == 0, NA, d$x)
d$sdx <- 0.2; d$sdy <- 0.3
d$ord <- factor(sample(1:4, n, TRUE), ordered = TRUE)
d$ordx <- sample(0:3, n, TRUE)
d$yb <- rbinom(n, 1, 0.5)
d$M <- cbind(d$x, d$z)
d$fc <- as.character(d$f)
d$dt <- as.Date("2020-01-01") + seq_len(n)
d$yp <- 2 * exp(-0.4 * abs(d$x)) + rnorm(n, 0, 0.1)
d$ych <- sample(c("p", "q", "r"), n, TRUE)
A <- diag(10); dimnames(A) <- list(levels(d$g), levels(d$g))
Wadj <- matrix(0, 6, 6); for (i in 1:5) Wadj[i, i + 1] <- Wadj[i + 1, i] <- 1
dimnames(Wadj) <- list(levels(d$loc), levels(d$loc))
xe <- rnorm(n)                      # a variable in the formula environment
kpow <- 2                           # a scalar constant in the environment

cases <- list()
C <- function(name, expr) cases[[name]] <<- substitute(expr)
F <- function(...) frm(..., dry_run = "frame")

# variables outside data
C("env_var_global", F(y ~ x + xe, data = d))
C("env_var_in_function", (function() { lv <- rnorm(60); F(y ~ lv, data = d) })())
C("env_var_calling_fn", (function(dd) { zz <- dd$x * 2; frm(y ~ zz, data = dd) })(d))
C("env_closure_fn", (function() { tr <- function(v) v^2; F(y ~ tr(x), data = d) })())
C("env_scalar_in_call", F(y ~ I(x * kpow), data = d))
C("env_scalar_poly", F(y ~ poly(x, kpow), data = d))
C("formula_built_elsewhere", (function() { fo <- (function() { off <- rnorm(60); y ~ off })(); F(fo, data = d) })())
C("ns_qualified", F(y ~ splines::ns(x, df = 3), data = d))
C("stats_poly_qualified", F(y ~ stats::poly(x, 2), data = d))
C("base_I_qualified", F(y ~ base::I(x^2), data = d))
C("poly", F(y ~ poly(x, 2), data = d))
C("ns_attached", { suppressWarnings(library(splines)); F(y ~ ns(x, 3), data = d) })
C("I_log", F(y ~ I(x^2) + log(expo), data = d))
C("pi_constant", F(y ~ I(x * pi), data = d))
C("dot_formula", F(y ~ ., data = d[, c("y", "x", "z")]))
C("matrix_column", F(y ~ M, data = d))
C("date_column", suppressMessages(F(y ~ dt, data = d)))
C("char_factor", F(y ~ fc, data = d))
# aterms and special terms
C("offset", F(cnt ~ x + offset(log(expo)), data = d, family = poisson()))
C("rate", F(cnt | rate(expo) ~ x, data = d, family = poisson()))
C("weights", F(y | weights(w) ~ x, data = d))
C("trials_col", F(k | trials(ntr) ~ x, data = d, family = binomial()))
C("trials_literal", F(k | trials(10) ~ x, data = d, family = binomial()))
C("trials_env_vector", (function() { nt <- rep(10L, 60); F(k | trials(nt) ~ x, data = d, family = binomial()) })())
C("cens_right", F(y | cens(cc) ~ x, data = d))
C("cens_interval", F(ylo | cens(cc, yhi) ~ x, data = d))
C("trunc", F(y | trunc(lb = -5) ~ x, data = d))
C("mi_response", F(ym | mi() ~ x, data = d))
C("mi_predictor", F(bf(y ~ mi(xm)) + bf(xm | mi() ~ 1) + set_rescor(FALSE), data = d))
C("se", F(y | se(sdy) ~ x, data = d))
C("me", F(y ~ me(x, sdx), data = d))
C("mo", F(y ~ mo(ordx), data = d))
C("cs_ordinal", F(ord ~ cs(x), data = d, family = acat()))
C("s_by", F(y ~ s(x, by = f), data = d))
C("t2", F(y ~ t2(x, z), data = d))
C("gp", F(y ~ gp(x), data = d))
C("gp_by", F(y ~ gp(x, by = f), data = d))
C("re_simple", F(y ~ x + (1 | g), data = d))
C("re_slope_factor", F(y ~ x + (f | g), data = d))
C("re_nested", F(y ~ x + (1 | g1:g2), data = d))
C("re_uncorr", F(y ~ x + (x || g), data = d))
C("re_id", F(bf(y ~ (1 | p | g), sigma ~ (1 | p | g)), data = d))
C("gr_cov_data2", F(y ~ x + (1 | gr(g, cov = A)), data = d, data2 = list(A = A)))
C("gr_cov_env", F(y ~ x + (1 | gr(g, cov = A)), data = d))
C("mm", F(y ~ x + (1 | mm(g1, g2)), data = d))
C("mm_weights", F(y ~ x + (1 | mm(g1, g2, weights = cbind(w1, w2))), data = d))
C("car_data2", F(y ~ car(Wadj, gr = loc), data = d, data2 = list(Wadj = Wadj)))
C("ar", F(y ~ x + ar(time = t, gr = g), data = d))
C("arma", F(y ~ x + arma(time = t, gr = g), data = d))
C("cosy", F(y ~ x + cosy(time = t, gr = g), data = d))
C("sigma_dpar", F(bf(y ~ x, sigma ~ z), data = d))
C("categorical_char", suppressMessages(F(ych ~ x, data = d, family = categorical())))
C("cbind_binomial", F(cbind(k, ntr - k) ~ x, data = d, family = binomial()))
C("mvbf", F(mvbf(bf(y ~ x), bf(z ~ x)), data = d))
C("mvbf_plus", F(bf(y ~ x) + bf(z ~ x) + set_rescor(TRUE), data = d))
# nonlinear
C("nl_basic", F(bf(yp ~ a * exp(-b * abs(x)), a ~ 1, b ~ 1, nl = TRUE), data = d))
C("nl_covariate_dpar", F(bf(yp ~ a * exp(-b * abs(x)), a ~ 1 + z, b ~ 1 + (1 | g), nl = TRUE), data = d))
C("nl_env_scalar", F(bf(yp ~ a * abs(x)^kpow, a ~ 1, nl = TRUE), data = d))
C("nl_env_function", (function() { hill <- function(v, e) v / (e + v); F(bf(yp ~ a * hill(abs(x), e), a ~ 1, e ~ 1, nl = TRUE), data = d) })())
C("nl_nlpar_named_like_env", (function() { a <- 5; F(bf(yp ~ a * abs(x), a ~ 1, nl = TRUE), data = d) })())
C("nlf_lf", F(bf(yp ~ a * abs(x), nl = TRUE) + lf(a ~ 1 + z), data = d))
# data containers
C("data_list", F(y ~ x + f, data = as.list(d[, c("y", "x", "f")])))
C("data_env", F(y ~ x + f, data = list2env(as.list(d[, c("y", "x", "f")]))))
C("data_tibble", if (requireNamespace("tibble", quietly = TRUE)) F(y ~ x + f, data = tibble::as_tibble(d[, c("y", "x", "f")])) else "no tibble")
C("data_table", if (requireNamespace("data.table", quietly = TRUE)) F(y ~ x + f, data = data.table::as.data.table(d[, c("y", "x", "f")])) else "no data.table")
C("data_list_posixlt", F(y ~ x + tm, data = list(y = d$y, x = d$x, tm = as.POSIXlt("2020-01-01") + 3600 * seq_len(n))))
C("data_df_posixlt_I", suppressMessages(F(y ~ x + tm, data = transform(d, tm = I(as.POSIXlt("2020-01-01") + 3600 * seq_len(n))))))
C("column_named_get", F(y ~ get, data = transform(d, get = x)))
C("column_named_c", F(y ~ c, data = transform(d, c = x)))
C("subset_expression", F(y ~ x, data = d[d$x > -1, ]))
C("na_in_predictor", suppressMessages(F(y ~ xm, data = d)))

# prior strings that parse on base
priors <- c("normal(0, 1)", "normal(0,1)", " normal( 0 , 10 ) ",
            "student_t(3, 0, 2.5)", "cauchy(0, 1)", "exponential(1)",
            "exponential(0.5)", "lkj(2)", "lkj(1)", "logistic(0, 1)",
            "gamma(2, 0.1)", "inv_gamma(1, 1)", "beta(1, 1)",
            "normal(0, 1e3)", "normal(-2, .5)", "student_t(3,0,2.5)",
            "constant(1)", "lkj_corr_cholesky(2)", "horseshoe(1)",
            "R2D2(0.5, 2)", "uniform(0, 1)", "normal(mu = 0, sigma = 1)",
            "", "flat", "exponential(2e-1)")
for (p in priors) cases[[paste0("prior_", p)]] <- bquote(set_prior(.(p)))
cases[["prior_bounds_sd"]] <- quote(frm(y ~ x + (1 | g), data = d, prior = set_prior("student_t(3, 0, 2)", class = "sd", lb = 0)))
cases[["prior_b_fit"]] <- quote(frm(y ~ x, data = d, prior = set_prior("normal(0, 5)", class = "b")))
cases[["prior_sigma_exp"]] <- quote(frm(y ~ x, data = d, prior = set_prior("exponential(1)", class = "sigma")))
cases[["prior_nl"]] <- quote(frm(bf(yp ~ a * exp(-b * abs(x)), a ~ 1, b ~ 1, nl = TRUE), data = d, prior = set_prior("normal(2, 1)", nlpar = "a") + set_prior("normal(0.5, 0.5)", nlpar = "b", lb = 0)))

# newdata paths, on fitted models
fits <- list(
  gm = quote(frm(y ~ x + f + (1 | g), data = d)),
  gs = quote(frm(y ~ x + f + (f | g), data = d)),
  bi = quote(frm(k | trials(ntr) ~ x + f, data = d, family = binomial())),
  po = quote(frm(cnt ~ x + offset(log(expo)), data = d, family = poisson())),
  wt = quote(frm(y | weights(w) ~ x, data = d)),
  sb = quote(frm(y ~ s(x, by = f), data = d)),
  mmf = quote(frm(y ~ x + (1 | mm(g1, g2)), data = d)),
  pl = quote(frm(y ~ poly(x, 2) + fc, data = d)),
  en = quote(frm(y ~ x + xe, data = d)),
  nl = quote(frm(bf(yp ~ a * exp(-b * abs(x)), a ~ 1 + f, b ~ 1, nl = TRUE), data = d, start = list(beta = c(a_Intercept = 2, b_Intercept = 0.4)))),
  mv = quote(frm(mvbf(bf(y ~ x + (1 | g)), bf(z ~ x + (1 | g))), data = d)),
  ce = quote(frm(y | cens(cc) ~ x + f, data = d)))
fitted_models <- list()
for (nm in names(fits)) {
  fitted_models[[nm]] <- tryCatch(suppressWarnings(eval(fits[[nm]])),
                                  error = function(e) e)
}
nd <- d[1:6, ]
nd_newg <- transform(nd, g = factor(c("99", "1", "2", "3", "4", "98")))
nd_nog <- nd[, setdiff(names(nd), "g")]
nd_noy <- nd[, setdiff(names(nd), c("y", "yp", "cnt", "k"))]
nd_char <- transform(nd, f = as.character(f))
nd_nalev <- transform(nd, f = factor(c("a", NA, "b", "c", "a", "b")))
nd_drop <- droplevels(d[d$f == "a", ][1:4, ])
nd_extra_unused <- transform(nd, f = factor(as.character(f), levels = c("a", "b", "c", "zz")))
nd_newf <- transform(nd, f = factor(c("a", "zz", "b", "c", "a", "b")))
nd_now <- nd[, setdiff(names(nd), "w")]
nd_notrials <- nd[, setdiff(names(nd), "ntr")]
nd_noexpo <- nd[, setdiff(names(nd), "expo")]
nd_nocc <- nd[, setdiff(names(nd), "cc")]
nd_mm_new <- transform(nd, g1 = factor(c("77", "1", "2", "3", "4", "5")))
P <- function(name, expr) cases[[paste0("nd_", name)]] <<- substitute(expr)
P("gm_base", predict(fitted_models$gm, newdata = nd))
P("gm_newg_allow", predict(fitted_models$gm, newdata = nd_newg, allow_new_levels = TRUE))
P("gm_newg_noallow", predict(fitted_models$gm, newdata = nd_newg))
P("gm_nog_reNA", predict(fitted_models$gm, newdata = nd_nog, re_formula = NA))
P("gm_nog_keep", predict(fitted_models$gm, newdata = nd_nog))
P("gm_noy", predict(fitted_models$gm, newdata = nd_noy))
P("gm_char", predict(fitted_models$gm, newdata = nd_char))
P("gm_nalev", predict(fitted_models$gm, newdata = nd_nalev))
P("gm_drop", predict(fitted_models$gm, newdata = nd_drop))
P("gm_extra_unused_level", predict(fitted_models$gm, newdata = nd_extra_unused))
P("gm_newf", predict(fitted_models$gm, newdata = nd_newf))
P("gm_newf_allow", predict(fitted_models$gm, newdata = nd_newf, allow_new_levels = TRUE))
P("gm_fitted_nd", fitted(fitted_models$gm, newdata = nd))
P("gm_simulate_nd", simulate(fitted_models$gm, newdata = nd, nsim = 2, seed = 1))
P("gm_ce", conditional_effects(fitted_models$gm))
P("gm_emmeans", as.data.frame(emmeans::emmeans(fitted_models$gm, ~ f)))
P("gm_mfx_pred", marginaleffects::avg_predictions(fitted_models$gm, by = "f"))
P("gm_mfx_slopes", marginaleffects::avg_slopes(fitted_models$gm, variables = "x"))
P("gm_mfx_comp", marginaleffects::avg_comparisons(fitted_models$gm, variables = "f"))
P("gs_newg_allow", predict(fitted_models$gs, newdata = nd_newg, allow_new_levels = TRUE))
P("gs_reNA_nog", predict(fitted_models$gs, newdata = nd_nog, re_formula = NA))
P("gs_newf", predict(fitted_models$gs, newdata = nd_newf))
P("gs_ce", conditional_effects(fitted_models$gs))
P("bi_base", predict(fitted_models$bi, newdata = nd))
P("bi_notrials_link", predict(fitted_models$bi, newdata = nd_notrials, type = "link"))
P("bi_notrials_resp", predict(fitted_models$bi, newdata = nd_notrials))
P("bi_noy", predict(fitted_models$bi, newdata = nd_noy))
P("bi_ce", conditional_effects(fitted_models$bi))
P("bi_emmeans", as.data.frame(emmeans::emmeans(fitted_models$bi, ~ f)))
P("bi_mfx", marginaleffects::avg_predictions(fitted_models$bi))
P("po_base", predict(fitted_models$po, newdata = nd))
P("po_noexpo", predict(fitted_models$po, newdata = nd_noexpo))
P("po_ce", conditional_effects(fitted_models$po))
P("wt_now", predict(fitted_models$wt, newdata = nd_now))
P("wt_ce", conditional_effects(fitted_models$wt))
P("sb_base", predict(fitted_models$sb, newdata = nd))
P("sb_newf", predict(fitted_models$sb, newdata = nd_newf))
P("sb_ce", conditional_effects(fitted_models$sb))
P("mm_base", predict(fitted_models$mmf, newdata = nd))
P("mm_new_allow", predict(fitted_models$mmf, newdata = nd_mm_new, allow_new_levels = TRUE))
P("mm_reNA_nog", predict(fitted_models$mmf, newdata = nd[, setdiff(names(nd), c("g1", "g2"))], re_formula = NA))
P("pl_base", predict(fitted_models$pl, newdata = nd))
P("pl_factor_for_char", predict(fitted_models$pl, newdata = transform(nd, fc = factor(fc))))
P("pl_ce", conditional_effects(fitted_models$pl))
P("en_base", predict(fitted_models$en, newdata = nd))
P("en_noxe_col", predict(fitted_models$en, newdata = transform(nd, xe = 1)))
P("nl_base", predict(fitted_models$nl, newdata = nd))
P("nl_ce", conditional_effects(fitted_models$nl))
P("mv_base", predict(fitted_models$mv, newdata = nd))
P("mv_newg_allow", predict(fitted_models$mv, newdata = nd_newg, allow_new_levels = TRUE))
P("ce_nocc", predict(fitted_models$ce, newdata = nd_nocc))
P("ce_ce", conditional_effects(fitted_models$ce))

res <- lapply(names(cases), function(nm) {
  out <- tryCatch({
    v <- suppressWarnings(suppressMessages(eval(cases[[nm]])))
    if (inherits(v, "error")) paste("FIT-ERROR", conditionMessage(v)) else "OK"
  }, error = function(e) paste0("ERROR [", class(e)[1L], "] ",
                               gsub("[[:space:]]+", " ", conditionMessage(e))))
  data.frame(case = nm, result = out, stringsAsFactors = FALSE)
})
out <- do.call(rbind, res)
fitstat <- vapply(fitted_models, function(v) if (inherits(v, "error"))
  paste("FIT-ERROR", conditionMessage(v)) else "fitted", "")
saveRDS(list(cases = out, fits = fitstat),
        sprintf("C:/Users/adf44/source/r/frmtmb-wt-conditions/dev/conditions-rev-log/falsealarm-%s.rds", arm))
cat("cases", nrow(out), " OK", sum(out$result == "OK"), "\n")
print(fitstat)
