# Reviewer: are the "refusal correct" / "compat row wrong (both refuse)"
# verdicts refusals of the FEATURE, or of the spec's spelling? Each case
# tries the brms spelling a user would write, builds it with
# make_stancode(), and compiles it through the copied cache.
.libPaths(c("C:/Users/adf44/source/r/rellib-r15",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(FRMTMB_STAN_CACHE = normalizePath("dev/review-equiv-cache",
                                             winslash = "/"),
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({library(frmtmb); library(testthat)})
E <- new.env(parent = asNamespace("frmtmb"))
for (f in c("dev/equiv-lib.R", "dev/equiv-specs.R", "dev/equiv-flow.R")) {
  sys.source(f, envir = E)
}
say <- function(lab, r) cat(sprintf("%-58s %s\n", lab,
  if (inherits(r, "error")) paste("ERROR:", substr(conditionMessage(r), 1,
                                                    150)) else "ok"))
tr <- function(expr) tryCatch(suppressWarnings(expr), error = function(e) e)
comp <- function(code) tr({E$brms_stan_model(code); TRUE})

# p158 fcor() x nl: brms reaches fcor() in an nl model through
# autocor = ~ fcor(M) (the compat note says so)
s <- E$eq_specs[["fcor() x nl"]]()
d <- s$data; d2 <- s$data2
cat("p158 data2 names:", names(d2), "\n")
b <- tr(brms::make_stancode(brms::bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1,
                                     nl = TRUE, autocor = ~ fcor(M)),
                            data = d, data2 = d2))
say("p158 brms bf(nl, autocor = ~ fcor(M)) make_stancode", b)
if (!inherits(b, "error")) say("p158 brms compile", comp(b))
b2 <- tr(brms::make_stancode(brms::bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1,
                                      nl = TRUE) + brms::fcor(M),
                             data = d, data2 = d2))
say("p158 brms bf(nl) + fcor(M) make_stancode", b2)
f1 <- tr(frm(bf(y ~ b1 * exp(b2 * x), b1 + b2 ~ 1, nl = TRUE,
                autocor = ~ fcor(M)), data = d, data2 = d2))
say("p158 frm bf(nl, autocor = ~ fcor(M))", f1)
f2 <- tr(frm(bf(y ~ b1 * exp(b2 * x) + fcor(M), b1 + b2 ~ 1, nl = TRUE),
             data = d, data2 = d2))
say("p158 frm fcor in body (lane spelling)", f2)

# p167 poisson x arma(): brms refuses MA under cov = FALSE only
s <- E$eq_specs[["poisson x arma()"]]()
d <- s$data
b <- tr(brms::make_stancode(brms::bf(y ~ x + arma(time = t, gr = g,
                                                    cov = TRUE)),
                            family = brms::poisson(), data = d))
say("p167 brms poisson arma(cov = TRUE) make_stancode", b)
if (!inherits(b, "error")) say("p167 brms compile", comp(b))
f1 <- tr(frm(bf(y ~ x + arma(time = t, gr = g, cov = TRUE)),
             family = poisson(), data = d))
say("p167 frm poisson arma(cov = TRUE)", f1)
b <- tr(brms::make_stancode(brms::bf(y ~ x + ar(time = t, gr = g)),
                            family = brms::poisson(), data = d))
say("p167 brms poisson ar() (cov = FALSE) make_stancode", b)
f1 <- tr(frm(bf(y ~ x + ar(time = t, gr = g)), family = poisson(), data = d))
say("p167 frm poisson ar()", f1)

# p014 cs_pred() x nl: cs() in a linear nlpar formula, not the body
s <- E$eq_specs[["cs_pred() x nl"]]()
d <- s$data
b <- tr(brms::make_stancode(brms::bf(y ~ b1 * x + eta, b1 ~ 1,
                                     eta ~ 0 + cs(z), nl = TRUE),
                            family = brms::sratio(), data = d))
say("p014 brms nlpar eta ~ 0 + cs(z) make_stancode", b)
if (!inherits(b, "error")) say("p014 brms compile", comp(b))
f1 <- tr(frm(bf(y ~ b1 * x + eta, b1 ~ 1, eta ~ 0 + cs(z), nl = TRUE),
             family = sratio(), data = d))
say("p014 frm nlpar eta ~ 0 + cs(z)", f1)
if (!inherits(f1, "error")) cat("   frm fit holds", length(unlist(
  f1$estimates$beta)), "coefficients; logLik", logLik(f1), "\n")

# p093 hurdle_cumulative residuals: what does brms say?
s <- E$eq_specs[["hurdle_cumulative x residuals"]]()
fit <- E$eq_frm_fit(s)
bf1 <- E$eq_fixed(s, brms::bf(y ~ x + z), brms::hurdle_cumulative(),
                  s$data, fit)
r <- tr(residuals(bf1, method = "posterior_epred", type = "ordinary",
                  summary = FALSE))
say("p093 brms residuals(hurdle_cumulative)", r)

# p154 cens() x cosy(): brms message names the reason?
s <- E$eq_specs[["cens() x cosy()"]]()
b <- tr(brms::make_stancode(brms::bf(y | cens(c) ~ x + cosy(gr = g)),
                            data = s$data))
say("p154 brms cens + cosy", b)
