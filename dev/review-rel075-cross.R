# Review: cross-lane cases no lane tested (brmsnames x sarfcor x bhaz x
# smallfix). Each case prints frmtmb's answer, and brms 2.23.0's where
# brms answers without a compile (default_prior, validate_prior).
#   Rscript dev/review-rel075-cross.R > dev/review-rel075-log/cross.txt
.libPaths(c("C:/Users/adf44/source/r/rellib-r13",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
mk <- "C:/Users/adf44/Documents/.R/Makevars.win"
Sys.setenv(R_MAKEVARS_USER = mk)
suppressMessages({
  library(frmtmb)
  library(frmtmb.spatial)
  library(frmtmb.survival)
})
cat("frmtmb", format(packageVersion("frmtmb")), "spatial",
    format(packageVersion("frmtmb.spatial")), "survival",
    format(packageVersion("frmtmb.survival")), "brms",
    format(packageVersion("brms")), "\n")
show <- function(label, expr) {
  cat("\n=====", label, "\n")
  r <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    cat("WARNING:", conditionMessage(w), "\n")
    invokeRestart("muffleWarning")
  }), error = function(e) {
    cat("ERROR [", paste(class(e), collapse = ","), "]:",
        conditionMessage(e), "\n")
    NULL
  })
  if (!is.null(r)) print(r)
  invisible(r)
}
ptab <- function(p) {
  p <- as.data.frame(p)
  keep <- intersect(c("prior", "class", "coef", "group", "resp", "dpar",
                      "nlpar", "lb", "ub", "source"), names(p))
  p[, keep]
}

# ---- data: a 6 x 6 rook lattice, row-standardized ----
k <- 6
g <- expand.grid(r = seq_len(k), c = seq_len(k))
A <- outer(seq_len(k^2), seq_len(k^2), function(i, j)
  as.numeric(abs(g$r[i] - g$r[j]) + abs(g$c[i] - g$c[j]) == 1))
W <- A / rowSums(A)
n <- nrow(W)
set.seed(20261024)
d <- data.frame(x = runif(n, 0, 1.5), z = rnorm(n))
IA <- diag(n) - 0.4 * W
d$y <- as.vector(solve(IA, 1 + 0.5 * exp(d$x) - 0.3 * d$z +
                          rnorm(n, 0, 0.5)))
V <- 0.5^abs(outer(1:n, 1:n, "-"))

# ---- 1. prior table: exp(x) + sar(W), frmtmb vs brms ----
f1 <- y ~ exp(x) + z + sar(W)
show("1a frmtmb default_prior(y ~ exp(x) + z + sar(W))",
     ptab(default_prior(f1, data = d, data2 = list(W = W))))
show("1b brms default_prior, same",
     ptab(brms::default_prior(f1, data = d, data2 = list(W = W))))
pr1 <- c(set_prior("normal(0, 0.5)", class = "lagsar"),
         set_prior("normal(0, 1)", class = "b", coef = "expx"),
         set_prior("normal(0, 2)", class = "b", coef = "exp(x)"))
show("1c frmtmb validate_prior: lagsar + coef expx (brms) ",
     ptab(validate_prior(c(set_prior("normal(0, 0.5)", class = "lagsar"),
                           set_prior("normal(0, 1)", class = "b",
                                     coef = "expx")),
                         f1, data = d, data2 = list(W = W))))
show("1d brms validate_prior, same",
     ptab(brms::validate_prior(c(brms::set_prior("normal(0, 0.5)",
                                                class = "lagsar"),
                                 brms::set_prior("normal(0, 1)", class = "b",
                                                 coef = "expx")),
                               f1, data = d, data2 = list(W = W))))
show("1e frmtmb validate_prior: both spellings of one coef (one slot)",
     ptab(validate_prior(pr1, f1, data = d, data2 = list(W = W))))
show("1f frmtmb validate_prior: lagsar on an error SAR (wrong class)",
     ptab(validate_prior(set_prior("normal(0, 0.5)", class = "lagsar"),
                         y ~ exp(x) + sar(W, type = "error"), data = d,
                         data2 = list(W = W))))
show("1g brms validate_prior: lagsar on an error SAR",
     ptab(brms::validate_prior(brms::set_prior("normal(0, 0.5)",
                                               class = "lagsar"),
                               y ~ exp(x) + sar(W, type = "error"), data = d,
                               data2 = list(W = W))))
# ---- 2. fit: exp(x) + sar(W) with priors, names ----
fit2 <- show("2a frm(y ~ exp(x) + z + sar(W)) fixef",
             frm(f1, data = d, data2 = list(W = W)))
if (!is.null(fit2)) {
  show("2b fixef", fixef(fit2))
  show("2c variables", variables(fit2))
  show("2d frmtmb_autocor_pars", frmtmb_autocor_pars(fit2))
  show("2e prior_summary", ptab(prior_summary(fit2)))
  show("2f drop1 sees exp(x)", drop1(fit2))
}
fit2p <- show("2g frm with lagsar prior (MAP) + expx prior", frm(
  f1, data = d, data2 = list(W = W),
  prior = c(set_prior("normal(0.4, 0.05)", class = "lagsar"),
            set_prior("normal(0, 1)", class = "b", coef = "expx"))))
if (!is.null(fit2p)) {
  show("2h prior_summary of MAP fit", ptab(prior_summary(fit2p)))
  show("2i rho with prior vs without",
       c(with = frmtmb_autocor_pars(fit2p)$coef,
         rho_prior = unname(fit2p$estimates$thetaac),
         rho_none = unname(fit2$estimates$thetaac)))
}
# ---- 3. bf(autocor = ~ sar(W)) with exp(x) ----
show("3a frmtmb default_prior(bf(y ~ exp(x), autocor = ~ sar(W)))",
     ptab(default_prior(bf(y ~ exp(x), autocor = ~ sar(W)), data = d,
                        data2 = list(W = W))))
show("3b brms default_prior, same",
     ptab(brms::default_prior(brms::bf(y ~ exp(x), autocor = ~ sar(W)),
                              data = d, data2 = list(W = W))))
f3 <- show("3c frm(bf(y ~ exp(x), autocor = ~ sar(W)))",
           frm(bf(y ~ exp(x), autocor = ~ sar(W)), data = d,
               data2 = list(W = W)))
if (!is.null(f3) && !is.null(fit2)) {
  f2b <- frm(y ~ exp(x) + sar(W), data = d, data2 = list(W = W))
  show("3d same logLik as y ~ exp(x) + sar(W)",
       c(as.numeric(logLik(f3)), as.numeric(logLik(f2b))))
  show("3e names", list(fixef = rownames(fixef(f3)),
                        vars = variables(f3)))
}
# ---- 4. fcor(V) with exp(x) and a predicted sigma ----
show("4a frmtmb default_prior(bf(y ~ exp(x) + fcor(V), sigma ~ exp(z)))",
     ptab(default_prior(bf(y ~ exp(x) + fcor(V), sigma ~ exp(z)), data = d,
                        data2 = list(V = V))))
show("4b brms default_prior, same",
     ptab(brms::default_prior(brms::bf(y ~ exp(x) + fcor(V),
                                       sigma ~ exp(z)),
                              data = d, data2 = list(V = V))))
f4 <- show("4c frm fcor + exp in both dpars", frm(
  bf(y ~ exp(x) + fcor(V), sigma ~ exp(z)), data = d, data2 = list(V = V)))
if (!is.null(f4)) show("4d fixef", fixef(f4))
# ---- 5. a dirichlet() on lagsar (smallfix timing x sarfcor bounds) ----
show("5a set_prior dirichlet on lagsar",
     set_prior("dirichlet(1)", class = "lagsar"))
show("5b validate_prior dirichlet on lagsar",
     ptab(validate_prior(set_prior("dirichlet(1)", class = "lagsar"), f1,
                         data = d, data2 = list(W = W))))
show("5c frm dirichlet on lagsar",
     frm(f1, data = d, data2 = list(W = W),
         prior = set_prior("dirichlet(1)", class = "lagsar")))
show("5d frm lagsar prior with lb/ub given",
     frm(f1, data = d, data2 = list(W = W),
         prior = set_prior("uniform(-1, 1)", class = "lagsar", lb = -1,
                           ub = 1)))
# ---- 6. cox + bhaz: brms names and order (bhaz x brmsnames x smallfix) --
set.seed(20261024)
ds <- data.frame(x = runif(150, 0, 1.5), z = rnorm(150))
ds$time <- rexp(150, exp(-1 + 0.5 * exp(ds$x) + 0.3 * ds$z))
ds$cens <- as.numeric(runif(150) < 0.2)
fs <- time | cens(cens) + bhaz(df = 4) ~ exp(x) + z
show("6a frmtmb default_prior cox + bhaz(df = 4) + exp(x)",
     ptab(default_prior(fs, data = ds, family = cox())))
show("6b brms default_prior, same",
     ptab(brms::default_prior(fs, data = ds, family = brms::cox())))
show("6c validate_prior: dirichlet on sbhaz + coef expx",
     ptab(validate_prior(c(set_prior("dirichlet(2)", class = "sbhaz"),
                           set_prior("normal(0, 1)", class = "b",
                                     coef = "expx")),
                         fs, data = ds, family = cox())))
show("6d brms validate_prior, same",
     ptab(brms::validate_prior(c(brms::set_prior("dirichlet(2)",
                                                class = "sbhaz"),
                                 brms::set_prior("normal(0, 1)", class = "b",
                                                 coef = "expx")),
                               fs, data = ds, family = brms::cox())))
fc <- show("6e frm cox + bhaz(df = 4) + exp(x)",
           frm(fs, data = ds, family = cox()))
if (!is.null(fc)) {
  show("6f fixef", fixef(fc))
  show("6g variables", variables(fc))
}
fcp <- show("6h frm cox with dirichlet(2) sbhaz prior",
            frm(fs, data = ds, family = cox(),
                prior = set_prior("dirichlet(2)", class = "sbhaz")))
if (!is.null(fcp)) show("6i prior_summary", ptab(prior_summary(fcp)))
show("6j dirichlet on b (not a simplex) refused at frm",
     frm(fs, data = ds, family = cox(),
         prior = set_prior("dirichlet(2)", class = "b")))
show("6k cox(bhaz = ) and bhaz() equal",
     {
       a <- frm(time | cens(cens) ~ exp(x) + z, data = ds,
                family = cox(df = 4))
       c(identical(logLik(a), logLik(fc)),
         identical(fixef(a), fixef(fc)))
     })
# ---- 7. the refusal: a saved sar() + exp() fit read in a session where
#         frmtmb.spatial is unregistered is handled in its own script ----
cat("\nDONE\n")
