# Review of lane priors2, items 3 and 4: the class theta density on
# frm() against the hand formula and against frm_sample()'s tape; the
# rows against brms 2.23.0; old spellings of the covariance class.
# Rscript priors2-review-theta.R
.libPaths(c("C:/Users/adf44/source/r/wt-priors2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.sample)
})
cat("lib:", find.package("frmtmb"), find.package("frmtmb.sample"), "\n")
set.seed(5)
n <- 150
d <- data.frame(x = rnorm(n))
k <- sample(1:3, n, TRUE, prob = c(0.5, 0.2, 0.3))
d$y <- c(-2, 1, 4)[k] + rnorm(n, sd = 0.7)
d$y2 <- rnorm(n)
fam3 <- mixture(gaussian(), gaussian(), gaussian())
a <- c(3, 1, 2)
pr <- set_prior("dirichlet(c(3, 1, 2))", class = "theta")
f0 <- frm(bf(y ~ 1), data = d, family = fam3, dry_run = "objective")
ri <- frmtmb:::resolve_prior_input(f0, pr)
e <- ri$entries
nlp <- frmtmb:::neg_log_prior_fn(e)
idx <- e[[1]]$idx
cat("entry idx", idx, "ref", e[[1]]$dist$ref, "\n")
pars <- f0$estimates
pts <- list(c(0, 0), c(0.4, -1.1), c(-2, 1.5), c(3, 3))
hand <- sapply(pts, function(z) {
  v <- c(z, 0)  # ref = K
  w <- exp(v) / sum(exp(v))
  lgamma(sum(a)) - sum(lgamma(a)) + sum(a * log(w))
})
frmv <- sapply(pts, function(z) {
  p <- pars
  p$betad[idx] <- z
  -nlp(p)
})
cat("frm prior log density at 4 points:", format(frmv, digits = 15), "\n")
cat("hand Dirichlet x softmax Jacobian: ", format(hand, digits = 15), "\n")
cat("max abs diff:", max(abs(frmv - hand)), "\n")
# the sampler's tape: the entry it adds where the user's row covers the
# simplex, and the default it drops
rp <- suppressMessages(frmtmb.sample:::sample_resolve_priors(f0, pr))
fs <- f0
fs$prior <- rp$effective
cat("mix_theta_nlp steps aside:",
    is.null(frmtmb.sample:::mix_theta_nlp(fs, rp$ri$entries)), "\n")
ents <- frmtmb.sample:::mo_sample_entries(rp$ri$entries, fs)
de <- Filter(function(z) identical(z$dist$kind, "dirichlet"), ents)
smpv <- sapply(pts, function(z) {
  p <- pars
  p$betad[idx] <- z
  -frmtmb:::neg_log_prior_fn(de)(p)
})
cat("sampler tape log density:          ", format(smpv, digits = 15), "\n")
cat("sampler - frm (constant?):", format(smpv - frmv, digits = 15), "\n")
# the MAP under a strong asymmetric prior moves theta1 up
f_ml <- frm(bf(y ~ 1), data = d, family = fam3)
f_map <- frm(bf(y ~ 1), data = d, family = fam3,
             prior = set_prior("dirichlet(c(200, 1, 1))", class = "theta"))
w_of <- function(f) {
  z <- f$estimates$betad[idx]
  v <- c(z, 0)
  exp(v) / sum(exp(v))
}
cat("weights ML:", format(w_of(f_ml), digits = 6), " code",
    f_ml$opt$convergence, "\n")
cat("weights dirichlet(200,1,1):", format(w_of(f_map), digits = 6), " code",
    f_map$opt$convergence, "\n")
# rows: 2-, 3-component, multivariate, against brms
show_rows <- function(lab, tab) {
  r <- tab[tab$class %in% c("theta", "thetacov") |
             startsWith(tab$class, "theta"), ]
  cat(lab, ": ", if (nrow(r)) paste(r$prior, r$class, r$resp,
                                     sep = "|", collapse = "; ") else "none",
      "\n", sep = "")
}
show_rows("frm fit 2", default_prior(bf(y ~ x), data = d,
                                     family = mixture(gaussian(),
                                                      gaussian())))
show_rows("frm sample 3", default_prior(bf(y ~ x), data = d, family = fam3,
                                        route = "sample"))
show_rows("frm fit mv", default_prior(
  bf(y ~ x, family = mixture(gaussian(), gaussian())) +
    bf(y2 ~ x, family = mixture(gaussian(), gaussian())), data = d))
show_rows("frm fit mv one mixture", default_prior(
  bf(y ~ x, family = mixture(gaussian(), gaussian())) +
    bf(y2 ~ x, family = gaussian()), data = d))
show_rows("frm fit theta1 ~ x", default_prior(
  bf(y ~ x, theta1 ~ x), data = d, family = mixture(gaussian(), gaussian())))
suppressPackageStartupMessages(library(brms))
bshow <- function(lab, tab) {
  r <- tab[startsWith(tab$class, "theta"), ]
  cat(lab, ": ", if (nrow(r)) paste(r$prior, r$class, r$resp, sep = "|",
                                    collapse = "; ") else "none", "\n",
      sep = "")
}
bshow("brms 2", default_prior(bf(y ~ x), data = d,
                              family = brms::mixture(gaussian, gaussian)))
bshow("brms 3", default_prior(bf(y ~ x), data = d,
                              family = brms::mixture(gaussian, gaussian,
                                                     gaussian)))
bshow("brms mv one mixture", default_prior(
  bf(y ~ x, family = brms::mixture(gaussian, gaussian)) +
    bf(y2 ~ x, family = gaussian()), data = d))
bshow("brms theta1 ~ x", default_prior(
  bf(y ~ x, theta1 ~ x), data = d, family = brms::mixture(gaussian, gaussian)))
detach("package:brms")
# old spellings of the covariance class, on models without a mixture
try_msg <- function(lab, expr) {
  r <- tryCatch(expr, error = function(e) paste("ERROR:",
                                                conditionMessage(e)))
  if (!is.character(r)) r <- "accepted"
  cat(lab, ":", substr(r, 1, 400), "\n")
}
d$g <- factor(rep(1:10, 15))
d$tf <- factor(rep(1:3, 50))
ft <- bf(y2 ~ x + toep(0 + tf | g))
try_msg("old theta coef theta_4 normal on toep", validate_prior(
  set_prior("normal(0, 1)", class = "theta", coef = "theta_4"), ft,
  data = d))
try_msg("old theta coef theta_4 bound-only on toep", validate_prior(
  set_prior("", class = "theta", coef = "theta_4", lb = -1, ub = 1), ft,
  data = d))
try_msg("old theta coef theta_4 empty on toep", validate_prior(
  set_prior("", class = "theta", coef = "theta_4"), ft, data = d))
try_msg("old theta class-wide empty on toep", validate_prior(
  set_prior("", class = "theta"), ft, data = d))
try_msg("frm old theta empty on toep", frm(ft, data = d,
  prior = set_prior("", class = "theta"), dry_run = "objective"))
try_msg("theta dirichlet on non-mixture", validate_prior(
  set_prior("dirichlet(2)", class = "theta"), ft, data = d))
try_msg("theta resp on mv, resp not mixture", validate_prior(
  set_prior("dirichlet(2)", class = "theta", resp = "y2"),
  bf(y ~ x, family = mixture(gaussian(), gaussian())) +
    bf(y2 ~ x, family = gaussian()), data = d))
try_msg("theta empty on mv mixture no resp", validate_prior(
  set_prior("", class = "theta"),
  bf(y ~ x, family = mixture(gaussian(), gaussian())) +
    bf(y2 ~ x, family = gaussian()), data = d))
try_msg("theta on univariate with resp", validate_prior(
  set_prior("dirichlet(2)", class = "theta", resp = "y"),
  bf(y ~ x), data = d, family = mixture(gaussian(), gaussian())))
# a prior table passed back: fit route and sample route
f2 <- frm(bf(y ~ x), data = d, family = mixture(gaussian(), gaussian()))
fb <- frm(bf(y ~ x), data = d, family = mixture(gaussian(), gaussian()),
          prior = default_prior(f2))
cat("fit-route table round trip bitwise:",
    identical(fb$opt$par, f2$opt$par), "\n")
fs2 <- suppressWarnings(frm(bf(y ~ x), data = d,
                            family = mixture(gaussian(), gaussian()),
                            prior = default_prior(f2, route = "sample")))
cat("sample-route table on frm(): code", fs2$opt$convergence,
    "theta row", default_prior(f2, route = "sample")$prior[
      default_prior(f2, route = "sample")$class == "theta"], "\n")
