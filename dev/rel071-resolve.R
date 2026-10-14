# Consolidation edits of the 0.71.0 release tree, applied after the
# textual merge (dev/rel071-merge.sh) and before the NEWS merge
# (dev/rel071-news-merge.R). Each edit checks that its old text is
# there, so a rerun on a fresh merge either applies all of them or
# stops on the first that no longer fits. A file no lane changed is
# restored from its HEAD blob (253b3ef7, 0.70.1) before it is edited,
# so every edit here is applied once however often this runs.
#
# REL071_LANES (set by the merge script) names the lanes merged; the
# edits that belong to lane illse apply only when it is among them.
#
#   Rscript dev/rel071-resolve.R
rel <- "C:/Users/adf44/source/r/frmtmb-wt-release"
root <- "C:/Users/adf44/source/r"
setwd(rel)
lanes <- strsplit(Sys.getenv("REL071_LANES",
                             "dblbar famsimplex famuni postfix optstart remeasure"),
                  " ")[[1]]
with_illse <- "illse" %in% lanes
cat("resolve for lanes:", lanes, "\n")
rd <- function(f) sub("\r$", "", readLines(f, warn = FALSE))
wr <- function(x, f) {
  force(x)  # the argument may read the file that file() truncates
  con <- file(f, "wb")
  writeLines(x, con, sep = "\n")
  close(con)
}
from_head <- function(f) {
  x <- system2("git", c("show", paste0("HEAD:", f)), stdout = TRUE)
  wr(sub("\r$", "", x), f)
}
edit_file <- function(f, old, new, n = 1L) {
  s <- paste(rd(f), collapse = "\n")
  hits <- gregexpr(old, s, fixed = TRUE)[[1]]
  k <- sum(hits > 0)
  if (k != n) stop(f, ": expected ", n, " match, found ", k, ": ",
                   substr(old, 1, 70))
  s <- gsub(old, new, s, fixed = TRUE)
  wr(strsplit(s, "\n", fixed = TRUE)[[1]], f)
  cat("edited", f, "\n")
}

# ---- the four union conflicts of famsimplex and famuni ---------------
# Both lanes append families after the same last entry; the union keeps
# famsimplex's (merged first), then famuni's.
edit_file("R/compat.R", paste0(
'<<<<<<< release
            "dirichlet", "dirichlet2", "logistic_normal")
=======
            "frechet", "gen_extreme_value", "discrete_weibull")
>>>>>>> famuni'), paste0(
'            "dirichlet", "dirichlet2", "logistic_normal",
            "frechet", "gen_extreme_value", "discrete_weibull")'))
edit_file("R/families.R", paste0(
"<<<<<<< release
  # the simplex families (R/families-simplex.R); dirichlet2 has no
  # constructor, as in brms, and is reached through brmsfamily()
  dirichlet                 = fam_dirichlet_deferred,
  dirichlet2                = fam_dirichlet2_deferred,
  logistic_normal           = fam_logistic_normal_deferred
=======
  # R/families-extreme.R
  frechet                   = fam_frechet,
  gen_extreme_value         = fam_gen_extreme_value,
  discrete_weibull          = fam_discrete_weibull
>>>>>>> famuni"), paste0(
"  # the simplex families (R/families-simplex.R); dirichlet2 has no
  # constructor, as in brms, and is reached through brmsfamily()
  dirichlet                 = fam_dirichlet_deferred,
  dirichlet2                = fam_dirichlet2_deferred,
  logistic_normal           = fam_logistic_normal_deferred,
  # R/families-extreme.R
  frechet                   = fam_frechet,
  gen_extreme_value         = fam_gen_extreme_value,
  discrete_weibull          = fam_discrete_weibull"))
edit_file("R/links.R", paste0(
"<<<<<<< release
#' | `categorical`, `multinomial`, `dirichlet` | `logit` |
#' | `logistic_normal` | `identity` |
#' | `dirichlet2`, through `brmsfamily()` | `log`, `softplus`, `squareplus`, `identity`, `logm1` |
=======
#' | `categorical`, `multinomial` | `logit` |
#' | `frechet` | `log`, `identity`, `inverse`, `softplus`, `squareplus` |
#' | `gen_extreme_value` | `identity`, `log`, `inverse`, `softplus`, `squareplus` |
#' | `discrete_weibull` | `logit`, `probit`, `probit_approx`, `cloglog`, `cauchit`, `softit`, `identity` |
>>>>>>> famuni"), paste0(
"#' | `categorical`, `multinomial`, `dirichlet` | `logit` |
#' | `logistic_normal` | `identity` |
#' | `dirichlet2`, through `brmsfamily()` | `log`, `softplus`, `squareplus`, `identity`, `logm1` |
#' | `frechet` | `log`, `identity`, `inverse`, `softplus`, `squareplus` |
#' | `gen_extreme_value` | `identity`, `log`, `inverse`, `softplus`, `squareplus` |
#' | `discrete_weibull` | `logit`, `probit`, `probit_approx`, `cloglog`, `cauchit`, `softit`, `identity` |"))
edit_file("dev/famlink-gen-links.R", paste0(
'<<<<<<< release
  dirichlet = "dirichlet", dirichlet2 = "dirichlet2",
  logistic_normal = "logistic_normal"
=======
  frechet = "frechet", gen_extreme_value = "gen_extreme_value",
  discrete_weibull = "discrete_weibull"
>>>>>>> famuni'), paste0(
'  dirichlet = "dirichlet", dirichlet2 = "dirichlet2",
  logistic_normal = "logistic_normal",
  frechet = "frechet", gen_extreme_value = "gen_extreme_value",
  discrete_weibull = "discrete_weibull"'))

# ---- R/se-check.R: ci701 (in HEAD) against illse -------------------------
# fit_outer_hessian(): ci701 added `pull` (the correlation past the
# clamp moved back to it before the Hessian is taken); illse keeps the
# last result on the objective, keyed on the point and the units. Both:
# the point is the pulled one, so the key is too, and a pulled Hessian
# is never served for the unpulled point (fit_outer_hessian_at(), which
# asks for the fit's own point, gets its own key). cov_from_hessian():
# ci701's `edge <- se_edge_sd(fit, H) | ec` stays ahead of illse's
# route, which runs only without an edge (the coordinator; illse
# findings section 10, item 8: with a correlation at the clamp, tier 3
# removes it and hands the rest to the route through `free`).
if (with_illse) {
  edit_file("R/se-check.R", paste0(
"<<<<<<< release
#' `pull` (logical, one per outer parameter) moves those parameters back
#' to the clamp before the Hessian is taken (se_edge_cor()).
=======
#' The last result is kept on the objective (`obj$env`), keyed on the
#' point and the units: conv_headroom() and grad_headroom() read it at
#' the fit's own point, as the SE check does, and on a Laplace objective
#' with 100 parameters and n = 1e4 each one cost 10 s (review of lane
#' illse, item 7).
>>>>>>> illse
#'
#' @noRd
fit_outer_hessian <- function(fit, pull = NULL) {
  obj <- fit$obj
  u <- fit$par_units
  scaled <- !(is.null(u) || all(u == 1))
  at <- if (scaled) fit$opt$par else sdr_outer_point(obj)
  env <- obj$env"), paste0(
"#' `pull` (logical, one per outer parameter) moves those parameters back
#' to the clamp before the Hessian is taken (se_edge_cor()).
#'
#' The last result is kept on the objective (`obj$env`), keyed on the
#' point (after `pull`) and the units: conv_headroom() and
#' grad_headroom() read it at the fit's own point, as the SE check does,
#' and on a Laplace objective with 100 parameters and n = 1e4 each one
#' cost 10 s (review of lane illse, item 7).
#'
#' @noRd
fit_outer_hessian <- function(fit, pull = NULL) {
  obj <- fit$obj
  u <- fit$par_units
  scaled <- !(is.null(u) || all(u == 1))
  at <- se_pull_back(if (scaled) fit$opt$par else sdr_outer_point(obj),
                     pull)
  env <- obj$env"))
  edit_file("R/se-check.R", paste0(
"<<<<<<< release
  u <- fit$par_units
  if (is.null(u) || all(u == 1)) {
    at <- se_pull_back(sdr_outer_point(obj), pull)
=======
  out <- if (!scaled) {
>>>>>>> illse
    h <- fd_hessian(at, obj$gr)
    list(H = h$H, at = at, Hq = NULL, E = h$E)
  } else {
    q0 <- fit$opt$par / u
    h <- fd_hessian(q0, function(q) obj$gr(q * u) * u)
    list(H = h$H / outer(u, u), at = fit$opt$par, Hq = h$H,
         E = h$E / outer(u, u))
  }"), paste0(
"  out <- if (!scaled) {
    h <- fd_hessian(at, obj$gr)
    list(H = h$H, at = at, Hq = NULL, E = h$E)
  } else {
    h <- fd_hessian(at / u, function(q) obj$gr(q * u) * u)
    list(H = h$H / outer(u, u), at = at, Hq = h$H,
         E = h$E / outer(u, u))
  }"))
  edit_file("R/se-check.R", paste0(
"<<<<<<< release
  at <- se_pull_back(fit$opt$par, pull)
  h <- fd_hessian(at / u, function(q) obj$gr(q * u) * u)
  list(H = h$H / outer(u, u), at = at, Hq = h$H,
       E = h$E / outer(u, u))
=======
  if (kept_only) {"), "  if (kept_only) {")
  edit_file("R/se-check.R", paste0(
"  h <- tryCatch(fit_outer_hessian(fit), error = function(e) NULL)
  H <- h$H
  if (is.matrix(H) && all(dim(H) == length(p))) H
>>>>>>> illse
}"), paste0(
"  h <- tryCatch(fit_outer_hessian(fit), error = function(e) NULL)
  H <- h$H
  if (is.matrix(H) && all(dim(H) == length(p))) H
}"))
  edit_file("R/se-check.R", paste0(
"<<<<<<< release
  edge <- se_edge_sd(fit, H) | ec
=======
  edge <- se_edge_sd(fit, H)
  # before the tiers"), paste0(
"  edge <- se_edge_sd(fit, H) | ec
  # before the tiers"))
  edit_file("R/se-check.R", paste0(
"    if (!is.null(rp)) return(rp)
  }
>>>>>>> illse
"), paste0(
"    if (!is.null(rp)) return(rp)
  }
"))
}

# ---- tests/testthat/test-se-conv.R: illse's case against optstart ---------
# Found at the consolidation. "an sd held near zero short of its maximum
# is not at the maximum" started only theta (-15), so on the merged tree
# lane optstart's least-squares start places beta and sigma, and nlminb
# converges at code 0 in 1 or 2 iterations to the maximum given sd = 0,
# 85 to 139 below the maximum, and 5 of its expectations fail
# (dev/rel071-seconv.R, dev/rel071-log/seconv.txt). The consolidation's
# first replacement started every parameter at 0; the release review
# showed that it no longer pins illse's R2-1 fix, because from beta =
# (0, 0) the intercept is short too (gradient 26.2), so the stop is
# called short with the fix removed (dev/relrev071-seconv.R). The
# review's replacement: start at lm()'s maximum given sd = 0 and take
# no step, so only the sd row is short. The large-gradient warning no
# longer applies there (every other gradient is at its maximum), so
# that assertion becomes "no warning calls the stop the maximum"
# (dev/rel071-seconv3.R, dev/rel071-log/seconv3-*.txt: fails on
# 0.70.0, fails with se_sd_gain_up() stubbed to 0, passes on the
# merged build).
if (with_illse) {
  ix <- rd("tests/testthat/test-se-conv.R")
  i0 <- grep("^test_that[(]\"an sd held near zero short of its maximum",
             ix)
  stopifnot(length(i0) == 1L)
  i1 <- i0 - 1L + which(ix[i0:length(ix)] == "})")[1L]
  stopifnot(any(grepl("start = list(theta = -15),", ix[i0:i1],
                      fixed = TRUE)))
  body <- c(
"test_that(\"an sd held near zero short of its maximum is not at the maximum\", {",
"  # The sd's row is empty and its gradient is 5.6e-10, so a Newton step",
"  # cannot see what a larger sd gains. The headroom asks the likelihood",
"  # (se_sd_gain_up()); without that, this stop 92 below the maximum was",
"  # called the maximum (review of lane illse, round 2, R2-1). The fit",
"  # starts at lm()'s maximum given sd = 0 and takes no step, so the sd",
"  # is the only direction that is short whatever start frm() would pick;",
"  # from a start where another coefficient is short too, the stop is",
"  # called short without R2-1 (the 0.71.0 release review)",
"  set.seed(1)",
"  n <- 300",
"  d <- data.frame(x = stats::rnorm(n),",
"                  g = factor(rep(1:20, length.out = n)))",
"  d$y <- 1 + 0.5 * d$x + stats::rnorm(20, 0, 0.5)[d$g] +",
"    stats::rnorm(n, 0, 0.4)",
"  ls0 <- stats::lm(y ~ x, data = d)",
"  start <- list(beta = unname(stats::coef(ls0)),",
"                betad = log(sqrt(mean(stats::residuals(ls0)^2))),",
"                theta = -15)",
"  r <- conv_conds(frm(y ~ x + (1 | g), data = d, start = start,",
"                      control = frmtmb_control(optCtrl = list(",
"                        iter.max = 0, eval.max = 1000))))",
"  expect_false(r$value$opt$convergence == 0L)",
"  w <- grep(conv_phrase, r$warnings, fixed = TRUE, value = TRUE)",
"  expect_length(w, 1L)",
"  expect_match(w, \"short of the maximum\", fixed = TRUE)",
"  expect_false(any(grepl(\"so the estimates are at the maximum\",",
"                         r$warnings, fixed = TRUE)))",
"})")
  wr(c(ix[seq_len(i0 - 1L)], body, ix[-seq_len(i1)]),
     "tests/testthat/test-se-conv.R")
  cat("edited tests/testthat/test-se-conv.R\n")
}

# ---- R/links-brms.R, generated once from both lanes' src entries -----
rs <- file.path(R.home("bin"), "Rscript.exe")
out <- system2(rs, "dev/famlink-gen-links.R", stdout = TRUE, stderr = TRUE)
cat(out, sep = "\n")
stopifnot(identical(attr(out, "status"), NULL))
wr(rd("R/links-brms.R"), "R/links-brms.R")
lb <- rd("R/links-brms.R")
for (nm in c("dirichlet", "dirichlet2", "logistic_normal", "frechet",
             "gen_extreme_value", "discrete_weibull")) {
  stopifnot(sum(grepl(paste0("^  ", nm, " = "), lb)) == 2L)
}
cat("links-brms: 6 new families, each in both tables\n")

# ---- tests/testthat/test-simulate-density.R: the count header --------
# famsimplex wrote 307 for its own rows; famuni's three sim_specs rows
# add two cells each of a chi-square, a mean and a variance. Counted on
# the merged build by dev/rel071-simcount.R (dev/rel071-log/simcount.txt:
# 109 chi-square, 117 mean, 117 variance, 18 P(Y = k) and 2 correlation;
# the script traces the bare expect_lt() of the correlation checks too,
# since the release review): 363. The header was already short at
# 0.70.0 (289 counted
# on HEAD's file against its 251, simcount-head-r8.txt; the P(Y = 1)
# atoms and ten cells were missing), so it now states the measured
# count. The increments agree: famsimplex 54 + 2, famuni 18. The
# arithmetic is 1 - (1 - a)^n.
edit_file("tests/testthat/test-simulate-density.R", paste0(
"# The trade depends on the COUNT, so it is written down here: 307 of
# this file's assertions are statistical (93 chi-square, 101 mean, 101
# variance, 10 P(Y = 0), 2 correlation); the simplex families added 56
# of them. At this alpha that is a 0.0003 chance that any reseed raises
# a permanent false alarm somewhere in the file; at 1e-3 it would be
# 0.26, and at a conventional 0.05 it would be a near-certainty."),
paste0(
"# The trade depends on the COUNT, so it is written down here: 363 of
# this file's assertions are statistical (109 chi-square, 117 mean, 117
# variance, 18 P(Y = 0) or P(Y = 1), 2 correlation; counted by
# dev/rel071-simcount.R). The simplex families added 56 of them, and
# frechet, gen_extreme_value and discrete_weibull 18. At this alpha
# that is a 0.0004 chance that any reseed raises a permanent false
# alarm somewhere in the file; at 1e-3 it would be 0.30, and at a
# conventional 0.05 it would be a near-certainty."))

# ---- com_poisson: frmtmb's mu is the mean (the user's decision) -------
edit_file("R/families.R", paste0(
"#' poisson family there. `threshold` is an argument of the ordinal
#' families alone, and `refcat` of `\"dirichlet\"` and
#' `\"logistic_normal\"`.
"), paste0(
"#' poisson family there. `threshold` is an argument of the ordinal
#' families alone, and `refcat` of `\"dirichlet\"` and
#' `\"logistic_normal\"`.
#'
#' `brmsfamily(\"com_poisson\")` gives [compois()], whose `mu` is the
#' MEAN of the count. brms 2.23.0 does not use the mean. Its Stan
#' program fits `mu` as the rate `lambda` of the Conway-Maxwell-Poisson
#' density, which is proportional to `lambda^y / (y!)^nu`, and its
#' `log_lik()` and `posterior_predict()` read `mu` as
#' `lambda^(1 / nu)`. So the coefficients of a brms `com_poisson()` fit
#' and of an frmtmb fit are on different scales. They agree only at
#' `nu = 1` (brms's `shape`), where all three are the Poisson mean.
"))
edit_file("R/families.R", paste0(
"#' @section Links:
#' Every constructor takes `link` for the mean and a `link_<dpar>` for
"), paste0(
"#' @section Conway-Maxwell-Poisson counts:
#' In `compois()`, `mu` is the mean of the count and `nu` is the
#' dispersion: `nu > 1` gives less variance than the Poisson, and
#' `nu < 1` more. The density is `RTMB::dcompois2()`, so
#' `exp(Intercept)` under the log link is the expected count. This is
#' not brms's parameterization. brms 2.23.0's `com_poisson()` fits `mu`
#' as the rate `lambda` of the density `lambda^y / (y!)^nu`, which is
#' not the mean unless `nu = 1`. On counts simulated with `lambda = 9`
#' and `nu = 2` (n = 400), frmtmb's `exp(Intercept)` is 2.78, the sample
#' mean, and brms's is 7.30. The coefficients of the two fits are not
#' comparable one by one (`dev/remeasure-findings.md`, R6).
#'
#' @section Links:
#' Every constructor takes `link` for the mean and a `link_<dpar>` for
"))

# ---- versions -----------------------------------------------------------
vfiles <- c("DESCRIPTION", "codemeta.json",
            "extensions/frmtmb.sample/DESCRIPTION",
            "extensions/frmtmb.eam/DESCRIPTION",
            "extensions/frmtmb.latent/DESCRIPTION",
            "extensions/frmtmb.spline/DESCRIPTION")
for (f in vfiles) from_head(f)
edit_file("DESCRIPTION", "Version: 0.70.1\n", "Version: 0.71.0\n")
edit_file("codemeta.json", '"version": "0.70.1",', '"version": "0.71.0",')
# frmtmb.sample reads the new sampling-API exports (lane postfix) and
# the GEV chart (lane famuni)
edit_file("extensions/frmtmb.sample/DESCRIPTION", "Version: 0.18.0",
          "Version: 0.19.0")
edit_file("extensions/frmtmb.sample/DESCRIPTION",
          "frmtmb (>= 0.70.0)", "frmtmb (>= 0.71.0)")
# frmtmb.eam: a breaking wiener() argument order (lane famuni); its
# code uses no new core function, so the floor stays (measured:
# dev/rel071-floor.sh)
edit_file("extensions/frmtmb.eam/DESCRIPTION", "Version: 0.11.3",
          "Version: 0.12.0")
# frmtmb.latent: a NEWS-only Breaking entry (lane famsimplex), no code
edit_file("extensions/frmtmb.latent/DESCRIPTION", "Version: 0.6.0",
          "Version: 0.6.1")
# frmtmb.spline: lane illse's curve SEs through the factor. No floor:
# the new version reads the factor only where frm_lp_basis() returns
# one, so on 0.70.0 it behaves as 0.10.1 (lane illse's
# dev/illse-spl-oldext.R), and a floor here would not stop 0.10.1 from
# meeting the new core, where it refuses such a fit loudly
if (with_illse) {
  edit_file("extensions/frmtmb.spline/DESCRIPTION", "Version: 0.10.1",
            "Version: 0.10.2")
}

# ---- the release harness reads rellib-r9 ------------------------------
for (f in c("dev/release/run-tests.R", "dev/release/run-check.ps1",
            "dev/release/run-docs.ps1", "dev/warnleak-scan.R",
            "dev/brmsport-run.R", "dev/brmsport-guards.R",
            "dev/ciharden-run1.R", "dev/ciharden-scan.sh")) {
  from_head(f)
  x <- rd(f)
  y <- gsub("rellib-r8", "rellib-r9", x, fixed = TRUE)
  if (!identical(x, y)) { wr(y, f); cat("harness:", f, "\n") }
  stopifnot(any(grepl("rellib-r9", rd(f), fixed = TRUE)))
}
# lane remeasure's bin-2 tooling (untracked, copied by the merge)
for (f in c("dev/brmsport-bin2-run.R", "dev/brmsport-bin2-record.sh",
            "dev/brmsport-bin2-all.sh")) {
  x <- rd(f)
  y <- gsub("rellib-r8", "rellib-r9", x, fixed = TRUE)
  if (!identical(x, y)) { wr(y, f); cat("harness:", f, "\n") }
}

# ---- _pkgdown.yml ---------------------------------------------------------
# Found at the consolidation: pkgdown::check_pkgdown() stops on lane
# famuni's new topic ("1 topic missing from index:
# frmtmb-extreme-families", dev/rel071-log/pkgdown.txt); it goes beside
# the other family pages
from_head("_pkgdown.yml")
edit_file("_pkgdown.yml", "  - frmtmb-families\n  - frmtmb-links\n",
          "  - frmtmb-families\n  - frmtmb-extreme-families\n  - frmtmb-links\n")

# ---- bin 2 of the ported suite: the verdict moves -------------------------
# Lane remeasure's verdict scripts (copied untracked by the merge) write
# the manual tables from R; each now passes its table through
# dev/rel071-bin2-moves.R, which says what moved and why.
for (part in c("2t", "2p")) {
  f <- paste0("dev/brmsport-bin2-verdicts-", part, ".R")
  edit_file(f, "V <- do.call(rbind, v)\n", paste0(
    "V <- do.call(rbind, v)\n",
    "# the verdict moves of the 0.71.0 consolidation\n",
    "V <- source(\"dev/rel071-bin2-moves.R\")$value(V, \"", part, "\")\n"))
}

# ---- dev/brms-port/families-coverage.R: generators for the new families ---
# The families instrument (lane remeasure) fits every family brms names
# on data from a generator, and counted the five families of this round
# (lanes famuni and famsimplex) as named but not fitted, "no data
# generator". Each new generator runs on its own seed and restores the
# RNG state, so the data of every other family are as before.
edit_file("dev/brms-port/families-coverage.R",
"bet <- function() rbeta(n, 2 * p01() * 5, 2 * (1 - p01()) * 5)\n",
paste0(
"bet <- function() rbeta(n, 2 * p01() * 5, 2 * (1 - p01()) * 5)
# the families added at 0.71.0 draw on their own seeds, so the data of
# every family before them in the table stay as they were
own_seed <- function(s, expr) {
  old <- get(\".Random.seed\", envir = globalenv())
  on.exit(assign(\".Random.seed\", old, envir = globalenv()))
  set.seed(20261013 + s)
  expr
}
# a composition of three parts with log ratios linear in x
simplex3 <- function() own_seed(4, {
  e <- cbind(0, 0.3 * x, -0.2 + 0.2 * x) + matrix(rnorm(3 * n, 0, 0.3), n)
  y <- exp(e) / rowSums(exp(e))
  colnames(y) <- c(\"a\", \"b\", \"c\")
  y
})
"))
edit_file("dev/brms-port/families-coverage.R",
"  cox = function() rexp(n, exp(0.3 * x)))\n",
paste0(
"  cox = function() rexp(n, exp(0.3 * x)),
  # mean mu = scale * gamma(1 - 1 / nu), nu = 5
  frechet = function() own_seed(1, exp(0.5 + 0.3 * x) / gamma(1 - 1 / 5) *
                                     (-log(runif(n)))^(-1 / 5)),
  gen_extreme_value = function() own_seed(2, {
    xi <- 0.1
    1 + 0.3 * x + ((-log(runif(n)))^(-xi) - 1) / xi
  }),
  # P(Y >= y) = q^(y^1.5), q = P(Y >= 1), brms's mu
  discrete_weibull = function() own_seed(3, {
    q <- plogis(1 + 0.3 * x)
    ceiling((log(runif(n)) / log(q))^(1 / 1.5)) - 1
  }))
"))
edit_file("dev/brms-port/families-coverage.R",
"    } else {\n      if (is.null(gen[[fam_name]])) stop(\"no data generator\")\n",
paste0(
"    } else if (fam_name %in% c(\"dirichlet\", \"logistic_normal\")) {
      d <- data.frame(x = x)
      d$Y <- simplex3()
      f <- suppressWarnings(frm(bf(Y ~ x),
                                family = frmtmb::brmsfamily(fam_name),
                                data = d))
    } else {
      if (is.null(gen[[fam_name]])) stop(\"no data generator\")
"))

# ---- dev/upstream-bugs.md: the brms entries the lanes drafted ----------
from_head("dev/upstream-bugs.md")
ub <- rd("dev/upstream-bugs.md")
at <- which(ub == "## RTMB")
stopifnot(length(at) == 1L, ub[at - 1L] == "")
ub <- append(ub, c(rd("dev/rel071-upstream-brms.md"), ""), after = at - 1L)
wr(ub, "dev/upstream-bugs.md")
cat("upstream-bugs: brms-24 to brms-27\n")

# ---- dev/test-backlog.md ----------------------------------------------
# The HEAD blob with each item this round closed marked in place, then
# one "Filed at the 0.71.0 release" section before "## Reference": the
# consolidation's own text (dev/rel071-backlog-section.md), then lane
# famuni's section as the lane wrote it, read from its tree.
bl <- sub("\r$", "",
          system2("git", c("show", "HEAD:dev/test-backlog.md"),
                  stdout = TRUE))
close_item <- function(x, start, note) {
  i <- which(startsWith(x, start))
  if (length(i) != 1L) stop("backlog item not unique: ", start)
  j <- i + 1L
  while (j <= length(x) && (startsWith(x[j], "  ") ||
           (x[j] == "" && j < length(x) && startsWith(x[j + 1L], "  ")))) {
    j <- j + 1L
  }
  wrapped <- strwrap(note, width = 78, initial = "  ", prefix = "  ")
  append(x, wrapped, after = j - 1L)
}
cl <- list(
  c("- **The wide high-R^2 design, seed 1, stops 9107.66 below lmer at",
    "CLOSED at 0.71.0 (lane optstart): mu starts at least squares and sigma at the spread it leaves; seeds 1 to 3 end within 4e-8 of lmer (seed 1 at code 1 with a warning, the convergence-warning policy's)."),
  c("- **A group sd 1000 times the residual sd gives \"false convergence\"",
    "CLOSED at 0.71.0 (lane optstart): a restart that neither gains nor moves keeps the converged run's code (restart_result()); 20 of 20 at code 0."),
  c("- **Read every `||` term as brms does** (moparity's r2, in its section",
    "CLOSED at 0.71.0 (lane dblbar): a `||` term is `(x | gr(g, cor = FALSE))`, one block coded with the formula's contrasts, column for column brms 2.23.0's on 29 spellings."),
  c("- **r2. Read every `||` term as brms does.** A `||` term with `mo()`",
    "CLOSED at 0.71.0 (lane dblbar), as above.")
)
if (with_illse) cl <- c(cl, source("dev/rel071-backlog-illse.R")$value)
for (c2 in cl) bl <- close_item(bl, c2[1], c2[2])
lane_sec <- function(l, head) {
  x <- rd(file.path(root, paste0("frmtmb-wt-", l), "dev/test-backlog.md"))
  i <- which(x == head)
  r <- which(x == "## Reference")
  stopifnot(length(i) == 1L, length(r) == 1L, i < r)
  x[i:(r - 1L)]
}
sec <- c(rd("dev/rel071-backlog-section.md"), "",
         lane_sec("famuni", "### Filed by lane famuni (round opened 2026-10-13)"))
if (with_illse && file.exists("dev/rel071-backlog-illse.md")) {
  sec <- c(sec, "", rd("dev/rel071-backlog-illse.md"))
}
ref <- which(bl == "## Reference")
stopifnot(length(ref) == 1L)
bl <- append(bl, sec, after = ref - 1L)
bl <- bl[!(c(FALSE, bl[-1] == "" & bl[-length(bl)] == ""))]
wr(bl, "dev/test-backlog.md")
cat("backlog: closed", length(cl), "items in place\n")
cat("RESOLVE DONE\n")
