# Lane famend: the new tests' assertions over 20 data sets perturbed by
# a few ulps, y * (1 + 2e-15 * k), so that none rests on one rounding.
#   Rscript dev/famend-perturb.R <lib> [k range]
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
ks <- if (length(a) > 1) eval(parse(text = a[2])) else 1:20
cat("lib", dirname(find.package("frmtmb")), "| R", R.home(),
    "| core type", Sys.getenv("OPENBLAS_CORETYPE"),
    as.character(packageVersion("frmtmb")), "| BLAS",
    extSoftVersion()[["BLAS"]], "\n")
conds <- function(expr) {
  w <- list(); m <- list()
  v <- withCallingHandlers(expr, warning = function(x) {
    w[[length(w) + 1L]] <<- x; invokeRestart("muffleWarning")
  }, message = function(x) {
    m[[length(m) + 1L]] <<- x; invokeRestart("muffleMessage")
  })
  list(value = v, warnings = w, messages = m)
}
has <- function(l, cls) any(vapply(l, function(x) inherits(x, cls), NA))
txt <- function(l) vapply(l, conditionMessage, "")
pert <- function(y, k) y * (1 + 2e-15 * k)
fr_data <- function(seed, n = 300, nu0 = 1.1) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n))
  d$y <- exp(0.5 + 0.3 * d$x) / gamma(1 - 1 / nu0) * rexp(n)^(-1 / nu0)
  d
}
mk <- function(seed, n = 300) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n))
  cls <- rbinom(n, 1, 0.4)
  lat <- ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) + rlogis(n)
  d$y <- 1L + (lat > -1 + 0.1 * d$z) + (lat > 0) + (lat > 1 - 0.1 * d$z)
  d
}
rg <- function(n, mu, s, xi) mu + s * ((-log(runif(n)))^(-xi) - 1) / xi
tally <- list()
add <- function(nm, ok) tally[[nm]] <<- c(tally[[nm]], isTRUE(ok))
d9031 <- fr_data(9031)
d9001 <- fr_data(9001)
c3 <- mk(3)
c4 <- mk(4)
set.seed(4)
g63 <- data.frame(y = rg(300, 0, 1, -0.63))
set.seed(4)
g45 <- data.frame(y = rg(300, 0, 1, -0.45))
set.seed(4)
g89 <- data.frame(y = rg(300, 0, 1, -0.89))
set.seed(4)
dn <- data.frame(x = rnorm(300))
dn$y <- rpois(300, exp(1 + 0.3 * dn$x))
for (k in ks) {
  d <- d9031
  d$y <- pert(d$y, k)
  r <- conds(frm(y ~ x, family = frechet(), data = d))
  add("frechet 9031 warns edge", has(r$warnings, "frmtmb_frechet_nu_edge"))
  add("frechet 9031 edge warning alone, state at the stop",
      length(r$warnings) == 1L &&
        identical(unname(r$value$obj$env$last.par.best),
                  unname(r$value$opt$par)) &&
        identical(r$value$obj$env$value.best, r$value$opt$objective))
  d <- d9001
  d$y <- pert(d$y, k)
  r <- conds(frm(y ~ x, family = frechet(), data = d))
  lp <- r$value$frame$linpreds[[ns$linpred_key(
    r$value$spec$responses[[1]]$resp_name, "nu")]]
  add("frechet 9001 silent, gain < 0",
      !length(r$warnings) && !length(r$messages) &&
        ns$frechet_ridge_gain(r$value, r$value$spec$responses[[1]]$resp_name,
                              lp) < 0)
  d <- g63
  d$y <- pert(d$y, k)
  r <- conds(frm(y ~ 1, family = gen_extreme_value(), data = d))
  add("gev -0.63 message", has(r$messages, "frmtmb_gev_xi_nonregular") &&
        !length(r$warnings))
  d <- g45
  d$y <- pert(d$y, k)
  r <- conds(frm(y ~ 1, family = gen_extreme_value(), data = d))
  add("gev -0.45 silent", !length(r$messages) && !length(r$warnings))
  for (s in c(3, 4)) {
    d <- if (s == 3) c3 else c4
    d$z <- pert(d$z, k)
    r <- conds(frm(bf(y ~ x + cs(z)), family = mixture(cumulative(), sratio()),
                   data = d))
    wall <- any(grepl("cs() effect has closed a category", txt(r$warnings),
                      fixed = TRUE))
    add(paste0("csmix seed ", s, if (s == 3) " wall" else " no wall"),
        if (s == 3) wall else !wall)
  }
  d <- dn
  d$x <- pert(d$x, k)
  r1 <- conds(frm(y ~ x, family = negbinomial(), data = d))
  r2 <- conds(frm(y ~ x, family = brmsfamily("negbinomial2"), data = d))
  lim <- function(r) has(r$messages, "frmtmb_limit_fit") ||
    has(r$messages, "frmtmb_boundary_fit")
  # the fix round's rule: the limit message and no warning, or a code-1
  # stop with the convergence warnings alone
  ok_end <- function(r) {
    w <- txt(r$warnings)
    if (lim(r)) return(!length(w))
    conv <- grepl("did not report convergence", w, fixed = TRUE) |
      grepl("Large maximum absolute gradient", w, fixed = TRUE)
    identical(r$value$opt$convergence, 1L) && length(w) && all(conv) &&
      any(grepl("did not report convergence", w, fixed = TRUE))
  }
  l1 <- as.numeric(logLik(r1$value))
  l2 <- as.numeric(logLik(r2$value))
  add("nb/nb2 both limit, logLik agree (round 1 rule)",
      lim(r1) && lim(r2) && abs(l1 - l2) / abs(l1) < 1e-8)
  add("nb/nb2 valid endings, logLik agree (fix rule)",
      ok_end(r1) && ok_end(r2) && abs(l1 - l2) / abs(l1) < 1e-8)
  add("nb/nb2 stops at code 1 (of 2 per k)",
      sum(c(r1$value$opt$convergence, r2$value$opt$convergence) == 1L) > 0)
  d <- g89
  d$y <- pert(d$y, k)
  r <- conds(frm(y ~ 1, family = gen_extreme_value(), data = d))
  add("gev -0.89 message", has(r$messages, "frmtmb_gev_xi_nonregular") &&
        !length(r$warnings))
}
for (nm in names(tally)) {
  cat(sprintf("PERTURB %-34s %d of %d\n", nm, sum(tally[[nm]]),
              length(tally[[nm]])))
}
cat("END\n")
