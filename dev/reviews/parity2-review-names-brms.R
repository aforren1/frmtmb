# Review of lane parity2: brms 2.23.0 side of parity2-review-names.R
# (5315 intercept place, 1560 nlf in a body, 5300 err names on more
# latent shapes), fixed_param fits at init 0.
#   Rscript dev/reviews/parity2-review-names-brms.R
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressMessages(library(brms))
options(mc.cores = 1)
cat("brms", as.character(packageVersion("brms")), "\n")
source("C:/Users/adf44/source/r/frmtmb-wt-parity2/dev/reviews/parity2-review-data.R")
try_ <- function(label, expr) {
  m <- character()
  r <- withCallingHandlers(
    tryCatch(expr, error = function(e) e),
    warning = function(cnd) invokeRestart("muffleWarning"),
    message = function(cnd) {
      m <<- c(m, conditionMessage(cnd))
      invokeRestart("muffleMessage")
    })
  if (inherits(r, "error")) cat(label, ": ERROR ", conditionMessage(r), "\n",
                                sep = "")
  for (x in unique(m)) {
    if (grepl("Replacing|Compiling|Start sampling", x)) {
      cat("   message: ", substr(gsub("\n", " ", x), 1, 160), "\n", sep = "")
    }
  }
  invisible(r)
}
bfix <- function(f, fam, data = rd, ...) {
  suppressWarnings(brm(f, family = fam, data = data,
                       algorithm = "fixed_param", chains = 1, iter = 20,
                       warmup = 0, init = 0, refresh = 0, seed = 1,
                       silent = 2, ...))
}
vars <- function(label, f, fam, data = rd, pat = "^(b_|bs|sd|sigma|Intercept|cor|r_)") {
  fit <- try_(label, bfix(f, fam, data))
  if (!inherits(fit, "error")) {
    cat(label, ": ", paste(grep(pat, variables(fit), value = TRUE),
                           collapse = " "), "\n", sep = "")
  }
  invisible(fit)
}
rd$yo <- cut(rd$y, quantile(rd$y, c(0, 0.3, 0.6, 1)), include.lowest = TRUE,
             labels = FALSE)

cat("\n#### 5315 intercept place\n")
vars("center_false", bf(y ~ x, center = FALSE), gaussian())
vars("both_0I", bf(y ~ 0 + Intercept + x, sigma ~ 0 + Intercept + z),
     gaussian())
vars("ord_cf", bf(yo ~ x, center = FALSE), cumulative())
vars("mv_0I", bf(y ~ 0 + Intercept + x) + bf(y2 ~ z) + set_rescor(FALSE),
     gaussian())
vars("re_0I", bf(y ~ 0 + Intercept + x + (1 | g), sigma ~ z), gaussian())
vars("nl_sig0I", bf(y ~ a * exp(b * x), a ~ 1, b ~ 1,
                    sigma ~ 0 + Intercept + z, nl = TRUE), gaussian())

cat("\n#### 1560 nlf() replacing a parameter used in another's body\n")
f1 <- try_("nlf over lf a", bf(y ~ a + b * x, a ~ 1 + w, b ~ 1, nl = TRUE) +
             nlf(a ~ exp(c0)) + lf(c0 ~ 1 + z))
if (!inherits(f1, "error")) {
  print(f1)
  vars("nlf over lf a fit", f1, gaussian(), pat = "^b_")
}
f3 <- try_("lf over mu", bf(y ~ x) + lf(mu ~ z))
if (!inherits(f3, "error")) print(f3)

cat("\n#### 5300 err names: mv latent and ma\n")
set.seed(4)
dp <- data.frame(t = rep(1:8, 6), g = factor(rep(1:6, each = 8)),
                 x = rnorm(48))
dp$y1 <- rpois(48, exp(0.5 + 0.3 * dp$x))
dp$y2 <- rpois(48, exp(0.2 + 0.1 * dp$x))
fe <- try_("mv latent", bfix(bf(y1 ~ x + ar(time = t, gr = g)) +
                               bf(y2 ~ x) + set_rescor(FALSE),
                             poisson(), data = dp))
if (!inherits(fe, "error")) {
  v <- variables(fe)
  cat("mv: n", length(v), " first err at", grep("^err", v)[1], " names",
      paste(head(grep("^(err|zerr|sderr|ar)", v, value = TRUE), 6),
            collapse = " "), "... last", tail(v, 1), "\n")
}
fma <- try_("ma latent", bfix(bf(y1 ~ x + ma(time = t, gr = g)), poisson(),
                              data = dp))
if (!inherits(fma, "error")) {
  v <- variables(fma)
  cat("ma: ", paste(v[!grepl("^err\\[([2-9]|[1-4][0-9])\\]", v)],
                    collapse = " "), "\n")
}
fcv <- try_("ar cov latent", bfix(bf(y1 ~ x + ar(time = t, gr = g, cov = TRUE)),
                                  poisson(), data = dp))
if (!inherits(fcv, "error")) {
  v <- variables(fcv)
  cat("ar cov: ", paste(v[!grepl("^err\\[([2-9]|[1-4][0-9])\\]", v)],
                        collapse = " "), "\n")
}
