# Reviewer probe 1: what prior = "flat" and "" specifications resolve to
# on models the lane did not probe. Usage:
# Rscript dev/review-reslot-p1.R before|after. Seed 1.
source("C:/Users/adf44/source/r/frmtmb-wt-reslot/dev/reslot-common.R")
suppressPackageStartupMessages(library(frmtmb.sample))
dd <- re_data(1L)
set.seed(1)
dd$w <- rnorm(nrow(dd))
dd$ord <- cut(dd$y, 4, labels = FALSE)
dd$yp <- abs(dd$y) + 0.1
dd$cnt <- rpois(nrow(dd), 2)
dd$yb <- plogis(dd$y / 3)
dd$y3 <- 0.3 * dd$y + rnorm(nrow(dd))
dd$yn <- 2 * exp(0.3 * dd$x) + rnorm(nrow(dd), 0, 0.3)
show_e <- function(ri) {
  for (e in ri$entries) {
    cat(sprintf("  entry %-8s %-10s %-8s %-9s link=%s eta=%s flat=%s innov=%s\n",
                e$comp, paste(e$idx, collapse = ","), e$dist$kind,
                e$scale, e$link$name %||% "-",
                paste(e$dist$eta, collapse = ","),
                isTRUE(e$dist$flat), !is.null(e$innov)))
  }
  if (length(ri$lower)) cat("  lower:", paste(names(ri$lower), ri$lower,
                                              sep = "=", collapse = " "), "\n")
  if (length(ri$upper)) cat("  upper:", paste(names(ri$upper), ri$upper,
                                              sep = "=", collapse = " "), "\n")
}
show <- function(label, f, prior = "flat", ...) {
  cat("--", label, "\n")
  fit <- tryCatch(suppressWarnings(suppressMessages(
    frm(f, data = dd, dry_run = "objective", ...))), error = function(e) e)
  if (inherits(fit, "error")) {
    cat("  FIT ERROR:", conditionMessage(fit), "\n")
    return(invisible())
  }
  r <- tryCatch(suppressWarnings(suppressMessages(
    frmtmb.sample:::sample_resolve_priors(fit, prior))),
    error = function(e) e)
  if (inherits(r, "error")) {
    cat("  ERROR:", conditionMessage(r), "\n")
    return(invisible())
  }
  if (is.null(r$ri)) cat("  no entries\n") else show_e(r$ri)
  if (!is.null(r$effective) && length(r$effective)) {
    print(as.data.frame(r$effective)[, c("prior", "class", "coef", "group",
                                         "dpar", "nlpar", "lb", "ub")],
          row.names = FALSE)
  }
}
show("student", bf(y ~ x), family = student())
show("sigma ~ x", bf(y ~ x, sigma ~ x))
show("beta phi", bf(yb ~ x), family = Beta())
show("zip", bf(cnt ~ x), family = zero_inflated_poisson())
show("negbinomial", bf(cnt ~ x), family = negbinomial())
show("lognormal", bf(yp ~ x), family = lognormal())
show("skew_normal", bf(y ~ x), family = skew_normal())
show("exgaussian", bf(y ~ x), family = exgaussian())
show("weibull", bf(yp ~ x), family = weibull())
show("hurdle_poisson", bf(cnt ~ x), family = hurdle_poisson())
show("cumulative flexible", bf(ord ~ x), family = cumulative())
show("mixture(gaussian, gaussian)", bf(y ~ x),
     family = mixture(gaussian(), gaussian()))
show("nonlinear", bf(yn ~ a * exp(k * x), a ~ 1 + (1 | g), k ~ 1, nl = TRUE))
show("mv rescor 3", bf(y ~ x) + bf(y2 ~ x) + bf(y3 ~ x) + set_rescor(TRUE))
show("(1 + x + z + w | g)", bf(y ~ x + (1 + x + z + w | g)))
show("ar(cov = FALSE)", bf(y ~ x + ar(t, g)))
show("ma(cov = TRUE)", bf(y ~ x + ma(t, g, cov = TRUE)))
show("arma(cov = TRUE)", bf(y ~ x + arma(t, g, cov = TRUE)))
show("unstr", bf(y ~ x + unstr(t, g)))
show("gp(x) + gp(z)", bf(y ~ gp(x) + gp(z)))
show("gp(x, z)", bf(y ~ gp(x, z)))
show("cs(0 + tf | g)", bf(y ~ x + cs(0 + factor(t) | g)))

cat("\n==== '' spellings on frm_sample's route ====\n")
show("cor ''", bf(y ~ x + (1 + x + z | g)), prior = set_prior("", class = "cor"))
show("cor '' group g", bf(y ~ x + (1 + x + z | g)),
     prior = set_prior("", class = "cor", group = "g"))
show("rescor ''", bf(y ~ x) + bf(y2 ~ x) + bf(y3 ~ x) + set_rescor(TRUE),
     prior = set_prior("", class = "rescor"))
show("cortime ''", bf(y ~ x + unstr(t, g)), prior = set_prior("", class = "cortime"))
show("sd ''", bf(y ~ x + (1 + x | g)), prior = set_prior("", class = "sd"))
show("sigma '' under ar(cov = TRUE)", bf(y ~ x + ar(t, g, cov = TRUE)),
     prior = set_prior("", class = "sigma"))
show("Intercept '' cumulative", bf(ord ~ x), family = cumulative(),
     prior = set_prior("", class = "Intercept"))
show("lscale ''", bf(y ~ gp(x)), prior = set_prior("", class = "lscale"))
show("nu '' student", bf(y ~ x), family = student(),
     prior = set_prior("", class = "nu"))
show("nu lb=2 student", bf(y ~ x), family = student(),
     prior = set_prior("", class = "nu", lb = 2))
