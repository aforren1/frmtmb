# Where a dirichlet() prior on a parameter that is not a simplex is
# refused, by spelling, in brms 2.23.0 and in frmtmb. brms is run, not
# recalled: set_prior(), validate_prior(), stancode() and, for the fit,
# brm() with the Stan compile (cache-free; one compile per spelling).
#   Rscript dev/smallfix-dirichlet-probe.R <lib;lib> <brmsfit: yes|no>
# Log: dev/smallfix-out/dirichlet-probe-<arm>.txt
a <- commandArgs(trailingOnly = TRUE)
# class = "sigma" is fitted on bf(y ~ x), where sigma is one number and
# brms has a class "sigma" (review finding 7; the first run used
# bf(y ~ x, sigma ~ 1), where brms has none).
.libPaths(c(strsplit(a[1], ";", fixed = TRUE)[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
fit_brms <- identical(a[2], "yes")
suppressMessages(library(frmtmb))
cat("frmtmb", format(packageVersion("frmtmb")), "from",
    dirname(find.package("frmtmb")), "| brms",
    format(packageVersion("brms")), "\n")
set.seed(1)
d <- data.frame(x = rnorm(40))
d$y <- 1 + d$x + rnorm(40)
spell <- list(
  `class = "Intercept", dpar = "sigma"` = list(class = "Intercept",
                                               dpar = "sigma"),
  `class = "sigma"` = list(class = "sigma", dpar = ""),
  `class = "b"` = list(class = "b", dpar = ""),
  `class = "Intercept"` = list(class = "Intercept", dpar = ""),
  `class = "b", dpar = "sigma"` = list(class = "b", dpar = "sigma"),
  `class = "sd"` = list(class = "sd", dpar = ""))
form <- function(sp) {
  if (identical(sp$dpar, "sigma")) {
    if (identical(sp$class, "b")) return(bf(y ~ x, sigma ~ x))
    return(bf(y ~ x, sigma ~ 1))
  }
  if (identical(sp$class, "sd")) return(bf(y ~ x + (1 | g)))
  bf(y ~ x)
}
d$g <- factor(rep(1:8, 5))
step <- function(lab, expr) {
  r <- tryCatch({ force(expr); "accepted" }, error = function(e) {
    paste("REFUSED:", substr(gsub("\\s+", " ", conditionMessage(e)), 1,
                             160))
  })
  cat(sprintf("    %-22s %s\n", lab, r))
}
for (nm in names(spell)) {
  sp <- spell[[nm]]
  cat("==", nm, "\n")
  cat("  brms\n")
  bp <- NULL
  step("set_prior", bp <- brms::set_prior("dirichlet(2)", class = sp$class,
                                          dpar = sp$dpar))
  bfm <- if (identical(sp$class, "sd")) brms::bf(y ~ x + (1 | g)) else
    if (identical(sp$class, "b") && identical(sp$dpar, "sigma"))
      brms::bf(y ~ x, sigma ~ x) else
    if (identical(sp$dpar, "sigma"))
      brms::bf(y ~ x, sigma ~ 1) else brms::bf(y ~ x)
  if (!is.null(bp)) {
    step("validate_prior", brms::validate_prior(bp, bfm, data = d))
    step("stancode", brms::stancode(bfm, data = d, prior = bp))
    if (fit_brms) {
      step("brm (compile)", suppressWarnings(suppressMessages(
        brms::brm(bfm, data = d, prior = bp, chains = 0, refresh = 0,
                  silent = 2, backend = "rstan"))))
    }
  }
  cat("  frmtmb\n")
  fp <- NULL
  step("set_prior", fp <- frmtmb::set_prior("dirichlet(2)", class = sp$class,
                                            dpar = sp$dpar))
  if (!is.null(fp)) {
    step("validate_prior", frmtmb::validate_prior(fp, form(sp), data = d))
    step("frm", frm(form(sp), family = gaussian(), data = d, prior = fp))
  }
}
