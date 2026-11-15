# Short check of lane priors2's fix round, item 4 (L1, L2): the old
# covariance spelling under class "theta" points to "thetacov"; the
# weights row with an empty density prints flat.
# Rscript priors2-short-l1.R
.libPaths(c("C:/Users/adf44/source/r/wt-priors2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
set.seed(1)
d <- data.frame(x = rnorm(150), g = factor(rep(1:15, 10)),
                t = factor(rep(1:5, 30)), tt = rep(1:10, 15))
d$y <- 1 + 0.4 * d$x + rnorm(15)[d$g] * 0.5 + rnorm(150, sd = 0.7)
k <- rbinom(150, 1, 0.4)
d$ym <- ifelse(k == 1, rnorm(150, -1), rnorm(150, 2))
try1 <- function(lab, expr) {
  r <- tryCatch(suppressWarnings(suppressMessages(expr)),
                error = function(e) paste("ERROR", conditionMessage(e)))
  if (inherits(r, "frmtmb_fit") || is.list(r) && !is.null(r$opt)) {
    r <- paste("fit, code", r$opt$convergence)
  }
  cat(lab, ":\n  ", substr(paste(format(r), collapse = " | "), 1, 260), "\n",
      sep = "")
}
for (cf in c("theta_4", "thetaac_1", "thetar_2")) {
  for (form in c("dens", "bound", "none")) {
    mk <- function() switch(form,
      dens = set_prior("normal(0, 1)", class = "theta", coef = cf),
      bound = set_prior("", class = "theta", coef = cf, lb = -1, ub = 1),
      none = set_prior("", class = "theta", coef = cf))
    try1(paste(cf, form, "set_prior"), mk())
    # the same row built past set_prior()'s check, as a brms table would
    # carry it
    p <- tryCatch({
      p <- set_prior("normal(0, 1)", class = "b")
      p$class <- "theta"
      p$coef <- cf
      if (form != "dens") p$prior <- ""
      if (form == "bound") {
        p$lb <- "-1"
        p$ub <- "1"
      }
      p
    }, error = function(e) {
      cat(cf, " ", form, " edited row:\n  ERROR ",
          substr(conditionMessage(e), 1, 260), "\n", sep = "")
      NULL
    })
    if (is.null(p)) next
    try1(paste(cf, form, "frm toep, row built by hand"),
         frm(y ~ toep(0 + t | g), data = d, prior = p))
  }
}
print(set_prior("", class = "theta", coef = "a"))
print(set_prior("", class = "theta"))
# L2: the refusal's wording
set.seed(1)
n <- 120
e <- data.frame(z = rnorm(n), g = factor(rep(1:10, each = 12)))
e$yneg <- 1 + e$z * 0.3 + rnorm(10)[e$g] * 0.3 + rnorm(n, sd = 0.3) -
  0.8 * e$z
r <- tryCatch(suppressWarnings(frm(yneg ~ z, data = e,
                                   prior = set_prior("gamma(0.5, 1)",
                                                     class = "b", coef = "z",
                                                     lb = 0))),
              error = conditionMessage)
cat("L2:", r, "\n")
