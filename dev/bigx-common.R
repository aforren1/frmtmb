# Shared setup for lane bigx's scripts. BIGX_ARM = "before" reads the
# base build (rellib-r18); "after" puts the lane's library first.
bigx_arm <- Sys.getenv("BIGX_ARM", "after")
bigx_libs <- c(
  if (bigx_arm == "after") "C:/Users/adf44/source/r/wt-bigx-lib",
  "C:/Users/adf44/source/r/rellib-r18",
  "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(bigx_libs)
suppressPackageStartupMessages(library(frmtmb))
cat("arm", bigx_arm, "frmtmb", format(packageVersion("frmtmb")),
    find.package("frmtmb"), "\n")

# A fit with its warnings caught (class and first 70 characters).
bigx_fit <- function(...) {
  w <- character(0)
  wc <- character(0)
  f <- withCallingHandlers(
    tryCatch(frm(...), error = function(e) e),
    warning = function(x) {
      w <<- c(w, substr(conditionMessage(x), 1, 70))
      wc <<- c(wc, paste(setdiff(class(x), c("warning", "condition")),
                         collapse = "/"))
      invokeRestart("muffleWarning")
    },
    message = function(m) invokeRestart("muffleMessage"))
  list(f = f, w = w, wc = wc, err = inherits(f, "error"))
}

# SEs of every outer parameter, in fit$opt$par order.
bigx_se <- function(f) {
  v <- tryCatch(suppressWarnings(frmtmb:::sdr_of(f)$cov.fixed),
                error = function(e) NULL)
  if (is.null(v)) return(rep(NA_real_, length(f$opt$par)))
  sqrt(diag(v))
}

# The data of the grid: family, offset m, data seed, with or without a
# group term. z is the centered covariate; w = z + m.
bigx_data <- function(family, m, seed, re = TRUE, n = 150L, G = 15L) {
  set.seed(seed)
  g <- factor(rep(seq_len(G), length.out = n))
  z <- rnorm(n)
  u <- if (re) rnorm(G, 0, 0.7)[g] else 0
  eta <- 0.5 + 0.4 * z + u
  y <- switch(family,
    poisson = rpois(n, exp(eta)),
    binomial = rbinom(n, 1, plogis(eta)),
    gaussian = eta + rnorm(n, 0, 0.5),
    Gamma = rgamma(n, shape = 3, rate = 3 / exp(eta)))
  data.frame(y = y, w = z + m, z = z, g = g)
}

bigx_family <- function(family) {
  switch(family, poisson = poisson(), binomial = bernoulli(),
         gaussian = gaussian(), Gamma = Gamma(link = "log"))
}
