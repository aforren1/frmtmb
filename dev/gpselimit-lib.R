# Lane gpselimit: shared helpers. Arm "lane" reads wt-gpselimit-lib over
# rellib-r22; arm "base" reads rellib-r22 alone (frmtmb 0.83.0 at
# c6383527).
GS_WT <- "C:/Users/adf44/source/r/frmtmb-wt-gpselimit"
GS_LIB <- "C:/Users/adf44/source/r/wt-gpselimit-lib"
GS_BASE <- "C:/Users/adf44/source/r/rellib-r22"
GS_USER <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
GS_OUT <- file.path(GS_WT, "dev/gpselimit-out")
dir.create(GS_OUT, showWarnings = FALSE)
gs_libs <- function(arm = c("lane", "base", "r21", "gpfield")) {
  arm <- match.arg(arm)
  # r21: 0.82.0 alone; gpfield: the gpfield lane's build over it, which
  # the gpfield review measured (both read only)
  r21 <- "C:/Users/adf44/source/r/rellib-r21"
  lp <- switch(arm, lane = c(GS_LIB, GS_BASE), base = GS_BASE, r21 = r21,
               gpfield = c("C:/Users/adf44/source/r/wt-gpfield-lib", r21))
  .libPaths(c(lp, GS_USER))
  Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
             FRMTMB_STAN_CACHE = file.path(GS_WT, "dev/stan-cache"),
             NOT_CRAN = "true")
  suppressPackageStartupMessages(library(frmtmb))
  cat("lib:", dirname(system.file(package = "frmtmb")), "frmtmb",
      as.character(packageVersion("frmtmb")), "arm", arm, "\n")
  invisible(arm)
}
# A fit with its warnings and messages kept (class and text), never an
# error thrown.
gs_fit <- function(...) {
  w <- character(0)
  m <- character(0)
  wc <- character(0)
  mc <- character(0)
  f <- tryCatch(withCallingHandlers(frm(...), warning = function(e) {
    w <<- c(w, conditionMessage(e))
    wc <<- c(wc, class(e)[1])
    invokeRestart("muffleWarning")
  }, message = function(e) {
    m <<- c(m, conditionMessage(e))
    mc <<- c(mc, class(e)[1])
    invokeRestart("muffleMessage")
  }), error = function(e) e)
  list(fit = f, warn = w, wclass = wc, msg = m, mclass = mc)
}
# Every condition of an expression, kept, never thrown on.
gs_quiet <- function(expr) {
  w <- character(0)
  v <- withCallingHandlers(expr, warning = function(e) {
    w <<- c(w, paste0(class(e)[1], ": ", conditionMessage(e)))
    invokeRestart("muffleWarning")
  }, message = function(e) invokeRestart("muffleMessage"))
  attr(v, "conds") <- w
  v
}
gs_put <- function(name, f) {
  ns <- asNamespace("frmtmb")
  environment(f) <- ns
  unlockBinding(name, ns); assign(name, f, envir = ns); lockBinding(name, ns)
}
# The data of row `row` of the gpfield review's 30-fit sweep
# (dev/reviews/gpfield-review-basis.R): f from an exact squared exponential
# GP, noise sd 0.3.
gs_sweep_sets <- list(c(n = 60, ls0 = 0.7), c(n = 100, ls0 = 1.5),
                      c(n = 100, ls0 = 3), c(n = 200, ls0 = 1),
                      c(n = 50, ls0 = 2))
gs_sweep_ks <- c(5, 8, 10, 15, 20, 30)
gs_sweep_row <- function(row, jitter = 0) {
  si <- (row - 1L) %/% length(gs_sweep_ks) + 1L
  k <- gs_sweep_ks[(row - 1L) %% length(gs_sweep_ks) + 1L]
  n <- gs_sweep_sets[[si]][["n"]]
  ls0 <- gs_sweep_sets[[si]][["ls0"]]
  set.seed(9600 + row)
  x <- sort(runif(n, 0, 10))
  K <- exp(-outer(x, x, "-")^2 / (2 * ls0^2)) + diag(1e-8, n)
  f <- as.vector(t(chol(K)) %*% rnorm(n))
  y <- f + 0.3 * rnorm(n)
  if (jitter) y <- y * (1 + 2e-15 * jitter)
  list(d = data.frame(x = x, y = y), k = k, n = n, ls0 = ls0)
}
# The standard errors of every outer parameter and the lost reasons, as
# the API reports them, with every condition kept.
gs_ses <- function(fit) {
  sdr <- gs_quiet(frmtmb:::sdr_of(fit))
  se <- sqrt(diag(sdr$cov.fixed))
  list(se = se, lost = sdr$se_lost, conds = attr(sdr, "conds"))
}
