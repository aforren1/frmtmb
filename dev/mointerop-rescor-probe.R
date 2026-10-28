## Probe for the bin-2 rescor rows (log_lik:107/109,
## posterior_predict:69/72): brms's prep of tests.log_lik.R and
## tests.posterior_predict.R moved onto a frmtmb set_rescor(TRUE) model,
## and its values against brms's own internal functions.
##   Rscript dev/mointerop-rescor-probe.R <lane|base>
arm <- commandArgs(TRUE)[1]
base <- c("C:/Users/adf44/source/r/rellib-r14",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (identical(arm, "lane")) {
  c("C:/Users/adf44/source/r/wt-mointerop-lib", base)
} else base)
suppressMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
`%||%` <- function(a, b) if (is.null(a)) b else a
src <- "C:/Users/adf44/source/r/frmtmb-wt-mointerop/dev/brmsport-bin2-adapters.R"
source(src)
set.seed(1)
ns <- 10
nvars <- 3
ncols <- 4
nobs <- nvars * ncols
prep <- structure(list(), class = "mvbrmsprep")
Sigma <- array(cov(matrix(rnorm(300), ncol = 3)), dim = c(3, 3, 10))
prep$mvpars <- list(
  Mu = array(rnorm(ns * nobs * nvars), dim = c(ns, nobs, nvars)),
  Sigma = aperm(Sigma, c(3, 1, 2))
)
prep$dpars <- list(nu = rgamma(ns, 5))
prep$ndraws <- ns
prep$data <- list(Y = matrix(rnorm(nobs), ncol = nvars))
for (f in c("gaussian_mv", "student_mv")) {
  i <- if (f == "gaussian_mv") 1L else 2L
  a <- tryCatch(brms_adapt_log_lik(f, i, prep), error = function(e) e)
  b <- get(paste0("log_lik_", f), asNamespace("brms"))(i, prep = prep)
  if (inherits(a, "error")) {
    cat(f, "ERROR:", conditionMessage(a), "\n")
  } else {
    cat(f, "length", length(a), "max abs diff to brms",
        format(max(abs(a - b))), "max rel",
        format(max(abs(a - b) / abs(b))), "\n")
  }
  # every row against brms's
  d <- vapply(seq_len(nrow(prep$data$Y)), function(r) {
    max(abs(brms_adapt_log_lik(f, r, prep) -
              get(paste0("log_lik_", f), asNamespace("brms"))(r, prep = prep)))
  }, 0)
  cat(f, "all", nobs, "rows: max abs diff", format(max(d)), "\n")
}
prep2 <- structure(list(ndraws = ns), class = "mvbrmsprep")
prep2$mvpars <- prep$mvpars
prep2$dpars <- prep$dpars
prep2$data <- list(N = nobs, N_trait = ncols)
for (f in c("gaussian_mv", "student_mv")) {
  i <- if (f == "gaussian_mv") 1L else 2L
  p <- tryCatch(brms_adapt_posterior_predict(f, i, prep2),
                error = function(e) e)
  if (inherits(p, "error")) {
    cat(f, "pp ERROR:", conditionMessage(p), "\n")
  } else {
    cat(f, "pp dim", dim(p), "\n")
  }
}
## distribution of the joint draw against brms's on one draw's Sigma:
## the same prep with every draw's Mu equal, many draws
nb <- 4000
prep3 <- structure(list(ndraws = nb), class = "mvbrmsprep")
mu1 <- c(1, -1, 0.5)
S1 <- prep$mvpars$Sigma[1, , ]
prep3$mvpars <- list(
  Mu = array(rep(mu1, each = nb * nobs), dim = c(nb, nobs, nvars)),
  Sigma = array(rep(S1, each = nb), dim = c(nb, 3, 3)))
prep3$dpars <- list(nu = rep(6, nb))
prep3$data <- list(N = nobs, N_trait = ncols)
for (f in c("gaussian_mv", "student_mv")) {
  set.seed(2)
  p <- brms_adapt_posterior_predict(f, 1L, prep3)
  set.seed(2)
  q <- get(paste0("posterior_predict_", f), asNamespace("brms"))(1L,
                                                                 prep = prep3)
  cat(f, "frmtmb mean", format(colMeans(p), digits = 4), "cov",
      format(c(cov(p)), digits = 4), "\n")
  cat(f, "brms   mean", format(colMeans(q), digits = 4), "cov",
      format(c(cov(q)), digits = 4), "\n")
  cat(f, "target cov", format(c(S1 * if (f == "student_mv") 6 / 4 else 1),
                              digits = 4), "\n")
}
