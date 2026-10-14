# Lane ci701: siblings of the "singular convergence (7)" escape. Each
# candidate is a degenerate ordinal-mixture fit in a test whose
# allow_warnings() list does not allow "did not report convergence".
# The fit is repeated on data with x moved by a few ulps (seed 4000 + k)
# and every warning is counted.
#
#   Rscript dev/ci701-sib-probe.R <lib|base> <candidate> <K>
#   candidate: cs   test-ordinal-mixture.R:401 (omx_data(20261012))
#              none frmtmb.sample test-ordinal-mixture-draws.R:47
args <- commandArgs(TRUE)
lib <- args[1]; cand <- args[2]; K <- as.integer(args[3])
base <- "C:/Users/adf44/source/r/rellib-r8"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (lib == "base") c(base, ul) else c(lib, base, ul))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), as.character(packageVersion("frmtmb")),
    "\n")
omx_data <- function(seed, n = 400) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n),
                  g = factor(sample(c("a", "b"), n, TRUE)))
  cls <- stats::rbinom(n, 1, stats::plogis(-0.4 + 0.5 * d$z))
  lat <- ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) +
    stats::rlogis(n) / exp(0.3 * d$z * cls)
  d$y <- 1L + (lat > -1.5) + (lat > 0) + (lat > 1.5)
  d$yh <- ifelse(stats::runif(n) < stats::plogis(-1.2 + 0.5 * d$z), 0L, d$y)
  d
}
omd_data <- function(seed, n = 200) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n),
                  g = factor(sample(c("a", "b"), n, TRUE)))
  cls <- stats::rbinom(n, 1, 0.4)
  lat <- ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) + stats::rlogis(n)
  d$y <- 1L + (lat > -1.5) + (lat > 0) + (lat > 1.5)
  d$yh <- ifelse(stats::runif(n) < 0.2, 0L, d$y)
  d
}
d0 <- switch(cand, cs = omx_data(20261012), none = omd_data(71))
fitit <- function(d) {
  switch(cand,
    cs = frm(bf(y ~ x, mu2 ~ cs(x)), family = mixture(cumulative(), sratio()),
             data = d),
    none = {
      st <- "student_t(3, 0, 2.5)"
      nb <- set_prior("normal(0, 2)", class = "b", dpar = "mu1") +
        set_prior("normal(0, 2)", class = "b", dpar = "mu2")
      frm(bf(y ~ x, disc1 ~ 0 + z),
          family = mixture(cumulative(), sratio()), data = d,
          prior = nb + set_prior(st, class = "Intercept", dpar = "mu1") +
            set_prior(st, class = "Intercept", dpar = "mu2"))
    })
}
for (k in 0:K) {
  d <- d0
  if (k > 0) {
    set.seed(4000 + k)
    d$x <- d$x * (1 + sample(0:4, nrow(d), TRUE) * 2^-52)
  }
  w <- character(0)
  fit <- withCallingHandlers(fitit(d), warning = function(cnd) {
    w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning")
  })
  cat(sprintf("RES cand=%s k=%d code=%d msg=%s ll=%.8f nwarn=%d\n", cand, k,
              fit$opt$convergence, gsub(" ", "_", fit$opt$message),
              as.numeric(logLik(fit)), length(w)))
  for (m in unique(w)) cat("   W:", substr(m, 1, 120), "\n")
}
