# Lane support, review finding 1: the two fits of the new test in
# test-families-extreme.R, importance = 500 at xi 1.5 (seed 2) and
# nested quadrature at xi 0.5 (seed 3), with the warnings each gives,
# then the same two on copies of y scaled by (1 + 2e-15 k), k = 1..K,
# to see whether completion depends on rounding.
#
#   Rscript dev/support-isnq-probe.R [K]   (FRMTMB_PORT_LIB picks the build)
lib <- Sys.getenv("FRMTMB_PORT_LIB", "C:/Users/adf44/source/r/rellib-r19")
.libPaths(c(strsplit(lib, ";")[[1]],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(system.file(package = "frmtmb")), "\n")
K <- as.integer(commandArgs(trailingOnly = TRUE)[1])
if (is.na(K)) K <- 0L
rg <- function(n, mu, s, k) mu + s * ((-log(stats::runif(n)))^(-k) - 1) / k
sim <- function(seed, xi) {
  set.seed(seed)
  g <- rep(seq_len(15), each = 8)
  u <- stats::rnorm(15, 0, 0.7)
  d <- data.frame(g = factor(g), y = rg(120, 1 + u[g], 1, xi))
  d$h <- factor(paste(d$g, rep(1:2, length.out = 120)))
  d
}
run <- function(label, expr) {
  w <- character()
  r <- tryCatch(withCallingHandlers(expr, warning = function(c) {
    w <<- c(w, conditionMessage(c))
    invokeRestart("muffleWarning")
  }), error = function(e) e)
  if (inherits(r, "error")) {
    cat(sprintf("%s ERROR: %s\n", label, substr(conditionMessage(r), 1, 160)))
    return(NA_real_)
  }
  ll <- as.numeric(logLik(r))
  cat(sprintf("%s logLik %.12g warnings %d\n", label, ll, length(w)))
  for (m in unique(w)) cat("   W:", substr(m, 1, 200), "\n")
  ll
}
fits <- function(k) {
  di <- sim(2, 1.5)
  dq <- sim(3, 0.5)
  di$y <- di$y * (1 + 2e-15 * k)
  dq$y <- dq$y * (1 + 2e-15 * k)
  run(sprintf("k %2d importance", k),
      frm(y ~ 1 + (1 | g), family = gen_extreme_value(), data = di,
          importance = 500L))
  run(sprintf("k %2d nested quad", k),
      frm(y ~ 1 + (1 | g) + (1 | h), family = gen_extreme_value(),
          data = dq, quadrature = TRUE))
}
for (k in 0:K) fits(k)
