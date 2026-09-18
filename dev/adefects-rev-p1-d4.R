source("C:/Users/adf44/source/r/frmtmb-wt-adefects/dev/adefects-rev-prelude.R")
suppressPackageStartupMessages(library(frmtmb))

addr <- function(x) {
  strsplit(trimws(utils::capture.output(.Internal(inspect(x)))[1L]),
           "[[:space:]]+")[[1L]][1L]
}

cat("== PIN 3 on the REAL class, which the lane's inverse only proxies\n")
cat("  before registering:",
    is.null(getS3method("$", "frmtmb_fit", optional = TRUE)), "(want TRUE)\n")
registerS3method("$", "frmtmb_fit", function(x, name) {
  .subset2(x, name, exact = FALSE)
})
cat("  after registering :",
    is.null(getS3method("$", "frmtmb_fit", optional = TRUE)),
    "(want FALSE)\n")

cat("\n== PIN 4's window, on the fit the test uses and on a wider one\n")
mk <- function(n, p = 0L, seed = 20260917) {
  set.seed(seed)
  d <- data.frame(g = factor(rep_len(1:6, n)), x = rnorm(n))
  for (k in seq_len(p)) d[[paste0("v", k)]] <- rnorm(n)
  d$y <- d$x + rnorm(n)
  d
}
probe <- function(n, p) {
  d <- mk(n, p)
  rhs <- paste(c("x", if (p) paste0("v", seq_len(p))), collapse = " + ")
  f <- stats::as.formula(paste("y ~", rhs, "+ (1 | g)"))
  fit <- frm(f, d)
  nd <- fit; nd$data <- NULL
  cols <- length(serialize(lapply(fit$data, identity), NULL))
  delta <- length(serialize(fit, NULL)) - length(serialize(nd, NULL))
  cat(sprintf("  n=%-6d p=%-2d cols=%-8d delta=%-8d ratio=%.4f in(0.9,3)=%s\n",
              n, p, cols, delta, delta / cols,
              delta > 0.9 * cols && delta < 3 * cols))
}
for (a in list(c(30, 0), c(30, 4), c(240, 0), c(2000, 0), c(20000, 0),
               c(20000, 8))) probe(a[1], a[2])

cat("\n== the terms-attribute environment claim\n")
d <- mk(30, 0)
fit <- frm(y ~ x + (1 | g), d)
mf <- fit$data
cat("  environment(attr(mf, 'terms')) is:",
    environmentName(environment(attr(mf, "terms"))), "\n")
cat("  bytes of serialize(mf)             :", length(serialize(mf, NULL)),
    "\n")
cat("  bytes of the columns alone         :",
    length(serialize(lapply(mf, identity), NULL)), "\n")
e <- new.env(parent = globalenv())
assign("ballast", runif(2e5), envir = e)
mf2 <- mf
attr(attr(mf2, "terms"), ".Environment") <- e
cat("  with 200k doubles in that env      :",
    length(serialize(mf2, NULL)), "\n")
cat("  -> the env IS serialized with the frame:",
    length(serialize(mf2, NULL)) > length(serialize(mf, NULL)) + 1e6, "\n")
cat("  but the DELTA is immune, because a second reference to the same\n")
cat("     environment is written as a reference:\n")
l1 <- list(mf2); l2 <- list(mf2, mf2)
cat("     one copy:", length(serialize(l1, NULL)),
    " two copies:", length(serialize(l2, NULL)),
    " delta:", length(serialize(l2, NULL)) - length(serialize(l1, NULL)),
    " cols:", length(serialize(lapply(mf, identity), NULL)), "\n")

cat("\n== the memory-sharing pin and its control\n")
cat("  addr(fit$data) == addr(frame$data_frame):",
    identical(addr(fit$data), addr(fit$frame[["data_frame"]])), "\n")
cat("  addr(fit$data) == addr(an equal copy)   :",
    identical(addr(fit$data), addr(data.frame(fit$data))), "\n")
cat("  identical(fit$data, an equal copy)      :",
    identical(fit$data, data.frame(fit$data)), "\n")
