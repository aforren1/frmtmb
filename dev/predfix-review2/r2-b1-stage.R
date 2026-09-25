source("C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-review2/r2-helpers.R")
# Which stage raises the default's error on seeds 533, 535, 537 of
# r2-b1-sep2.R (b)? Count entries/exits of the pre-fit, and rerun with
# verbose = TRUE to see the stage lines.
.r2$entered <- 0L
suppressMessages(trace("autoscale_prefit", where = asNamespace("frmtmb"),
  tracer = quote(.r2$entered <- .r2$entered + 1L),
  exit = quote(assign("tpl", returnValue(), envir = .r2)), print = FALSE))
for (seed in c(533, 535, 537)) {
  set.seed(seed)
  g <- factor(rep(1:20, each = 12))
  d2 <- data.frame(g, x = rnorm(240) * 0.03, z = rnorm(240))
  d2$yb <- as.integer(d2$z > 0)
  .r2$entered <- 0L; .r2$tpl <- NULL
  e <- tryCatch(frm(yb ~ z + x + (1 + x | g), family = bernoulli(), data = d2),
                error = function(e) e)
  cat(sprintf("seed %d: pre-fit entered %d, returned a template %s, result %s\n",
              seed, .r2$entered, !is.null(.r2$tpl),
              if (inherits(e, "error")) "ERROR" else "fit"))
  if (!is.null(.r2$tpl)) cat("   template beta:", signif(.r2$tpl$beta, 4),
                             " theta:", signif(.r2$tpl$theta, 4), "\n")
  if (seed == 533) {
    v <- capture.output(tryCatch(frm(yb ~ z + x + (1 + x | g), family = bernoulli(),
                      data = d2, verbose = TRUE), error = function(e) cat("ERR\n")),
                      type = "message")
    cat(v, sep = "\n")
  }
}
