# Review of lane slopescale, attack item 3: autoscale_block_exact() is a
# check at one point (theta = start + 0.25 cos(k)), not a proof. A
# registered structure whose sd is exp(theta) * (1 + exp(5 * (theta - 7)))
# is a valid one-to-one sd parameterization, scale equivariant to 1e-15
# near theta = 0.25 and not at all near theta = 8.5. Does it pass the
# check, and what does the default then do with a slope on x * 1e-4,
# whose sd in theta lands near 8.5? Reference: the same model on x (sd
# near 0.5, where the bump is 1e-17), which is the same maximum.
# Usage: Rscript slopescale-review-blockexact.R
.libPaths(c("C:/Users/adf44/source/r/wt-slopescale-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb", find.package("frmtmb"), "\n")
sdf <- function(theta) exp(theta) * (1 + exp(5 * (theta - 7)))
frmtmb_register_covstruct(
  "mykink",
  lpdf = function(b, theta, blk) {
    d <- blk[["dim"]]
    B <- RTMB::matrix(b, d)
    s <- sdf(theta[seq_len(d)])
    sum(RTMB::dnorm(as.vector(B / s), 0, 1, log = TRUE)) -
      ncol(B) * sum(log(s))
  },
  vcov = function(theta, blk) {
    V <- diag(sdf(theta[seq_len(blk[["dim"]])])^2, blk[["dim"]])
    dimnames(V) <- list(blk[["cnms"]], blk[["cnms"]])
    V
  },
  npar = function(blk) blk[["dim"]],
  sd_idx = function(blk) seq_len(blk[["dim"]]),
  package = "reviewpkg")
res <- list()
for (seed in 1:6) {
  set.seed(seed)
  g <- factor(rep(1:20, each = 15))
  x0 <- rnorm(300)
  y <- 0.5 + rnorm(20, 0, 0.5)[g] * x0 + rnorm(300)
  d1 <- data.frame(g = g, x = x0, y = y)
  ds <- d1
  ds$x <- x0 * 1e-4
  fo <- y ~ 1 + mykink(0 + x | g)
  fr <- frm(fo, data = ds, dry_run = "frame")
  pl <- frmtmb:::autoscale_plan(fr)
  pass <- length(pl) && length(pl[[1]]$z) > 0
  nw <- 0L
  q <- function(expr) withCallingHandlers(expr, warning = function(w) {
    nw <<- nw + 1L
    invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage"))
  f1 <- q(frm(fo, data = d1))
  nw <- 0L
  fs <- q(frm(fo, data = ds))
  nws <- nw
  nw <- 0L
  ff <- q(frm(fo, data = ds, control = frmtmb_control(autoscale = FALSE)))
  res[[seed]] <- data.frame(
    seed = seed, check_pass = pass, autoscaled = isTRUE(fs$autoscaled),
    sd_ref = signif(sdf(f1$opt$par[["theta"]]), 4),
    sd_default = signif(sdf(fs$opt$par[["theta"]]) * 1e-4, 4),
    short_default = signif(as.numeric(logLik(f1)) - as.numeric(logLik(fs)), 3),
    warn_default = nws,
    short_false = signif(as.numeric(logLik(f1)) - as.numeric(logLik(ff)), 3))
}
r <- do.call(rbind, res)
print(r, row.names = FALSE)
cat("\ncheck passed:", sum(r$check_pass), "of", nrow(r),
    "; default short > 1e-4:", sum(r$short_default > 1e-4),
    "(silent", sum(r$short_default > 1e-4 & r$warn_default == 0), ")",
    "; max short", max(r$short_default), "\n")
