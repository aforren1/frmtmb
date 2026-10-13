# Every covariance structure's log density with each log sd (and all of
# them at once) put where exp() underflows, with the field off zero.
# A finite value or -Inf is a density a sampler rejects; +Inf or NaN is
# the defect (dev/test-backlog.md, "A dense block with several scales").
# Usage: Rscript dev/moparity-floor-sweep.R [lib]
args <- commandArgs(TRUE)
lib <- if (length(args) && nzchar(args[[1]])) args[[1]] else NULL
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
set.seed(3)
n_g <- 12
dd <- data.frame(g = factor(rep(seq_len(n_g), each = 5)),
                 tim = factor(rep(1:5, n_g)),
                 tn = num_factor(rep(c(0, 0.3, 0.4, 1.1, 1.5), n_g)),
                 x = stats::runif(5 * n_g), z = stats::rnorm(5 * n_g))
dd$y <- stats::rnorm(5 * n_g)
A <- diag(n_g) * 0.5 + 0.5
dimnames(A) <- list(levels(dd$g), levels(dd$g))
Q <- solve(A)
forms <- list(
  us1 = y ~ 1 + (1 | g),
  us2 = y ~ 1 + (1 + x | g),
  us3 = y ~ 1 + (1 + x + z | g),
  diag2 = y ~ 1 + diag(1 + x | g),
  homdiag = y ~ 1 + homdiag(1 + x | g),
  ust1 = y ~ 1 + (1 | gr(g, dist = "student")),
  ust2 = y ~ 1 + (1 + x | gr(g, dist = "student")),
  diagt2 = y ~ 1 + (1 + x || gr(g, dist = "student")),
  ar1 = y ~ 1 + ar1(tim + 0 | g),
  hetar1 = y ~ 1 + hetar1(tim + 0 | g),
  cs = y ~ 1 + cs(tim + 0 | g),
  homcs = y ~ 1 + homcs(tim + 0 | g),
  toep = y ~ 1 + toep(tim + 0 | g),
  homtoep = y ~ 1 + homtoep(tim + 0 | g),
  ou = y ~ 1 + ou(tn + 0 | g),
  grcov1 = y ~ 1 + (1 | gr(g, cov = A)),
  grcov2 = y ~ 1 + (1 + x | gr(g, cov = A)),
  grprec1 = y ~ 1 + (1 | gr(g, prec = Q)),
  grprec2 = y ~ 1 + (1 + x | gr(g, prec = Q)),
  gp = y ~ gp(x),
  hsgp = y ~ gp(x, k = 5))
row <- function(nm, which, v) {
  if (inherits(v, "error")) {
    cat(sprintf("%-8s %-10s %-9s %s
", nm, which, "ERROR", conditionMessage(v)))
    return(invisible())
  }
  cat(sprintf("%-8s %-10s %-9s %s\n", nm, which,
              if (is.nan(v)) "NaN" else if (identical(v, Inf)) "+Inf"
              else if (identical(v, -Inf)) "-Inf" else "finite",
              format(v, digits = 6)))
}
for (nm in names(forms)) {
  fr <- tryCatch(frm(bf(forms[[nm]]), family = gaussian(), data = dd,
                     data2 = list(A = A, Q = Q), dry_run = "frame"),
                 error = function(e) e)
  if (inherits(fr, "error")) {
    cat(nm, "FRAME ERROR:", conditionMessage(fr), "\n")
    next
  }
  blk <- fr$re_blocks[[1L]]
  reg <- frmtmb:::covstruct_registry[[blk$covstruct]]
  tpl <- fr$par_template
  b <- tpl$b + seq_along(tpl$b) / length(tpl$b)
  b0 <- b
  b0[seq(1L, length(b0), by = 2L)] <- 0
  th <- tpl$theta[blk$theta_idx]
  sdi <- reg$sd_idx(blk$dim)
  cat(sprintf("== %s covstruct %s dim %d theta %d sd_idx %s\n", nm,
              blk$covstruct, blk$dim, length(th),
              paste(sdi, collapse = ",")))
  row(nm, "estimates", tryCatch(reg$nll(b, th, blk), error = identity))
  for (lsd in c(-1137.64)) {
    for (k in sdi) {
      t2 <- th
      t2[k] <- lsd
      row(nm, paste0("sd", k), tryCatch(reg$nll(b, t2, blk), error = identity))
      row(nm, paste0("sd", k, "b0"), tryCatch(reg$nll(b0, t2, blk), error = identity))
    }
    if (length(sdi) > 1L) {
      t2 <- th
      t2[sdi] <- lsd
      row(nm, "all", tryCatch(reg$nll(b, t2, blk), error = identity))
      row(nm, "all_b0", tryCatch(reg$nll(b0, t2, blk), error = identity))
    }
  }
}

# The same on RTMB's tape, which is what a sampler evaluates: a numeric
# solve() refuses a near-singular matrix that the tape's solve takes.
cat("\n==== on the tape\n")
for (nm in names(forms)) {
  fr <- tryCatch(frm(bf(forms[[nm]]), family = gaussian(), data = dd,
                     data2 = list(A = A, Q = Q), dry_run = "frame"),
                 error = function(e) e)
  if (inherits(fr, "error")) next
  blk <- fr$re_blocks[[1L]]
  reg <- frmtmb:::covstruct_registry[[blk$covstruct]]
  tpl <- fr$par_template
  b <- tpl$b + seq_along(tpl$b) / length(tpl$b)
  th <- tpl$theta[blk$theta_idx]
  sdi <- reg$sd_idx(blk$dim)
  obj <- tryCatch(RTMB::MakeADFun(function(p) {
    -reg$nll(p$b, p$theta, blk)
  }, list(b = b, theta = th), silent = TRUE), error = identity)
  if (inherits(obj, "error")) {
    cat(nm, "TAPE ERROR:", conditionMessage(obj), "\n")
    next
  }
  ld <- function(t2) tryCatch(-obj$fn(c(b, t2)), error = identity)
  row(nm, "t:estim", ld(th))
  for (k in sdi) {
    t2 <- th
    t2[k] <- -1137.64
    row(nm, paste0("t:sd", k), ld(t2))
  }
  if (length(sdi) > 1L) {
    t2 <- th
    t2[sdi] <- -1137.64
    row(nm, "t:all", ld(t2))
  }
}
