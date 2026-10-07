# Reviewer of lane surface, claim 9: the blocks the floor does not
# cover, at the log sd the gp chain reached (-1137.64), with b != 0.
# +Inf here is the same defect the lane fixed for six structures.
#   Rscript dev/surface-rev-otherblocks.R lane|base
source("dev/surface-rev-env.R")
arm <- rev_env(commandArgs(TRUE)[1])
suppressPackageStartupMessages(library(frmtmb))
set.seed(2)
n_g <- 30
times <- c(0, 1, 2, 3, 4, 5)
d <- data.frame(y = rnorm(n_g * 6), g = factor(rep(seq_len(n_g), each = 6)),
                tim = factor(rep(times, n_g)), x = rnorm(n_g * 6))
d$y <- d$y + rnorm(n_g, 0, 0.8)[d$g]
fits <- list(
  us1 = frm(bf(y ~ 1 + (1 | g)), family = gaussian(), data = d),
  us2 = frm(bf(y ~ 1 + (1 + x | g)), family = gaussian(), data = d),
  ar1 = suppressWarnings(frm(bf(y ~ 1 + ar1(tim + 0 | g)),
                             family = gaussian(), data = d)),
  diag = suppressWarnings(frm(bf(y ~ 1 + diag(tim + 0 | g)),
                              family = gaussian(), data = d)),
  cs = suppressWarnings(frm(bf(y ~ 1 + cs(tim + 0 | g)),
                            family = gaussian(), data = d)),
  homcs = suppressWarnings(frm(bf(y ~ 1 + homcs(tim + 0 | g)),
                               family = gaussian(), data = d)))
for (nm in names(fits)) {
  f <- fits[[nm]]
  blk <- f$frame$re_blocks[[1L]]
  reg <- frmtmb:::covstruct_registry[[blk$covstruct]]
  b <- f$estimates$b + 0.1
  th <- f$estimates$theta
  th2 <- th; th2[1] <- -1137.64
  v0 <- tryCatch(reg$nll(b, th, blk), error = function(e) NA)
  v1 <- tryCatch(reg$nll(b, th2, blk), error = function(e) conditionMessage(e))
  cat(sprintf("%-6s covstruct %-7s theta %d | log density at estimates %s | at log sd[1] -1137.64: %s\n",
              nm, blk$covstruct, length(th), format(v0, digits = 8),
              format(v1, digits = 8)))
}
