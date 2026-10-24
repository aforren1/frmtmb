# Reviewer: the OLRE check on a car() block with one row per area, on base.
.libPaths(c("C:/Users/adf44/source/r/rellib-r11",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", as.character(packageVersion("frmtmb")), dirname(find.package("frmtmb")), "\n")
lattice_W <- function(r, c) {
  g <- expand.grid(r = seq_len(r), c = seq_len(c))
  n <- nrow(g)
  W <- 1 * (abs(outer(g$r, g$r, "-")) + abs(outer(g$c, g$c, "-")) == 1)
  dimnames(W) <- list(paste0("L", seq_len(n)), paste0("L", seq_len(n)))
  W
}
set.seed(42)
W <- lattice_W(4, 4); for (i in 1:96) rnorm(1)
invisible(rnorm(16)); invisible(rnorm(96))
W6 <- lattice_W(6, 6)
d1 <- data.frame(loc = factor(rownames(W6), levels = rownames(W6)), x = rnorm(36))
Q <- diag(rowSums(W6)) - 0.9 * W6
phi <- drop(crossprod(chol(solve(Q)), rnorm(36)))
d1$y <- 0.5 * d1$x + phi + rnorm(36, 0, 0.5)
withCallingHandlers(
  f <- frm(bf(y ~ x + car(W6, gr = loc, type = "escar")) + gaussian(),
           data = d1, data2 = list(W6 = W6)),
  warning = function(w) {cat("[warning]", substr(conditionMessage(w), 1, 200), "\n"); invokeRestart("muffleWarning")})
cat("logLik", format(as.numeric(logLik(f)), digits = 10), "\n")
