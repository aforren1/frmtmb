# Punch round 1, m4: does the gap floor of ord_log_interior() change a
# plain fit's numbers, on the tape as well as in doubles? 100000 random
# threshold pairs (a > b, on the log-odds scale of the link): the value
# and the gradient taped through RTMB, floor on against floor off,
# compared as doubles. If every pair is identical, keeping the floor
# out of plain fits is not needed for bitwise identity, and the
# comment that says so is wrong.
.libPaths(c("C:/Users/adf44/source/r/wt-ordmix-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
oli <- frmtmb:::ord_log_interior
set.seed(20261005)
n <- 100000
b <- runif(n, -40, 40)
# gaps from 1e-12 to 50, log-uniform, as fits reach both ends
a <- b + exp(runif(n, log(1e-12), log(50)))
x <- c(a, b)
tape <- function(on) {
  RTMB::MakeTape(function(p) oli(p[seq_len(n)], p[n + seq_len(n)], on), x)
}
ton <- tape(TRUE)
toff <- tape(FALSE)
v_on <- ton(x)
v_off <- toff(x)
# the function is separable, so the gradient of the sum is each pair's
# own two partial derivatives
gsum <- function(on) {
  RTMB::MakeTape(function(p) sum(oli(p[seq_len(n)], p[n + seq_len(n)], on)),
                 x)
}
g_on <- gsum(TRUE)$jacobian(x)
g_off <- gsum(FALSE)$jacobian(x)
cat(sprintf("values identical on %d of %d pairs (finite on %d)\n",
            sum(v_on == v_off | (is.na(v_on) & is.na(v_off))), n,
            sum(is.finite(v_off))))
cat(sprintf("gradients identical on %d of %d entries\n",
            sum(g_on == g_off | (is.na(g_on) & is.na(g_off))), 2 * n))
d <- which(v_on != v_off)
if (length(d)) {
  cat("first differing pairs (a, b, on, off):\n")
  print(head(cbind(a[d], b[d], v_on[d], v_off[d])))
}
dg <- which(g_on != g_off)
if (length(dg)) {
  cat(sprintf("max relative gradient difference %.3g\n",
              max(abs(g_on[dg] - g_off[dg]) / abs(g_off[dg]))))
}
if (length(dg)) {
  gap <- (a - b)[(dg - 1L) %% n + 1L]
  rel <- abs(g_on[dg] - g_off[dg]) / pmax(abs(g_off[dg]), abs(g_on[dg]))
  cat(sprintf("differing entries: gap range %.3g to %.3g; b range %.3g to %.3g\n",
              min(gap), max(gap), min(b[(dg - 1L) %% n + 1L]),
              max(b[(dg - 1L) %% n + 1L])))
  cat("relative difference quantiles (of the larger):\n")
  print(quantile(rel, c(0, 0.5, 0.9, 0.99, 1)))
  cat("off-gradient zero where on is not:", sum(g_off[dg] == 0), "\n")
  k <- dg[order(-rel)][1:5]
  i <- (k - 1L) %% n + 1L
  print(data.frame(a = a[i], b = b[i], gap = a[i] - b[i], wrt = ifelse(k > n, "b", "a"),
                   on = g_on[k], off = g_off[k]))
}
