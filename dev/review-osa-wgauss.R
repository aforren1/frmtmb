# Review lane osa, short check: where the 0.0018 between a weight-2
# gaussian (1 | g) fit and its duplicated rows comes from. Same data as
# dev/review-osa-weights2.R item 1. Exact: quadrature over the group
# effect, rows of weight w entering the conditioning as p^w, the residual
# of the last copy as int F p^(w-1) / int (F + S) p^(w-1).
.libPaths(c("C:/Users/adf44/source/r/wt-osa-lib",
            "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib: frmtmb", as.character(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "\n")
cat("BEGIN review-osa-wgauss\n")
G <- 20L
K <- 10L
set.seed(801)
d <- data.frame(g = factor(rep(seq_len(G), each = K)), t = rep(1:K, G),
                x = rnorm(G * K))
u <- rnorm(G, 0, 1)[d$g]
d$yg <- 1 + 0.5 * d$x + u + rnorm(G * K)
d$w <- ifelse(d$t <= 3, 2, 1)
dd <- d[rep(seq_len(nrow(d)), d$w), ]
last <- !duplicated(rep(seq_len(nrow(d)), d$w), fromLast = TRUE)
fw <- frm(yg | weights(w) ~ x + (1 | g), family = gaussian(), data = d)
fd <- frm(yg ~ x + (1 | g), family = gaussian(), data = dd)
r <- function(f, blk) {
  op <- options(frmtmb.osa_block = blk)
  on.exit(options(op))
  residuals(f, type = "osa")[, "Estimate"]
}
rw_b <- r(fw, 0L)
rw_1 <- r(fw, Inf)
rd_b <- r(fd, 0L)[last]
rd_1 <- r(fd, Inf)[last]
mu <- fw$estimates$beta[[1]] + fw$estimates$beta[[2]] * d$x
sg <- exp(fw$estimates$betad[[1]] %||% NA)
dp <- lapply(frmtmb:::eval_dpars(fw)[[1]], as.numeric)
sg <- dp$sigma[1]
sdu <- exp(fw$estimates$theta[[1]])
s <- seq(-10, 10, length.out = 4001) * sdu
ez <- numeric(nrow(d))
for (gi in split(seq_len(nrow(d)), d$g)) {
  lw <- dnorm(s, 0, sdu, log = TRUE)
  for (i in gi) {
    lp <- dnorm(d$yg[i], mu[i] + s, sg, log = TRUE)
    lF <- pnorm(d$yg[i], mu[i] + s, sg, log.p = TRUE)
    lS <- pnorm(d$yg[i], mu[i] + s, sg, lower.tail = FALSE, log.p = TRUE)
    a <- lw + (d$w[i] - 1) * lp
    m <- max(a)
    Fi <- sum(exp(a - m + lF))
    Si <- sum(exp(a - m + lS))
    ez[i] <- qnorm(Fi / (Fi + Si))
    lw <- lw + d$w[i] * lp
  }
}
w2 <- d$w == 2
f <- function(a, b) sprintf("max %.3g (weight-2 rows %.3g)", max(abs(a - b)),
                            max(abs(a - b)[w2]))
cat("weighted blocked vs weighted one call:  ", f(rw_b, rw_1), "\n")
cat("duplicated blocked vs one call:         ", f(rd_b, rd_1), "\n")
cat("weighted blocked vs duplicated blocked: ", f(rw_b, rd_b), "\n")
cat("weighted blocked vs exact:              ", f(rw_b, ez), "\n")
cat("duplicated blocked vs exact:            ", f(rd_b, ez), "\n")
cat("END review-osa-wgauss\n")
