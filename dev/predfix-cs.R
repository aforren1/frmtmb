# cs() in predict(): the simulated category proportions against the
# exact category probabilities fitted() gives, in sample and at
# newdata, for sratio, acat and cumulative. z = (proportion - p) /
# sqrt(p (1 - p) / ndraws) per cell; |z| above 5 is not Monte Carlo.
#   PREDFIX_ARM=base|lane Rscript dev/predfix-cs.R > dev/predfix-log/cs-<arm>.txt
source("C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-prelude.R")
ND <- 20000L
set.seed(11)
n <- 400
x <- rnorm(n)
# category-specific effect: the first threshold moves with x, the second
# the other way, so the effect is not proportional
eta1 <- -0.3 + 1.5 * x
eta2 <- 0.8 - 1.2 * x
u <- runif(n)
p1 <- plogis(eta1)
p2 <- (1 - p1) * plogis(eta2)
yo <- ifelse(u < p1, 1L, ifelse(u < p1 + p2, 2L, 3L))
d <- data.frame(x = x, yo = yo, yf = factor(yo, ordered = TRUE))
maxz <- function(prop, p) {
  se <- sqrt(pmax(p * (1 - p), 1e-12) / ND)
  max(abs(prop - p) / se)
}
fams <- list(sratio = function() sratio(), acat = function() acat(),
             cumulative = function() cumulative())
nd <- data.frame(x = c(3, 3, -3, 0))
for (nm in names(fams)) {
  f <- tryCatch(frm(bf(yo ~ cs(x)), family = fams[[nm]](), data = d),
                error = function(e) e)
  if (inherits(f, "error")) {
    cat(nm, ": frm() refuses cs(x):", conditionMessage(f), "\n")
    next
  }
  ex_in <- fitted(f)[, "Estimate", ]
  set.seed(1)
  pr_in <- predict(f, ndraws = ND)
  i <- which.max(d$x)
  cat(sprintf("%-10s in sample, max |z| over all rows %.1f; row x = %.3f: predict (%s) fitted (%s)\n",
              nm, maxz(pr_in, ex_in), d$x[i],
              paste(sprintf("%.3f", pr_in[i, ]), collapse = ", "),
              paste(sprintf("%.3f", ex_in[i, ]), collapse = ", ")))
  ex_nd <- fitted(f, newdata = nd)[, "Estimate", ]
  set.seed(1)
  pr_nd <- predict(f, newdata = nd, ndraws = ND)
  cat(sprintf("%-10s newdata, max |z| %.1f\n", nm, maxz(pr_nd, ex_nd)))
  for (r in seq_len(nrow(nd))) {
    cat(sprintf("    x = %4.1f: predict (%s) fitted (%s)\n", nd$x[r],
                paste(sprintf("%.3f", pr_nd[r, ]), collapse = ", "),
                paste(sprintf("%.3f", ex_nd[r, ]), collapse = ", ")))
  }
}
