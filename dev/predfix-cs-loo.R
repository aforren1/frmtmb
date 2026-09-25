# Does frmtmb.sample's pointwise log-likelihood (loo.R:124, the nested
# with_cs_offsets() call) carry cs()? Its row sum at the ML estimates
# must equal logLik(fit).
#   PREDFIX_ARM=base|lane Rscript dev/predfix-cs-loo.R
source("C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-prelude.R")
suppressPackageStartupMessages(library(frmtmb.sample))
set.seed(11)
n <- 400
x <- rnorm(n)
p1 <- plogis(-0.3 + 1.5 * x)
p2 <- (1 - p1) * plogis(0.8 - 1.2 * x)
u <- runif(n)
d <- data.frame(x = x, yo = ifelse(u < p1, 1L, ifelse(u < p1 + p2, 2L, 3L)))
for (fam in list(sratio(), acat(), cratio())) {
  f <- frm(bf(yo ~ cs(x)), family = fam, data = d)
  ll <- frmtmb.sample:::draws_row_loglik(f, "yo")
  cat(sprintf("%-7s sum pointwise %.9f  logLik %.9f  difference %.3e\n",
              fam$family, sum(ll), as.numeric(logLik(f)),
              sum(ll) - as.numeric(logLik(f))))
}
