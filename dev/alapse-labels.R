# Lane alapse: the refusal's words for the structures added in 5449.
#   Rscript dev/alapse-labels.R <base|lane>
lib_arg <- commandArgs(TRUE)[1]
source("dev/alapse-common.R")
set.seed(20261116)
d <- data.frame(g = factor(rep(1:20, each = 8)),
                tt = factor(rep(rep(1:4, each = 2), 20)))
d$x <- stats::rnorm(160)
d$tn <- num_factor(c(0, 0.3, 0.9, 1.5)[as.integer(d$tt)])
d$y <- 1 + 0.5 * d$x + stats::rnorm(160)
Q <- diag(20) * 2
for (i in 1:19) Q[i, i + 1] <- Q[i + 1, i] <- -0.9
dimnames(Q) <- list(levels(d$g), levels(d$g))
for (tm in c("hetar1(tt + 0 | g)", "ou(tn + 0 | g)", "toep(tt + 0 | g)",
             "homtoep(tt + 0 | g)", "homcs(tt + 0 | g)", "exp(tn + 0 | g)",
             "gau(tn + 0 | g)", "mat(tn + 0 | g)", "(1 | gr(g, prec = Q))")) {
  fo <- stats::as.formula(paste("y ~ x +", tm))
  r <- tryCatch(frm(bf(fo), family = asym_laplace(), data = d,
                    data2 = list(Q = Q)), error = function(e) e)
  cat(tm, "->", if (inherits(r, "error")) substr(conditionMessage(r), 1, 70)
      else "FIT", "\n")
}
