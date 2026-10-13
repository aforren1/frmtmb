# Lane optrap, item 1: what the optimizer does on the ri20 design with a
# rescaled response. Per seed and scale: the default fit's code, message,
# iterations and logLik gap to lme4, and the same with the group log sd
# started at log(sd(y)) (the start on the response's scale).
#   Rscript dev/optrap-diag1.R <lib> [seeds]
args <- commandArgs(TRUE)
.libPaths(c(if (nzchar(args[1])) args[1],
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
seeds <- if (length(args) > 1) eval(parse(text = args[2])) else 1:40
suppressMessages({library(frmtmb); library(lme4)})
cat("lib", find.package("frmtmb"), "\n")
q <- function(expr) suppressMessages(suppressWarnings(expr))
X <- list()
for (s in seeds) for (sc in c(1e-3, 1, 1e3)) {
  set.seed(s)
  d <- data.frame(g = factor(rep(1:20, each = 5)), x = rnorm(100))
  d$y <- (1 + 0.5 * d$x + rnorm(100)) * sc
  l4 <- q(lmer(y ~ x + (1 | g), data = d, REML = FALSE))
  ll4 <- as.numeric(logLik(l4))
  f0 <- q(frm(y ~ x + (1 | g), data = d))
  f1 <- q(frm(y ~ x + (1 | g), data = d,
              start = list(theta = log(sd(d$y)))))
  X[[length(X) + 1L]] <- data.frame(
    seed = s, scale = sc,
    code0 = f0$opt$convergence, it0 = f0$opt$iterations %||% NA,
    msg0 = substr(f0$opt$message, 1, 28),
    d0 = as.numeric(logLik(f0)) - ll4,
    code1 = f1$opt$convergence, d1 = as.numeric(logLik(f1)) - ll4,
    rsd0 = exp(f0$estimates$theta[1] - f0$estimates$betad[1]),
    rsd4 = attr(VarCorr(l4)$g, "stddev") / sigma(l4))
}
X <- do.call(rbind, X)
dir.create("dev/optrap-log", showWarnings = FALSE)
write.table(X, sprintf("dev/optrap-log/diag1-%s.tsv",
                       basename(dirname(find.package("frmtmb")))),
            sep = "\t", quote = FALSE, row.names = FALSE)
for (sc in unique(X$scale)) {
  x <- X[X$scale == sc, ]
  cat(sprintf("scale %-6g default: short %d (min %.4g) | theta=log sd(y): short %d (min %.4g)\n",
              sc, sum(x$d0 < -1e-4), min(x$d0), sum(x$d1 < -1e-4), min(x$d1)))
}
print(X[X$d0 < -1e-4 | X$d1 < -1e-4, ], row.names = FALSE, digits = 4)
