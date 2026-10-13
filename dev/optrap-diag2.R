# Lane optrap, item 1: does nlminb's step scale on the mu coefficients
# explain the short fits? A custom optimizer passes scale = 1/sd(y) for
# beta (arm "unit"), with and without the group log sd started at
# log(sd(y)) (arm "unit+theta").
#   Rscript dev/optrap-diag2.R <lib> [seeds]
args <- commandArgs(TRUE)
.libPaths(c(if (nzchar(args[1])) args[1],
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
seeds <- if (length(args) > 1) eval(parse(text = args[2])) else 1:40
suppressMessages({library(frmtmb); library(lme4)})
cat("lib", find.package("frmtmb"), "\n")
q <- function(expr) suppressMessages(suppressWarnings(expr))
mk_opt <- function(u) {
  function(par, fn, gr, lower, upper, control) {
    sc <- ifelse(names(par) == "beta", 1 / u, 1)
    stats::nlminb(par, fn, gr, control = control, scale = sc,
                  lower = lower, upper = upper)
  }
}
X <- list()
for (s in seeds) for (sc in c(1e-3, 1, 1e3)) {
  set.seed(s)
  d <- data.frame(g = factor(rep(1:20, each = 5)), x = rnorm(100))
  d$y <- (1 + 0.5 * d$x + rnorm(100)) * sc
  l4 <- q(lmer(y ~ x + (1 | g), data = d, REML = FALSE))
  ll4 <- as.numeric(logLik(l4))
  ctl <- frmtmb_control(optimizer = mk_opt(sd(d$y)))
  fa <- q(frm(y ~ x + (1 | g), data = d, control = ctl))
  fb <- q(frm(y ~ x + (1 | g), data = d, control = ctl,
              start = list(theta = log(sd(d$y)))))
  X[[length(X) + 1L]] <- data.frame(
    seed = s, scale = sc,
    codea = fa$opt$convergence, da = as.numeric(logLik(fa)) - ll4,
    codeb = fb$opt$convergence, db = as.numeric(logLik(fb)) - ll4)
}
X <- do.call(rbind, X)
for (sc in unique(X$scale)) {
  x <- X[X$scale == sc, ]
  cat(sprintf("scale %-6g unit: short %d (min %.4g) | unit+theta: short %d (min %.4g)\n",
              sc, sum(x$da < -1e-4), min(x$da), sum(x$db < -1e-4), min(x$db)))
}
print(X[X$da < -1e-4 | X$db < -1e-4, ], row.names = FALSE, digits = 4)
