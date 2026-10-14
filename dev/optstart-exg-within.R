# Lane optstart: the battery's exgaussian rebig design (group sd 30,
# sigma = beta = 1.5; dev/optstart-battery.R) from starts given through
# `start =`: lane (mu at least squares, sigma = beta = rms / sqrt(2),
# the group sd at the unit, the lane's default given explicitly),
# within (sigma and beta from the residuals less their group means),
# anova (as within, and the group sd from the group means), and old
# (0.70.0's). Per arm: fits that end in an error, and the log likelihood
# of the rest against the best arm.
#   Rscript dev/optstart-exg-within.R <lib> [seeds]
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
seeds <- if (length(a) > 1) eval(parse(text = a[2])) else 1:10
cat("lib", dirname(find.package("frmtmb")), "\n")
out <- list()
for (s in seeds) {
  set.seed(s)
  n <- 400
  g <- factor(rep(seq_len(20), length.out = n))
  x <- rnorm(n)
  bg <- rnorm(20, 0, 30 * 1.5)
  mu <- 0.3 + 0.5 * 1.5 * x + bg[g]
  y <- frmtmb:::family_registry$exgaussian()$sim(
    list(mu = mu, sigma = rep(1.5, n), beta = rep(1.5, n)), list(), n)
  d <- data.frame(y = y, x = x, g = g)
  X <- cbind(1, x)
  q <- qr(X)
  r <- qr.resid(q, y)
  rw <- r - ave(r, g)
  ng <- mean(table(g))
  sb <- sqrt(max(var(tapply(r, g, mean)) - mean(rw^2) / ng, 1e-8))
  b <- unname(qr.coef(q, y))
  th <- unclass(par_template(bf(y ~ x + (1 | g)), data = d,
                             family = exgaussian()))$theta
  bd <- function(s1) c(`sigma_(Intercept)` = log(s1), `beta_(Intercept)` = log(s1))
  arms <- list(
    lane = list(beta = b, betad = bd(sqrt(mean(r^2) / 2)), theta = th),
    within = list(beta = b, betad = bd(sqrt(mean(rw^2) / 2)), theta = th),
    anova = list(beta = b, betad = bd(sqrt(mean(rw^2) / 2)), theta = log(sb)),
    old = list(beta = c(mean(y), 0), betad = bd(sd(y) / 2), theta = th))
  for (nm in names(arms)) {
    f <- tryCatch(suppressWarnings(frm(y ~ x + (1 | g), data = d,
      family = exgaussian(), start = arms[[nm]])), error = function(e) NULL)
    out[[length(out) + 1L]] <- data.frame(seed = s, arm = nm,
      ll = if (is.null(f)) NA else as.numeric(logLik(f)),
      code = if (is.null(f)) NA else f$opt$convergence)
  }
}
X <- do.call(rbind, out)
best <- ave(X$ll, X$seed, FUN = function(v) max(c(v, -Inf), na.rm = TRUE))
X$gap <- X$ll - best
print(aggregate(cbind(error = is.na(ll), short = !is.na(ll) & gap < -1e-4,
                      code = !is.na(code) & code != 0) ~ arm, data = X,
                FUN = sum, na.action = na.pass))
