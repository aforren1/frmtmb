# Reviewer of lane optstart: the exgaussian rebig regression, by arm.
#   Rscript dev/optstart-rev-exg.R <lib> <label>
# Design: dev/optstart-exg-within.R (group sd 30 * 1.5, sigma = beta =
# 1.5, n 400 in 20 groups), seeds 1:10. Arms isolate which part of the
# start moves the fit into the failing region.
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat(a[2], "lib", dirname(find.package("frmtmb")), "\n")
rows <- list()
for (s in 1:10) {
  set.seed(s)
  n <- 400
  g <- factor(rep(seq_len(20), length.out = n))
  x <- rnorm(n)
  bg <- rnorm(20, 0, 30 * 1.5)
  mu <- 0.3 + 0.5 * 1.5 * x + bg[g]
  y <- frmtmb:::family_registry$exgaussian()$sim(
    list(mu = mu, sigma = rep(1.5, n), beta = rep(1.5, n)), list(), n)
  d <- data.frame(y = y, x = x, g = g)
  q <- qr(cbind(1, x))
  r <- qr.resid(q, y)
  b_ls <- unname(qr.coef(q, y))
  tpl <- unclass(par_template(bf(y ~ x + (1 | g)), data = d,
                              family = exgaussian()))
  bd <- function(s1) c(`sigma_(Intercept)` = log(s1),
                       `beta_(Intercept)` = log(s1))
  s_lane <- sqrt(mean(r^2) / 2)
  s_old <- sd(y) / 2
  arms <- list(
    default = NULL,
    old = list(beta = c(mean(y), 0), betad = bd(s_old)),
    mu_old_sig_lane = list(beta = c(mean(y), 0), betad = bd(s_lane)),
    mu_lane_sig_old = list(beta = b_ls, betad = bd(s_old)),
    truth = list(beta = c(0.3, 0.75), betad = bd(1.5), theta = log(45)),
    truth_b = list(beta = c(0.3, 0.75), betad = bd(1.5),
                   theta = log(45), b = bg))
  for (nm in names(arms)) {
    msg <- NA_character_
    f <- tryCatch(suppressWarnings(frm(y ~ x + (1 | g), data = d,
      family = exgaussian(), start = arms[[nm]])),
      error = function(e) {
        msg <<- substr(conditionMessage(e), 1, 70)
        NULL
      })
    rows[[length(rows) + 1L]] <- data.frame(
      seed = s, arm = nm,
      ll = if (is.null(f)) NA else as.numeric(logLik(f)),
      code = if (is.null(f)) NA else f$opt$convergence,
      err = msg, s_lane = s_lane, s_old = s_old,
      tpl_sigma = exp(tpl$betad[[1]]), tpl_b1 = tpl$beta[[1]])
  }
}
X <- do.call(rbind, rows)
best <- ave(X$ll, X$seed, FUN = function(v) max(c(v, -Inf), na.rm = TRUE))
X$gap <- X$ll - best
print(X[, c("seed", "arm", "ll", "code", "gap", "err")], digits = 8)
print(aggregate(cbind(error = is.na(ll), short = !is.na(ll) & gap < -1e-4)
                ~ arm, data = X, FUN = sum, na.action = na.pass))
print(unique(X[, c("seed", "s_lane", "s_old", "tpl_sigma", "tpl_b1")]))
