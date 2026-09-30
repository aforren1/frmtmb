# Punch 2, n3: the kappa check on the re-check's cases
# (dev/fams2-rev2-guards.R): the nine Beta() data sets it stayed quiet
# on, kappa ~ x at seed 31, the A1 and 3b data (kappa placed), and the
# timing of the fit with and without the check's sdreport().
.libPaths(c("C:/Users/adf44/source/r/wt-fams2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
capw <- function(expr) {
  w <- character()
  v <- withCallingHandlers(expr, warning = function(cw) {
    w <<- c(w, conditionMessage(cw)); invokeRestart("muffleWarning")
  })
  list(v = v, w = w)
}
kw <- function(w) {
  k <- w[grepl("^xbeta:", w)]
  if (length(k)) sub(".*does not place it: ", "", sub("[.] Its.*", "", k)) else "-"
}
cat("== Beta() data, the nine the round-2 check was quiet on ==\n")
for (cfg in list(c(400, 5, 9), c(2000, 5, 2), c(2000, 5, 7), c(400, 50, 4),
                 c(400, 50, 5), c(400, 50, 7), c(400, 50, 8), c(400, 50, 9),
                 c(400, 50, 10))) {
  set.seed(cfg[3]); n <- cfg[1]; phi <- cfg[2]
  x <- rnorm(n); mu <- plogis(0.2 + 0.5 * x)
  y <- rbeta(n, mu * phi, (1 - mu) * phi)
  f <- capw(frm(y ~ x, family = xbeta(), data = data.frame(y = y, x = x)))
  se <- frmtmb:::xbeta_kappa_se(f$v, f$v$frame$linpreds[[
    frmtmb:::linpred_key("y", "kappa")]])
  cat(sprintf("  n %d phi %g seed %d: kappa %.3g, SE %.3g | %s\n", n, phi,
              cfg[3], max(frmtmb:::eval_dpars(f$v)[[1]]$kappa), se, kw(f$w)))
}
cat("== kappa ~ x, seed 31, n 400 ==\n")
set.seed(31)
x <- rnorm(400); mu <- plogis(0.2 + 0.5 * x)
d31 <- data.frame(x = x, y = rbeta(400, mu * 5, (1 - mu) * 5))
f <- capw(frm(bf(y ~ x, kappa ~ x), family = xbeta(), data = d31))
cat("  ", kw(f$w), "\n")
cat("== kappa placed, no row at 0 or 1 (3b) ==\n")
for (cfg in list(c(1, 200, 2000, 4101), c(0.3, 50, 1000, 4102),
                 c(0.5, 100, 1000, 4103))) {
  set.seed(cfg[4])
  z <- rbeta(cfg[3], 0.5 * cfg[2], 0.5 * cfg[2])
  y <- (1 + 2 * cfg[1]) * z - cfg[1]
  f <- capw(frm(y ~ 1, family = xbeta(), data = data.frame(y = y)))
  se <- frmtmb:::xbeta_kappa_se(f$v, f$v$frame$linpreds[[
    frmtmb:::linpred_key("y", "kappa")]])
  cat(sprintf("  kappa %g phi %g n %g: kappa-hat %.3g, SE %.3g | %s | other: %s\n",
              cfg[1], cfg[2], cfg[3], max(frmtmb:::eval_dpars(f$v)[[1]]$kappa),
              se, kw(f$w), paste(substr(f$w[!grepl("^xbeta:", f$w)], 1, 50),
                                 collapse = "; ")))
}
cat("== cost: seed 31 (kappa ~ 1), n 400; minimum of 5 fits, each on data ==\n")
cat("== moved by a relative 1e-9 so that no fit is served from a cache ==\n")
tm <- function(d) {
  min(vapply(1:5, function(i) {
    di <- d; di$y <- di$y * (1 - 1e-9 * i)
    system.time(suppressWarnings(frm(y ~ x, family = xbeta(),
                                     data = di)))[["elapsed"]]
  }, 0))
}
t_warn <- tm(d31)
d0 <- d31; d0$y[1] <- 0
t_end <- tm(d0)
cat(sprintf("  no end (check runs sdreport): %.3f s; one row at 0 (no sdreport): %.3f s\n",
            t_warn, t_end))
