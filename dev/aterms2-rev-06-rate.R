# Reviewer, claim 4: rate() post-fit semantics and refusals. Seed 606.
# Log: dev/aterms2-rev-log-06-rate.txt
.libPaths(c("C:/Users/adf44/source/r/wt-aterms2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
q <- function(expr) suppressWarnings(suppressMessages(expr))
show <- function(label, expr) {
  r <- tryCatch(q(expr), error = function(e) paste("ERROR:",
                                                   conditionMessage(e)))
  cat(sprintf("%-40s ", label))
  if (is.character(r) && length(r) == 1L) cat(substr(r, 1, 200), "\n")
  else cat("ok", class(r)[1], "\n")
  invisible(r)
}
set.seed(606)
n <- 200
d <- data.frame(x = rnorm(n), time = runif(n, 0.5, 4))
d$y <- rpois(n, exp(0.2 + 0.5 * d$x) * d$time)
d$yn <- rnbinom(n, mu = exp(0.2 + 0.5 * d$x) * d$time, size = 2 * d$time)
fr <- q(frm(y | rate(time) ~ x, data = d, family = poisson()))
fo <- q(frm(y ~ x + offset(log(time)), data = d, family = poisson()))

cat("== rate() against offset(log(time)) on poisson\n")
cat("  fixef identical:", identical(fixef(fr), fixef(fo)), "\n")
cat("  fitted identical:", identical(fitted(fr), fitted(fo)), " max rel",
    max(abs(fitted(fr)[, 1] - fitted(fo)[, 1]) / fitted(fo)[, 1]), "\n")
nd <- d[1:8, ]
nd$time <- c(0.1, 1, 2, 5, 10, 50, 100, 1000)
cat("  fitted(newdata) max rel:",
    max(abs(fitted(fr, newdata = nd)[, 1] - fitted(fo, newdata = nd)[, 1]) /
          fitted(fo, newdata = nd)[, 1]), "\n")
set.seed(11)
pr <- predict(fr, newdata = nd, ndraws = 4000)
set.seed(11)
po <- predict(fo, newdata = nd, ndraws = 4000)
cat("  predict(newdata) Estimate rate/offset:",
    format(pr[, 1] / po[, 1], digits = 4), "\n")
s1 <- simulate(fr, nsim = 2, seed = 13)
s2 <- simulate(fo, nsim = 2, seed = 13)
cat("  simulate identical:", identical(s1, s2), " rows differing",
    sum(s1 != s2), "of", length(unlist(s1)), "\n")
cat("  pearson residuals max rel:",
    max(abs(residuals(fr, type = "pearson")[, 1] -
              residuals(fo, type = "pearson")[, 1])) /
      max(abs(residuals(fo, type = "pearson")[, 1])), "\n")
ce_r <- q(conditional_effects(fr))[[1]]
# the offset model's conditional_effects() fails on both arms
# (dev/aterms2-rev-log-06b-ceoffset.txt), so the rate model is checked
# against its own mu at the held exposure
mu_ce <- fitted(fr, newdata = ce_r, dpar = "mu")[, 1]
cat("  conditional_effects rate: estimate / (mu * held time), range:",
    format(range(ce_r$estimate__ / (mu_ce * ce_r$time)), digits = 8), "\n")
cat("  conditional_effects rate: time held at", unique(ce_r$time),
    " mean(time)", mean(d$time), "\n")

cat("\n== negbinomial and geometric: draws follow newdata's denom\n")
fn <- q(frm(yn | rate(time) ~ x, data = d, family = negbinomial()))
fg <- q(frm(yn | rate(time) ~ x, data = d, family = geometric()))
nd2 <- d[c(1, 1, 1), ]
nd2$time <- c(1, 10, 100)
for (f in list(fn, fg)) {
  ep <- fitted(f, newdata = nd2)[, 1]
  mu <- fitted(f, newdata = nd2, dpar = "mu")[, 1]
  set.seed(5)
  pp <- predict(f, newdata = nd2, ndraws = 20000, summary = FALSE)
  sh <- if (f$spec$responses[[1]]$family$family == "negbinomial") {
    exp(f$opt$par[names(f$opt$par) == "betad"][1])
  } else 1
  cat(sprintf("  %-12s epred/(mu*time) %s  draw mean/epred %s\n",
              f$spec$responses[[1]]$family$family,
              paste(format(ep / (mu * nd2$time), digits = 6), collapse = " "),
              paste(format(colMeans(pp) / ep, digits = 4), collapse = " ")))
  # brms's variance: mu d (1 + mu / shape) with shape * d
  vtheo <- ep * (1 + mu / sh)
  cat(sprintf("  %-12s draw var / (mu d (1 + mu/shape)) %s\n", "",
              paste(format(apply(pp, 2, var) / vtheo, digits = 4),
                    collapse = " ")))
  ss <- simulate(f, nsim = 400, seed = 3)
  cat(sprintf("  %-12s simulate mean / fitted over rows: %.4f\n", "",
              mean(rowMeans(as.matrix(ss))) / mean(fitted(f)[, 1])))
}

cat("\n== refusals\n")
d0 <- d
d0$time[3] <- 0
show("denom 0", frm(y | rate(time) ~ x, data = d0, family = poisson()))
d0$time[3] <- -1
show("denom negative", frm(y | rate(time) ~ x, data = d0, family = poisson()))
d0$time[3] <- NA
f_na <- show("denom NA (row dropped?)",
             frm(y | rate(time) ~ x, data = d0, family = poisson()))
if (inherits(f_na, "frmtmb_fit")) cat("    n_obs", f_na$frame$n_obs, "\n")
cat("  brms denom 0:", tryCatch({q(brms::standata(
  brms::bf(y | rate(time) ~ x), transform(d, time = replace(time, 3, 0)),
  family = poisson())); "ok"}, error = function(e) conditionMessage(e)),
  "\n")
ndz <- nd
ndz$time[2] <- 0
show("fitted(newdata denom 0)", fitted(fr, newdata = ndz))
r0 <- q(fitted(fr, newdata = ndz))
if (is.matrix(r0)) cat("    value at denom 0:", r0[2, 1], "\n")
ndz$time[2] <- -2
show("fitted(newdata denom -2)", fitted(fr, newdata = ndz))
r0 <- q(fitted(fr, newdata = ndz))
if (is.matrix(r0)) cat("    value at denom -2:", r0[2, 1], "\n")
show("predict(newdata denom -2)", predict(fr, newdata = ndz, ndraws = 5))

cat("  brms newdata denom -2 (standata newdata):",
    tryCatch({q(brms::standata(brms::bf(y | rate(time) ~ x), ndz,
                                 family = poisson())); "ok"},
             error = function(e) conditionMessage(e)), "\n")
show("zero_inflated_poisson + rate",
     frm(y | rate(time) ~ x, data = d, family = zero_inflated_poisson()))
show("hurdle_poisson + rate",
     frm(y | rate(time) ~ x, data = d, family = hurdle_poisson()))
show("gaussian + rate", frm(y | rate(time) ~ x, data = d))
show("binomial + rate", frm(y | rate(time) + trials(200) ~ x, data = d,
                            family = binomial()))
show("rate(time * 2) expression",
     frm(y | rate(time * 2) ~ x, data = d, family = poisson()))
cat("  brms rate(time * 2):",
    tryCatch({q(brms::standata(brms::bf(y | rate(time * 2) ~ x), d,
                                 family = poisson())); "ok"},
             error = function(e) conditionMessage(e)), "\n")
show("rate(2) constant", frm(y | rate(2) ~ x, data = d, family = poisson()))
cat("  brms rate(2):",
    tryCatch({s <- q(brms::standata(brms::bf(y | rate(2) ~ x), d,
                                    family = poisson())); length(s$denom)},
             error = function(e) conditionMessage(e)), "\n")
