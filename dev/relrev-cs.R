# Reviewer, item 2: cs() on cumulative().
#  (a) on <lib>, every spelling of cs() on cumulative(): refused or fit;
#  (b) brms's warning count with two cs() predictors and cs(x + z);
#  (c) frmtmb's count on the same shapes;
#  (d) crossing rows: fitted() NaN and simulate() NA stay on the crossing
#      rows only; simulate() return types; predict().
# Seed 2026100601.
#   Rscript dev/relrev-cs.R <r5|r6> > dev/relrev-log/cs-<arm>.txt
arm <- commandArgs(TRUE)[1]
lib <- if (arm == "r6") "C:/Users/adf44/source/r/rellib-r6" else
  "C:/Users/adf44/source/r/rellib-r5"
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(brms); library(frmtmb)})
cat("arm", arm, find.package("frmtmb"), format(packageVersion("frmtmb")), "\n")
set.seed(2026100601)
n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n), w = rnorm(n),
                g = factor(sample(letters[1:4], n, TRUE)))
u <- rlogis(n) + 0.5 * d$x + 0.4 * d$z
d$y <- 1L + (u > -1 + 0.15 * d$x) + (u > 0.3) + (u > 1.5 - 0.15 * d$x)
d$y2 <- 1L + (rlogis(n) > 0)  + (rlogis(n) > 0.5)
d$yh <- ifelse(runif(n) < 0.2, 0L, d$y)
count <- function(expr) {
  w <- character()
  r <- tryCatch(withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage")),
  error = function(e) structure(conditionMessage(e), class = "err"))
  list(r = r, n = sum(grepl("Category specific effects for this family",
                            w, fixed = TRUE)), w = w)
}
show <- function(tag, res) {
  if (inherits(res$r, "err")) {
    cat(sprintf("%-50s REFUSED [%d warn] %s\n", tag, res$n,
                substr(gsub("\n", " ", res$r), 1, 120)))
  } else {
    cat(sprintf("%-50s FIT     [%d cs-warn, %d other warn]\n", tag, res$n,
                length(res$w) - res$n))
  }
}
cat("\n== (a) spellings of cs() on cumulative()\n")
sp <- list(
  "y ~ cs(x)" = quote(frm(y ~ cs(x), family = cumulative(), data = d)),
  "bf(y ~ cs(x))" = quote(frm(bf(y ~ cs(x)), family = cumulative(), data = d)),
  "family = 'cumulative'" = quote(frm(y ~ cs(x), family = "cumulative", data = d)),
  "bf(y ~ cs(x)) + cumulative()" = quote(frm(bf(y ~ cs(x)) + cumulative(), data = d)),
  "y ~ z + cs(x), probit" = quote(frm(y ~ z + cs(x), family = cumulative("probit"), data = d)),
  "y ~ cs(x) equidistant" = quote(frm(y ~ cs(x), family = cumulative(threshold = "equidistant"), data = d)),
  "y ~ cs(x) sum_to_zero" = quote(frm(y ~ cs(x), family = cumulative(threshold = "sum_to_zero"), data = d)),
  "y ~ cs(x) + (1 | g)" = quote(frm(y ~ cs(x) + (1 | g), family = cumulative(), data = d)),
  "y ~ (cs(1) | g)" = quote(frm(y ~ x + (cs(1) | g), family = cumulative(), data = d)),
  "y | thres(gr = g) ~ cs(x)" = quote(frm(y | thres(gr = g) ~ cs(x), family = cumulative(), data = d)),
  "mvbf cumulative cs" = quote(frm(mvbf(bf(y ~ cs(x), family = cumulative()), bf(y2 ~ x, family = cumulative())), data = d)),
  "hurdle_cumulative cs" = quote(frm(yh ~ cs(x), family = hurdle_cumulative(), data = d)),
  "mixture(cumulative, cumulative) cs" = quote(frm(y ~ cs(x), family = mixture(cumulative(), cumulative()), data = d)),
  "gaussian cs (control)" = quote(frm(z ~ cs(x), family = gaussian(), data = d)),
  "poisson (cs(1) | g) (control)" = quote(frm(y ~ x + (cs(1) | g), family = poisson(), data = d)),
  "sratio cs (control)" = quote(frm(y ~ cs(x), family = sratio(), data = d))
)
fits <- list()
for (nm in names(sp)) {
  res <- count(eval(sp[[nm]]))
  show(nm, res)
  if (!inherits(res$r, "err")) fits[[nm]] <- res$r
}
cat("\n== (b) brms 2.23.0 warning count (stancode)\n")
bsh <- list(
  "y ~ cs(x)" = list(brms::bf(y ~ cs(x)), brms::cumulative()),
  "y ~ cs(x) + cs(z)" = list(brms::bf(y ~ cs(x) + cs(z)), brms::cumulative()),
  "y ~ cs(x + z)" = list(brms::bf(y ~ cs(x + z)), brms::cumulative()),
  "y ~ cs(x) + cs(z) + cs(w)" = list(brms::bf(y ~ cs(x) + cs(z) + cs(w)), brms::cumulative()),
  "y ~ cs(x) + (cs(1) | g)" = list(brms::bf(y ~ cs(x) + (cs(1) | g)), brms::cumulative()),
  "mix(cum,cum) y ~ cs(x) + cs(z)" = list(brms::bf(y ~ cs(x) + cs(z)), brms::mixture(brms::cumulative(), brms::cumulative())),
  "mvbf two cumulative cs" = list(brms::mvbf(brms::bf(y ~ cs(x)), brms::bf(y2 ~ cs(z))), brms::cumulative()),
  "y ~ cs(x), gaussian" = list(brms::bf(z ~ cs(x)), gaussian()),
  "y ~ (cs(1) | g), poisson" = list(brms::bf(y ~ x + (cs(1) | g)), poisson())
)
for (nm in names(bsh)) {
  res <- count(as.character(brms::stancode(bsh[[nm]][[1]], data = d,
                                           family = bsh[[nm]][[2]])))
  show(paste("brms", nm), res)
}
if (arm == "r6") {
  cat("\n== (c) frmtmb 0.68.0 warning count, same shapes\n")
  fsh <- list(
    "y ~ cs(x)" = quote(frm(y ~ cs(x), family = cumulative(), data = d)),
    "y ~ cs(x) + cs(z)" = quote(frm(y ~ cs(x) + cs(z), family = cumulative(), data = d)),
    "y ~ cs(x + z)" = quote(frm(y ~ cs(x + z), family = cumulative(), data = d)),
    "y ~ cs(x) + cs(z) + cs(w)" = quote(frm(y ~ cs(x) + cs(z) + cs(w), family = cumulative(), data = d)),
    "mix(cum,cum) y ~ cs(x) + cs(z)" = quote(frm(y ~ cs(x) + cs(z), family = mixture(cumulative(), cumulative()), data = d)),
    "mvbf two cumulative cs" = quote(frm(mvbf(bf(y ~ cs(x)), bf(y2 ~ cs(z))), family = cumulative(), data = d))
  )
  for (nm in names(fsh)) show(paste("frmtmb", nm), count(eval(fsh[[nm]])))

  cat("\n== (d) crossing rows\n")
  fit <- fits[["y ~ cs(x)"]]
  est <- fit$estimates
  lp <- Filter(function(l) identical(l$dpar, "mu"), fit$frame$linpreds)[[1L]]
  b <- est[[lp$cs[[1L]]$par]]
  tau <- frmtmb:::ord_threshold_values(family(fit), est$tau_raw)
  cat("tau", format(tau, digits = 5), " bcs", format(b, digits = 5), "\n")
  # thresholds k and k+1 cross where tau_k - b_k x = tau_{k+1} - b_{k+1} x
  xc <- c((tau[2] - tau[1]) / (b[2] - b[1]), (tau[3] - tau[2]) / (b[3] - b[2]))
  cat("crossing x:", format(xc, digits = 5), "\n")
  xs <- c(-2, 0, 2, 3 * xc[1], 3 * xc[2], 1, -1)
  nd <- data.frame(x = xs)
  thr <- outer(rep(1, length(xs)), tau) - outer(xs, b)
  crossing <- apply(thr, 1, function(t) any(diff(t) < 0))
  cat("rows crossing:", which(crossing), "\n")
  P <- fitted(fit, newdata = nd)
  E <- P[, "Estimate", ]
  cat("fitted Estimate NaN by row:", rowSums(is.nan(E)), "\n")
  cat("fitted any non-finite in non-crossing rows:",
      any(!is.finite(E[!crossing, ])), "\n")
  cat("fitted row sums non-crossing:", format(rowSums(E[!crossing, ]), digits = 15), "\n")
  cat("fitted other columns NaN by row (Est.Error, Q2.5, Q97.5):\n")
  for (cl in setdiff(dimnames(P)[[2]], "Estimate")) {
    cat("  ", cl, rowSums(is.nan(P[, cl, ])), "\n")
  }
  pr <- tryCatch(predict(fit, newdata = nd), error = function(e) conditionMessage(e))
  cat("predict class:", class(pr), "\n"); print(pr)
  s <- simulate(fit, nsim = 200, seed = 3, newdata = nd)
  cat("simulate class:", class(s), " typeof:", typeof(as.matrix(s)), " dim:",
      dim(as.matrix(s)), "\n")
  sm <- as.matrix(s)
  cat("simulate NA count by row:", rowSums(is.na(sm)), "\n")
  s0 <- simulate(fit, nsim = 2, seed = 3)
  cat("in-sample simulate class:", class(s0), typeof(s0[[1]]),
      " NA:", sum(is.na(as.matrix(s0))), "\n")
  # the same fit with no crossing rows in sample (the fitted thresholds):
  thr_s <- outer(rep(1, n), tau) - outer(d$x, b)
  cat("in-sample crossing rows:", sum(apply(thr_s, 1, function(t) any(diff(t) < 0))), "\n")
  # a plain cumulative() fit's simulate type (no cs), for comparison
  f0 <- frm(y ~ x, family = cumulative(), data = d)
  s1 <- simulate(f0, nsim = 2, seed = 3)
  cat("plain cumulative simulate class:", class(s1), typeof(s1[[1]]), "\n")
  saveRDS(list(fit = fit, nd = nd, crossing = crossing),
          "C:/Users/adf44/source/r/frmtmb-wt-release/dev/relrev-log/cs-fit.rds")
}
