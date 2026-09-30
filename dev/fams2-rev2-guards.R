# Reviewer, punch round 1, items 3 and 4: the kappa warning (xbeta) and
# the disc-intercept warning (hurdle_cumulative), each where it should
# fire and where it should not.
.libPaths(c("C:/Users/adf44/source/r/wt-fams2-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
capw <- function(expr) {
  w <- character()
  v <- withCallingHandlers(
    tryCatch(expr, error = function(e) structure(conditionMessage(e), class = "err")),
    warning = function(cw) { w <<- c(w, conditionMessage(cw)); invokeRestart("muffleWarning") })
  list(v = v, w = w)
}
kw <- function(w) any(grepl("kappa ran to 0", w))
cat("== 3a. kappa really 0: Beta() data, 30 fits ==\n")
res <- NULL
for (cfg in list(c(100, 5), c(400, 5), c(2000, 5), c(400, 50))) {
  for (s in 1:10) {
    set.seed(s)
    n <- cfg[1]; phi <- cfg[2]
    x <- rnorm(n)
    mu <- plogis(0.2 + 0.5 * x)
    y <- rbeta(n, mu * phi, (1 - mu) * phi)
    f <- capw(frm(y ~ x, family = xbeta(), data = data.frame(y = y, x = x)))
    kmax <- if (inherits(f$v, "err")) NA else max(frmtmb:::eval_dpars(f$v)[[1]]$kappa)
    res <- rbind(res, data.frame(n = n, phi = phi, seed = s, err = inherits(f$v, "err"),
                                 kappa = kmax, kwarn = kw(f$w),
                                 other = paste(substr(setdiff(f$w, f$w[grepl("kappa ran", f$w)]), 1, 40),
                                               collapse = " | ")))
  }
}
print(res, digits = 3, row.names = FALSE)
cat(sprintf("fits %d, errors %d, kappa warning %d, kappa-hat >= 1e-6 without warning %d\n",
            nrow(res), sum(res$err), sum(res$kwarn),
            sum(!res$err & !res$kwarn & res$kappa >= 1e-6, na.rm = TRUE)))
set.seed(31)
x <- rnorm(400); mu <- plogis(0.2 + 0.5 * x)
d31 <- data.frame(x = x, y = rbeta(400, mu * 5, (1 - mu) * 5))
f <- capw(frm(y ~ x, family = xbeta(), data = d31))
cat("\nwording (seed 31, the worker's case):\n  ", f$w, sep = "\n  ")
f <- capw(frm(bf(y ~ x, kappa ~ x), family = xbeta(), data = d31))
cat("kappa ~ x on the same data: warnings", length(f$w), "; kappa range",
    if (!inherits(f$v, "err")) signif(range(frmtmb:::eval_dpars(f$v)[[1]]$kappa), 3) else f$v, "\n")
if (length(f$w)) cat("  ", substr(f$w, 1, 120), sep = "\n  ")

cat("\n== 3b. kappa not 0 though no row is at 0 or 1 ==\n")
for (cfg in list(c(1, 200, 2000, 4101), c(0.3, 50, 1000, 4102), c(0.5, 100, 1000, 4103))) {
  set.seed(cfg[4])
  z <- rbeta(cfg[3], 0.5 * cfg[2], 0.5 * cfg[2])
  y <- (1 + 2 * cfg[1]) * z - cfg[1]
  if (any(y <= 0 | y >= 1)) { cat("  (a row reached an end; skipped)\n"); next }
  f <- capw(frm(y ~ 1, family = xbeta(), data = data.frame(y = y)))
  fb <- frm(y ~ 1, family = Beta(), data = data.frame(y = y))
  if (inherits(f$v, "err")) {
    cat(sprintf("  kappa %g phi %g n %g: ERROR %s\n", cfg[1], cfg[2], cfg[3], substr(f$v, 1, 90)))
  } else {
    cat(sprintf("  kappa %g phi %g n %g: kappa-hat %.3g, logLik xbeta %.3f vs Beta %.3f, kappa warning %s, other warnings: %s\n",
                cfg[1], cfg[2], cfg[3], max(frmtmb:::eval_dpars(f$v)[[1]]$kappa),
                as.numeric(logLik(f$v)), as.numeric(logLik(fb)), kw(f$w),
                paste(substr(f$w[!grepl("kappa ran", f$w)], 1, 60), collapse = " | ")))
  }
}

cat("\n== 4. disc-intercept warning on hurdle_cumulative ==\n")
set.seed(4400)
n <- 1500
dD <- data.frame(x = rnorm(n), z = rnorm(n), g = gl(30, n / 30))
disc <- exp(0.4 * dD$z)
u <- rlogis(n) / disc + 0.8 * dD$x
yc <- 1L + (u > -1) + (u > 0.3) + (u > 1.5)
dD$y <- ifelse(runif(n) < plogis(-0.6 + 0.5 * dD$x), 0L, yc)
dD$y3 <- pmin(dD$y, 3L)
dw <- function(w) any(grepl("disc has an intercept", w))
chk <- function(lab, expr, se = FALSE) {
  r <- capw(expr)
  ee <- if (!inherits(r$v, "err")) {
    fe <- fixef(r$v)
    sprintf("max Est.Error %.3g", max(fe[, "Est.Error"]))
  } else paste("ERROR", substr(r$v, 1, 80))
  cat(sprintf("  %-55s disc warning %-5s | %s | other: %s\n", lab, dw(r$w), ee,
              paste(substr(r$w[!grepl("disc has an intercept", r$w)], 1, 70), collapse = " | ")))
  invisible(r)
}
r <- chk("disc ~ 1 + z", frm(bf(y ~ x, disc ~ 1 + z), data = dD, family = hurdle_cumulative()))
cat("  wording:", r$w[grepl("disc has", r$w)], "\n")
chk("disc ~ 1 (intercept only)", frm(bf(y ~ x, disc ~ 1), data = dD, family = hurdle_cumulative()))
chk("disc ~ z (implicit intercept)", frm(bf(y ~ x, disc ~ z), data = dD, family = hurdle_cumulative()))
chk("disc ~ 0 + z", frm(bf(y ~ x, disc ~ 0 + z), data = dD, family = hurdle_cumulative()))
chk("disc not modeled", frm(bf(y ~ x), data = dD, family = hurdle_cumulative()))
chk("disc ~ 1 + z, prior normal(0,1) Intercept dpar disc",
    frm(bf(y ~ x, disc ~ 1 + z), data = dD, family = hurdle_cumulative(),
        prior = set_prior("normal(0, 1)", class = "Intercept", dpar = "disc")))
chk("disc ~ 1 + z, prior on b dpar disc only",
    frm(bf(y ~ x, disc ~ 1 + z), data = dD, family = hurdle_cumulative(),
        prior = set_prior("normal(0, 1)", class = "b", dpar = "disc")))
chk("disc ~ 1 + z, prior on class Intercept (thresholds)",
    frm(bf(y ~ x, disc ~ 1 + z), data = dD, family = hurdle_cumulative(),
        prior = set_prior("normal(0, 3)", class = "Intercept")))
chk("disc ~ 1 + z, prior on b coef Intercept dpar disc",
    frm(bf(y ~ x, disc ~ 1 + z), data = dD, family = hurdle_cumulative(),
        prior = set_prior("normal(0, 1)", class = "b", coef = "Intercept", dpar = "disc")))
chk("disc ~ 0 + z + (1 | g)", frm(bf(y ~ x, disc ~ 0 + z + (1 | g)), data = dD,
                                  family = hurdle_cumulative()))
chk("thres(4) on 0..3, no prior (unplaced threshold check)",
    frm(bf(y3 | thres(4) ~ x), data = dD, family = hurdle_cumulative()))
chk("thres(4) on 0..3, disc ~ 1 + z, no prior",
    frm(bf(y3 | thres(4) ~ x, disc ~ 1 + z), data = dD, family = hurdle_cumulative()))
chk("thres(4) on 0..3 with an Intercept prior",
    frm(bf(y3 | thres(4) ~ x), data = dD, family = hurdle_cumulative(),
        prior = set_prior("normal(0, 3)", class = "Intercept")))
chk("cumulative() thres(4) on 1..3, no prior (sibling, unchanged?)",
    frm(bf(y3 | thres(4) ~ x), data = transform(dD, y3 = pmax(y3, 1L)), family = cumulative()))
