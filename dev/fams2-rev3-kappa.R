# Reviewer, punch round 2, item 3: the kappa check's new rules (a row
# below 1e-6; a log-scale SE above 10 or not finite) on fits where kappa
# IS placed, and the sdreport it reads on other fit paths.
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
sim <- function(seed, n, mu, phi, lk, x) {
  set.seed(seed)
  kap <- exp(lk)
  z <- rbeta(n, mu * phi, (1 - mu) * phi)
  y <- (1 + 2 * kap) * z - kap
  data.frame(x = x, y = y)
}
report <- function(lab, r, d) {
  if (inherits(r$v, "err")) { cat(sprintf("%-50s ERROR %s\n", lab, substr(r$v, 1, 90))); return() }
  kap <- frmtmb:::eval_dpars(r$v)[[1]]$kappa
  ci <- confint(r$v)
  kr <- grep("^kappa", rownames(ci))
  se <- sqrt(diag(vcov(r$v, full = TRUE)))[rownames(ci)[kr]]
  fb <- frm(y ~ 1, family = Beta(), data = d)
  cat(sprintf("%-50s ends %d | kappa rows %.2g..%.2g | kappa coefs %s (SE %s) | gain over Beta() %.1f | kappa warning: %s\n",
              lab, sum(d$y <= 0 | d$y >= 1), min(kap), max(kap),
              paste(signif(ci[kr, "est"], 3), collapse = ", "),
              paste(signif(se, 3), collapse = ", "),
              as.numeric(logLik(r$v) - logLik(fb)),
              if (any(grepl("does not place it", r$w))) "YES" else "no"))
  oth <- r$w[!grepl("does not place it", r$w)]
  if (length(oth)) cat("    other warnings:", paste(substr(oth, 1, 90), collapse = " | "), "\n")
  kw <- r$w[grepl("does not place it", r$w)]
  if (length(kw)) cat("    text:", kw, "\n")
}
cat("== 3a. kappa ~ x with a real effect and no row at 0 or 1 ==\n")
for (s in 1:3) {
  set.seed(100 + s); n <- 2000; x <- rnorm(n)
  d <- sim(200 + s, n, 0.5, 1000, -0.3 + 0.5 * x, x)
  report(sprintf("lk = -0.3 + 0.5 x, phi 1000, seed %d", 200 + s),
         capw(frm(bf(y ~ 1, kappa ~ x), family = xbeta(), data = d)), d)
}
cat("\n== 3b. a real effect that takes true kappa below 1e-6 at some rows ==\n")
for (s in 1:3) {
  set.seed(300 + s); n <- 2000; x <- runif(n, -4.7, 0)
  d <- sim(400 + s, n, 0.5, 200, -1 + 3 * x, x)
  cat(sprintf("  truth: kappa from %.2g to %.2g\n", exp(-1 + 3 * min(x)), exp(-1 + 3 * max(x))))
  report(sprintf("lk = -1 + 3 x, x in (-4.7, 0), phi 200, seed %d", 400 + s),
         capw(frm(bf(y ~ 1, kappa ~ x), family = xbeta(), data = d)), d)
}
cat("\n== 3c. a single well-placed kappa, larger n (SE should be small) ==\n")
set.seed(500); d <- sim(501, 5000, 0.5, 100, rep(log(0.5), 5000), rep(0, 5000))
report("kappa 0.5, phi 100, n 5000", capw(frm(y ~ 1, family = xbeta(), data = d)), d)

cat("\n== 3d. other fit paths: does reading the sdreport fail or leak? ==\n")
set.seed(600); n <- 800
g <- gl(40, n / 40)
d <- sim(601, n, 0.5, 5, rep(-25, n), rnorm(n))   # Beta() data: kappa should run to 0
d$g <- g
d$y <- plogis(qlogis(d$y) + rnorm(40, 0, 0.3)[g])
d2 <- sim(602, 2000, 0.5, 200, rep(0, 2000), rnorm(2000))  # kappa truly 1, placed weakly
paths <- list(
  "(1 | g), Beta() data" = quote(frm(y ~ 1 + (1 | g), family = xbeta(), data = d)),
  "(1 | g), profile = TRUE" = quote(frm(y ~ 1 + (1 | g), family = xbeta(), data = d,
                                         control = frmtmb_control(profile = TRUE))),
  "REML = TRUE, (1 | g)" = quote(frm(y ~ 1 + (1 | g), family = xbeta(), data = d, REML = TRUE)),
  "kappa link identity" = quote(frm(y ~ 1, family = xbeta(link_kappa = "identity"), data = d2)),
  "kappa link softplus" = quote(frm(y ~ 1, family = xbeta(link_kappa = "softplus"), data = d2)),
  "multivariate, two xbeta responses" = quote(frm(bf(y ~ 1) + xbeta() + bf(y2 ~ 1) + xbeta(),
                                                  data = transform(d2, y2 = d$y[seq_len(nrow(d2)) %% nrow(d) + 1])))
)
for (nm in names(paths)) {
  t0 <- proc.time()[["elapsed"]]
  r <- capw(eval(paths[[nm]]))
  el <- proc.time()[["elapsed"]] - t0
  cat(sprintf("%-38s %s | %.2f s | warnings: %s\n", nm,
              if (inherits(r$v, "err")) paste("ERROR", substr(r$v, 1, 80)) else "fits",
              el, if (length(r$w)) paste(substr(r$w, 1, 110), collapse = " || ") else "none"))
}
