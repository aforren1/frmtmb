# Quick checks of fix round 1 on the installed build: each item prints what the
# review measured, on the review's own designs.
.libPaths(c("C:/Users/adf44/source/r/wt-latres-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
options(warn = 1, width = 110)
catch <- function(expr) {
  w <- m <- character(0)
  r <- withCallingHandlers(
    tryCatch(list(val = expr, err = NULL),
             error = function(e) list(val = NULL, err = conditionMessage(e))),
    warning = function(c) {
      w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning")
    },
    message = function(c) {
      m <<- c(m, conditionMessage(c)); invokeRestart("muffleMessage")
    })
  c(r, list(warn = w, msg = m))
}
mk <- function(seed, G = 10, K = 6) {
  set.seed(seed)
  d <- expand.grid(t = 1:K, g = factor(1:G))
  d <- d[order(d$g, d$t), ]
  d$x <- rnorm(nrow(d))
  d$ft <- factor(d$t)
  e <- as.vector(replicate(G, arima.sim(list(ar = .4), K, sd = .25)))
  d$y <- rpois(nrow(d), exp(0.5 + 0.3 * d$x + e))
  d
}
cat("== M2 predict against ar1(), review's seeds\n")
for (s in c(1, 3, 7, 9, 13, 17, 19)) {
  d <- mk(s)
  fl <- suppressWarnings(frm(y ~ x + ar(t, g, cov = TRUE), family = poisson(), data = d))
  fr <- suppressWarnings(frm(y ~ x + ar1(ft + 0 | g), family = poisson(), data = d))
  pl <- suppressWarnings(predict(fl, ndraws = 500))[, 1]
  pr <- suppressWarnings(predict(fr, ndraws = 500))[, 1]
  nl <- suppressWarnings(predict(fl, newdata = d, ndraws = 500))[, 1]
  nr <- suppressWarnings(predict(fr, newdata = d, ndraws = 500))[, 1]
  cat(sprintf("seed %2d maxy %3d lat %.4g ar1 %.4g | new lat %.4g ar1 %.4g | sd %.3g\n",
              s, max(d$y), max(pl), max(pr), max(nl), max(nr),
              exp(fl$estimates$theta[1])))
}

cat("\n== M3 boundary, review's design\n")
for (s in 1:4) {
  set.seed(s)
  d <- expand.grid(t = 1:6, g = factor(1:10))
  d$x <- rnorm(nrow(d))
  d$y <- rpois(nrow(d), exp(0.5 + 0.3 * d$x + rep(rnorm(10, 0, .5), each = 6)))
  r <- catch(frm(y ~ x + ar(t, g, cov = TRUE), family = poisson(), data = d))
  f <- r$val
  cat("seed", s, "ar", summary(f)$cor_pars[1, 1], "se", summary(f)$cor_pars[1, 2], "\n")
  cat(" msg:", r$msg, "\n warn:", r$warn, "\n")
  if (any(ns$latres_at_limit(f$estimates$theta, ns$frame_latres(f$frame)[[1]]$aux))) {
    print(summary(f)$cor_pars)
    print(diagnose(f)$singular)
  }
}

cat("\n== M5 nlf(shape ~ ...) beside linear mu\n")
d <- mk(1)
d$ypos <- rgamma(nrow(d), 3, 3 / exp(0.5 + 0.3 * d$x))
r <- catch(frm(bf(ypos ~ x + ar(t, g, cov = TRUE), nlf(shape ~ exp(a)), a ~ 1),
               family = Gamma("log"), data = d))
print(r$err); print(r$warn); if (!is.null(r$val)) print(logLik(r$val))
r <- catch(frm(bf(ypos ~ exp(a) + ar(t, g, cov = TRUE), a ~ 1, nl = TRUE),
               family = Gamma("log"), data = d, dry_run = "frame"))
print(substr(r$err, 1, 120))

cat("\n== M6 fitted() at new rows, exp(0.5 + s^2/2)\n")
d <- mk(167, 30, 5)
fp <- suppressWarnings(frm(y ~ x + ar(t, g, cov = TRUE), family = poisson(), data = d))
nd <- data.frame(t = 1:5, g = factor(rep("new", 5)), x = 0)
for (s in c(0.25, 1, 2)) {
  fr <- fp
  fr$estimates$beta[] <- c(0.5, 0.4)
  fr$estimates$theta[] <- c(log(s), 0.5 / sqrt(1 - 0.25))
  ms <- s / sqrt(1 - 0.25)
  fv <- fitted(fr, newdata = nd, allow_new_levels = TRUE)
  lp <- frm_linpred(fr, newdata = nd, type = "response", allow_new_levels = TRUE)
  cat(sprintf("s %.2f marg %.3f exact %.5f fitted %.5f linpred %.5f se %s\n", s, ms,
              exp(0.5 + ms^2 / 2), fv[1, 1], lp[1], format(fv[1, 2], digits = 3)))
}

cat("\n== M7 skip\n")
set.seed(1)
d7 <- expand.grid(t = 1:40, g = factor(1:100)); d7$x <- rnorm(nrow(d7))
d7$y <- rpois(nrow(d7), exp(.5 + .3 * d7$x))
f7 <- suppressWarnings(frm(y ~ x + ar(t, g, cov = TRUE), family = poisson(), data = d7))
print(system.time(r <- catch(fitted(f7))))
cat(substr(r$msg, 1, 200), "\n"); print(head(r$val, 2))
print(system.time(r <- catch(frm_linpred(f7, se.fit = TRUE))))
cat(substr(r$msg, 1, 200), "\n"); print(head(r$val$se.fit, 2))

cat("\n== m1 coef\n")
d <- mk(1)
fm <- suppressWarnings(frm(y ~ x + ar(t, g, cov = TRUE), family = poisson(), data = d))
print(substr(catch(coef(fm))$err %||% "no error", 1, 120))
str(catch(coef(fm))$val)
cat("\n== m2 no time column\n")
print(catch(simulate(fm, newdata = d[, c("x", "g")]))$err)
print(catch(predict(fm, newdata = d[, c("x", "t")], ndraws = 10))$err)

cat("\n== m3 asym_laplace\n")
d$yc <- rnorm(nrow(d))
print(catch(frm(yc ~ x + ar(t, g, cov = TRUE), family = asym_laplace(), data = d, dry_run = "frame"))$err)

cat("\n== m4 all single-row groups\n")
d4 <- data.frame(t = 1:40, g = factor(1:40), x = rnorm(40)); d4$y <- rpois(40, 2)
print(catch(frm(y ~ x + ar(t, g, cov = TRUE), family = poisson(), data = d4, dry_run = "frame"))$err)
print(catch(frm(y ~ x + cosy(t, g), family = poisson(), data = d4, dry_run = "frame"))$err)
print(catch(frm(y ~ x + ar(t, g, cov = FALSE), family = poisson(), data = d4, dry_run = "frame"))$err)

cat("\n== m8 autocor_matrix under cov = FALSE\n")
fu <- suppressWarnings(frm(y ~ x + ar(t, g, cov = FALSE), family = poisson(), data = d))
print(round(catch(autocor_matrix(fu))$val, 3))
print(substr(catch(autocor_matrix(fu))$err %||% "", 1, 200))

cat("\n== m9 steps on one series\n")
d9 <- data.frame(t = 1:50000, x = 0)
a <- ns$latres_aux(list(struct = "ar", cov = TRUE, p = 1L, q = 0L, fn = "ar", label = "ar(t)"),
                   rep(1L, 50000L), 1:50000, 50000L, FALSE)
print(system.time(st <- ns$latres_steps(a)))
print(length(st))

cat("\n== n1 dirichlet bypass\n")
d$un <- rbeta(nrow(d), 2, 2)
print(substr(catch(frm(bf(cbind(un, 1 - un) ~ x + ar(t, g, cov = TRUE)), family = dirichlet(), data = d, dry_run = "frame"))$err, 1, 200))

cat("\n== n2 marginaleffects\n")
if (requireNamespace("marginaleffects", quietly = TRUE)) {
  r <- catch(marginaleffects::predictions(fu, newdata = d[1:3, ]))
  print(r$err); print(r$warn)
}
cat("\n== mi on Gamma\n")
d$ypos[3] <- NA
print(substr(catch(frm(bf(ypos | mi() ~ x + ar(t, g, cov = TRUE)), family = Gamma("log"), data = d, dry_run = "frame"))$err, 1, 200))
