# Each new test assertion over 20 data sets moved by a few ulps,
# v * (1 + 2e-15 k), k = 1..20 (lane-env-r083.md). Cases:
#   limit  item 4867, the covariate x moved
#   sdesc  items 4836/4861, the response moved
#   mix    item 5222, the response moved
#   latln  item 5620, the lognormal response moved (both spellings)
#   latbern item 5620's guard-present case, the covariate moved
#   Rscript dev/quadseam2-ulp.R <libpath|base> <case> [k range]
a <- commandArgs(TRUE)
.libPaths(c(if (a[1] != "base") a[1], "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
case <- a[2]
ks <- if (length(a) > 2) eval(parse(text = a[3])) else 1:20
mv <- function(v, k) v * (1 + 2e-15 * k)
conds <- function(expr) {
  w <- character(); m <- character()
  v <- withCallingHandlers(expr, warning = function(x) {
    w <<- c(w, paste(class(x)[1], conditionMessage(x)))
    invokeRestart("muffleWarning")
  }, message = function(x) {
    m <<- c(m, class(x)[1]); invokeRestart("muffleMessage")
  })
  list(value = v, w = w, m = m)
}
ok_all <- 0L
for (k in ks) {
  ok <- switch(case,
    limit = {
      set.seed(5)
      d <- data.frame(g = factor(rep(1:20, each = 9)),
                      f = factor(rep(c("a", "b", "c"), 60)), x = rnorm(180))
      d$y <- rpois(180, exp(0.5 + 0.3 * d$x + rnorm(20, 0, 0.5)[d$g]))
      d$y[d$f == "c"] <- 0L
      d$x <- mv(d$x, k)
      r <- conds(frm(y ~ f + x + (1 | g), family = poisson(), data = d,
                     quadrature = TRUE))
      v <- diag(sdr_of(r$value)$cov.fixed)
      length(r$w) == 1L && grepl("toward -Inf", r$w) &&
        identical(unname(sdr_of(r$value)$se_lost), "limit") &&
        is.nan(v[[3L]]) && all(is.finite(v[-3L]))
    },
    sdesc = {
      set.seed(6); n <- 200
      xs <- -abs(rnorm(n)) * 3
      y <- xs + (abs(rnorm(n)) - sqrt(2 / pi)) * 1.5 + rnorm(n, 0, 0.3)
      y <- mv(y, k)
      dd <- data.frame(y = y, xs = xs, g = factor(rep(1:20, length.out = n)))
      m3 <- mean((y - mean(y))^3) / sd(y)^3
      bad <- list(betad = c(log(sd(y)), 2 * sign(m3) + 0.5 * m3))
      fo <- bf(y ~ xs + (1 | g), sigma ~ 1, alpha ~ 1)
      r0 <- conds(frm(fo, family = skew_normal(), data = dd,
                      quadrature = TRUE))
      r1 <- conds(frm(fo, family = skew_normal(), data = dd,
                      quadrature = TRUE, start = bad))
      ll0 <- as.numeric(logLik(r0$value))
      g <- r1$value$opt[["sd_escape"]][["gain"]]
      cat(sprintf("   k %d gain %s gap %.3e w %d m %d\n", k,
                  format(g, digits = 4),
                  abs(as.numeric(logLik(r1$value)) - ll0),
                  length(r1$w), sum(r1$m == "frmtmb_boundary_fit")))
      isTRUE(g > 1e-2) &&
        abs(as.numeric(logLik(r1$value)) - ll0) / abs(ll0) < 1e-7 &&
        !any(r1$m == "frmtmb_boundary_fit") &&
        is.null(r0$value$opt[["sd_escape"]]) && !length(r0$w) && !length(r1$w)
    },
    mix = {
      sn <- function(n, mu, sigma, alpha) {
        dl <- alpha / sqrt(1 + alpha^2)
        mu + sigma * (dl * abs(rnorm(n)) + sqrt(1 - dl^2) * rnorm(n))
      }
      set.seed(101); n <- 400
      z <- rbinom(n, 1, 0.5)
      dd <- data.frame(y = mv(ifelse(z == 1, sn(n, -3, 1, 5), sn(n, 3, 1, -5)),
                              k))
      mf <- mixture(skew_normal(), skew_normal())
      u <- frm(bf(y ~ 1), family = mf, data = dd, dry_run = "objective")
      free <- conds(frm(bf(y ~ 1), family = mf, data = dd))
      bd <- u$estimates$betad
      al <- grep("^alpha", names(bd))
      bd[al] <- 0
      st <- conds(frm(bf(y ~ 1), family = mf, data = dd,
                      start = list(betad = bd)))
      ll <- as.numeric(logLik(free$value))
      all(abs(st$value$estimates$betad[al]) > 1) &&
        abs(as.numeric(logLik(st$value)) - ll) < abs(ll) * 1e-8 &&
        identical(unname(st$value$opt$stationary_escape[["starts"]]), 6) &&
        !length(st$w)
    },
    latln = {
      set.seed(40401)
      d <- expand.grid(t = 1:6, g = factor(1:10))
      d <- d[order(d$g, d$t), ]
      d$x <- rnorm(nrow(d))
      e <- as.vector(replicate(10, arima.sim(list(ar = 0.6), 6, sd = 0.5)))
      d$y <- mv(exp(0.6 + 0.4 * d$x + e + rnorm(nrow(d), 0, 0.3)), k)
      d$y[sort(sample(nrow(d), 8))] <- NA
      dw <- d; dw$w <- ifelse(is.na(d$y), 0, 1); dw$y[is.na(dw$y)] <- 1
      r1 <- conds(frm(bf(y | mi() ~ x + ar(t, g, cov = TRUE)),
                      family = lognormal(), data = d))
      r2 <- conds(frm(bf(y | weights(w) ~ x + ar(t, g, cov = TRUE)),
                      family = lognormal(), data = dw))
      !any(grepl("^frmtmb_latres_laplace", c(r1$w, r2$w)))
    },
    latbern = {
      set.seed(502)
      d <- expand.grid(t = 1:10, g = factor(1:20))
      d <- d[order(d$g, d$t), ]
      d$x <- rnorm(nrow(d))
      e <- as.vector(replicate(20, arima.sim(list(ar = 0.5), 10, sd = 0.8)))
      d$yp <- rpois(nrow(d), exp(0.3 * d$x + e))
      d$yb <- rbinom(nrow(d), 1, plogis(0.3 * d$x + e))
      d$x <- mv(d$x, k)
      r <- conds(frm(bf(yb ~ x + ar(t, g, cov = TRUE)), family = bernoulli(),
                     data = d))
      sum(grepl("^frmtmb_latres_laplace", r$w)) == 1L
    })
  cat(sprintf("K %d %s\n", k, if (isTRUE(ok)) "PASS" else "FAIL"))
  ok_all <- ok_all + isTRUE(ok)
}
cat(sprintf("ULP %s pass %d of %d\n", case, ok_all, length(ks)))
