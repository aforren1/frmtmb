# Review of lane arspat: "nothing else reads" detection under the ar()
# family (mvbf with mi(y) as a predictor, partial reads, subset()), the
# latent-residual families, post-fit methods on an mi() fit, and the gap
# message's false positives and noise. Run on the lane build, and with
# ARSPAT_LIB="" on rellib-r20 where a "before" is meaningful.
source("dev/reviews/arspat-review-common.R")

set.seed(23)
G <- 10; K <- 6
d <- expand.grid(t = seq_len(K), g = factor(seq_len(G)))
n <- nrow(d)
d$x <- rnorm(n)
e <- as.vector(vapply(seq_len(G), function(i)
  as.vector(arima.sim(list(ar = 0.6), K, sd = 0.8)), numeric(K)))
d$y <- 1 + 0.5 * d$x + e
d$z <- 0.5 + 0.8 * d$y + rnorm(n, 0, 0.5)
d$sub <- d$g %in% as.character(1:7)
d$xm <- d$x
d$xm[c(5, 17, 40)] <- NA
d$yp <- exp(d$y / 2)
miss_y <- c(3, 4, 14, 31, 36, 55)
miss_z <- c(4, 22, 36)
dm <- d
dm$y[miss_y] <- NA
dmz <- dm
dmz$z[miss_z] <- NA

# dense gaussian reference for y | mi() ~ x + ar() and z ~ mi(y)
joint_ref <- function(fit, dd) {
  fe <- fixef(fit)[, 1]
  print(fe)
  dp <- ns$eval_dpars(fit)
  nat <- ac_nat(fit)
  phi <- unname(nat[1])
  sy <- dp[[1]]$sigma[1]; sz <- dp[[2]]$sigma[1]
  b <- fe[grep("^y_", names(fe))]; a <- fe[grep("^z_", names(fe))]
  muy <- b[1] + b[2] * dd$x
  Sy <- matrix(0, n, n)
  for (i in split(seq_len(n), dd$g)) {
    Sy[i, i] <- sy^2 / (1 - phi^2) * phi^abs(outer(dd$t[i], dd$t[i], "-"))
  }
  muz <- a[1] + a[2] * muy
  mu <- c(muy, muz)
  S <- rbind(cbind(Sy, a[2] * Sy), cbind(a[2] * Sy, a[2]^2 * Sy + diag(sz^2, n)))
  v <- c(dd$y, dd$z)
  ok <- !is.na(v)
  mvtnorm::dmvnorm(v[ok], mu[ok], S[ok, ok], log = TRUE)
}

chk("mvbf: z ~ mi(y), every missing y read", {
  f <- frm(bf(y | mi() ~ x + ar(t, g, cov = TRUE)) + bf(z ~ mi(y)) +
             set_rescor(FALSE), data = dm)
  cat("unread:", f$frame$mi_map[[1]]$unread, "\n")
  ll <- as.numeric(logLik(f)); ref <- joint_ref(f, dm)
  cat(sprintf("logLik %.8f ref %.8f diff %.3e conv %d\n", ll, ref, ll - ref,
              f$opt$convergence))
  NULL
})
chk("mvbf: z | mi() ~ mi(y), some missing y unread", {
  f <- frm(bf(y | mi() ~ x + ar(t, g, cov = TRUE)) + bf(z | mi() ~ mi(y)) +
             set_rescor(FALSE), data = dmz)
  cat("unread y:", f$frame$mi_map[["y"]]$unread, " unread z:",
      f$frame$mi_map[["z"]]$unread, "\n")
  ll <- as.numeric(logLik(f)); ref <- joint_ref(f, dmz)
  cat(sprintf("logLik %.8f ref %.8f diff %.3e conv %d\n", ll, ref, ll - ref,
              f$opt$convergence))
  NULL
})
chk("mvbf student y read by z", {
  f <- frm(bf(y | mi() ~ x + ar(t, g, cov = TRUE), family = student()) +
             bf(z ~ mi(y)) + set_rescor(FALSE), data = dm)
  cat("logLik", as.numeric(logLik(f)), "conv", f$opt$convergence, "\n")
  cat("summary mentions Laplace/approx? ",
      any(grepl("Laplace|approx", capture.output(summary(f)))), "\n")
  NULL
})
chk("subset(): y on a subset with mi()", {
  f <- frm(bf(y | mi() + subset(sub) ~ x + ar(t, g, cov = TRUE)) +
             bf(z ~ x) + set_rescor(FALSE), data = dm)
  f1 <- frm(bf(y | mi() ~ x + ar(t, g, cov = TRUE)), data = dm[dm$sub, ])
  f2 <- frm(z ~ x, data = dm)
  cat(sprintf("joint %.8f separate %.8f diff %.3e\n", as.numeric(logLik(f)),
              as.numeric(logLik(f1)) + as.numeric(logLik(f2)),
              as.numeric(logLik(f)) - as.numeric(logLik(f1)) -
                as.numeric(logLik(f2))))
  NULL
})
chk("mi() on a predictor, ar on y (before and after must agree)", {
  f <- frm(bf(y ~ x + mi(xm) + ar(t, g, cov = TRUE)) + bf(xm | mi() ~ 1) +
             set_rescor(FALSE), data = d)
  cat(sprintf("logLik %.10f\n", as.numeric(logLik(f))))
  NULL
})
chk("lognormal ar cov=TRUE with mi() (before and after must agree)", {
  dl <- d; dl$yp[miss_y] <- NA
  f <- frm(bf(yp | mi() ~ x + ar(t, g, cov = TRUE)), family = lognormal(),
           data = dl)
  cat(sprintf("logLik %.10f conv %d\n", as.numeric(logLik(f)),
              f$opt$convergence))
  NULL
})
chk("plain ar/ma/arma/cosy/unstr gaussian (before and after must agree)", {
  for (tm in c("ar(t, g, cov = TRUE)", "ma(t, g, cov = TRUE)",
               "arma(t, g, cov = TRUE)", "cosy(t, g)", "unstr(t, g)")) {
    f <- frm(bf(stats::as.formula(paste("y ~ x +", tm))), data = d)
    ft <- frm(bf(stats::as.formula(paste("y ~ x +", tm))), data = d,
              family = student())
    cat(sprintf("%-24s gaussian %.10f student %.10f\n", tm,
                as.numeric(logLik(f)), as.numeric(logLik(ft))))
  }
  NULL
})

## post-fit methods on an ar() + mi() fit
fm <- tryCatch(frm(bf(y | mi() ~ x + ar(t, g, cov = TRUE)), data = dm),
               error = function(e) NULL)
if (!is.null(fm)) {
  chk("fitted rows", dim(fitted(fm)))
  chk("predict head", head(predict(fm, ndraws = 200), 5))
  chk("residuals at missing rows", residuals(fm)[miss_y, ])
  chk("simulate dims", dim(simulate(fm, nsim = 3, seed = 1)))
  chk("nobs", nobs(fm))
  chk("vcov_cluster refused under mi()", vcov_cluster(fm, ~ g))
  chk("predict newdata", head(predict(fm, newdata = d[1:12, ], ndraws = 100), 3))
}

## gap message: false positives and its class
gapcase <- function(label, dd, tm = "ar(t, g, cov = TRUE)", ...) {
  chk(label, {
    cls <- NULL
    f <- withCallingHandlers(
      frm(bf(stats::as.formula(paste("y ~ x +", tm))), data = dd, ...),
      message = function(m) {
        cls <<- c(cls, paste(class(m), collapse = "/"))
      })
    cat("message classes:", if (is.null(cls)) "none" else cls, "\n")
    NULL
  })
}
d_stag <- d[!(d$g == "1" & d$t %in% 1:2) & !(d$g == "2" & d$t %in% 5:6), ]
gapcase("staggered start/end, no interior gap (expect none)", d_stag)
d_gap <- d[!(d$g == "3" & d$t == 3), ]
gapcase("interior gap in group 3 (expect one)", d_gap)
gapcase("interior gap under cov = FALSE (expect none)", d_gap, "ar(t, g)")
gapcase("interior gap under cosy (expect none)", d_gap, "cosy(t, g)")
d_f <- d_gap; d_f$t <- factor(letters[d_f$t])
gapcase("factor time, interior gap (expect one)", d_f)
d_r <- d_gap; d_r$t <- d_r$t * 0.37 + 0.1
gapcase("non-integer time, interior gap (expect one)", d_r)
d_one <- d
d_one$t[d_one$g == "1" & d_one$t == 6] <- 5.5
gapcase("one off-grid time in group 1 (others skip 5.5?)", d_one)
d_na <- d; d_na$x[c(15, 16)] <- NA
gapcase("NA predictor drops rows (expect one)", d_na)
d_mi <- dm
chk("mi() rows keep the series whole (expect none)", {
  cls <- NULL
  withCallingHandlers(frm(bf(y | mi() ~ x + ar(t, g, cov = TRUE)), data = d_mi),
                      message = function(m) cls <<- c(cls, class(m)[1]))
  cat("message classes:", if (is.null(cls)) "none" else cls, "\n")
  NULL
})
chk("gap message on predict/simulate/bootstrap calls", {
  cnt <- 0L
  withCallingHandlers({
    f <- frm(y ~ x + ar(t, g, cov = TRUE), data = d_gap)
    cat("after fit:", cnt, "\n")
    invisible(predict(f, ndraws = 20))
    cat("after predict:", cnt, "\n")
    invisible(predict(f, newdata = d_gap, ndraws = 20))
    cat("after predict(newdata):", cnt, "\n")
    invisible(simulate(f, nsim = 2, newdata = d_gap, seed = 1))
    cat("after simulate(newdata):", cnt, "\n")
    invisible(update(f, data = d_gap))
    cat("after update:", cnt, "\n")
    b <- tryCatch(frm_bootstrap(f, R = 5, seed = 1), error = function(e)
      conditionMessage(e))
    cat("after frm_bootstrap(R = 5):", cnt, "\n")
  }, frmtmb_autocor_gap = function(m) {
    cnt <<- cnt + 1L
    invokeRestart("muffleMessage")
  })
  NULL
})
