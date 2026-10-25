# A brms 2.23.0 model with bhaz() and frmtmb's spelling of it are the
# same model. For three bases:
#   1. frmtmb's Zbhaz and Zcbhaz (its frame, dry_run = "frame") against
#      brms's standata();
#   2. the log likelihood at each of brms's posterior draws: brms's
#      log_lik() summed over rows, against frmtmb's cox() lpdf and
#      lccdf at the same draw, on frmtmb's own basis;
#   3. frmtmb's maximum likelihood fit: the bhaz() spelling and the
#      cox() spelling are identical(), and no brms draw has a higher
#      log likelihood than frmtmb's maximum.
# Needs a Stan compile:
#   Rscript dev/bhaz-brms-fit.R > dev/bhaz-log/brms-fit.txt 2>&1
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
.libPaths(c("C:/Users/adf44/source/r/wt-bhaz-lib",
            "C:/Users/adf44/source/r/rellib-r12",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(brms)
  library(frmtmb.survival)
})
cat("frmtmb", format(packageVersion("frmtmb")), "at",
    find.package("frmtmb"), "\nfrmtmb.survival",
    format(packageVersion("frmtmb.survival")), "at",
    find.package("frmtmb.survival"), "\nbrms",
    format(packageVersion("brms")), " rstan",
    format(packageVersion("rstan")), "\n")
set.seed(11); n <- 200; x <- rnorm(n)
tt <- rexp(n, exp(-0.5 + 0.7 * x)); ct <- rexp(n, 0.2)
dd <- data.frame(time = pmin(tt, ct), cens = as.numeric(tt > ct), x = x)
cat("events", sum(dd$cens == 0), "of", n, "\n")

models <- list(
  df6 = list(
    fo = time | cens(cens) + bhaz(df = 6) ~ x,
    frm_cox = frmtmb.survival::cox(df = 6)),
  knots = list(
    fo = time | cens(cens) + bhaz(knots = c(0.5, 1.1),
                                  Boundary.knots = c(0, 16)) ~ x,
    frm_cox = frmtmb.survival::cox(knots = c(0.5, 1.1),
                                   Boundary.knots = c(0, 16))),
  deg2 = list(
    fo = time | cens(cens) + bhaz(df = 5, degree = 2,
                                  intercept = FALSE) ~ x,
    frm_cox = frmtmb.survival::cox(df = 5, degree = 2, intercept = FALSE))
)

strip <- function(m) matrix(as.numeric(m), nrow(m))
for (nm in names(models)) {
  m <- models[[nm]]
  cat("\n===== ", nm, ": ", deparse1(m$fo), " =====\n", sep = "")
  sd <- brms::standata(brms::bf(m$fo), data = dd, family = brms::cox())
  fr <- frm(frmtmb::bf(m$fo), family = frmtmb.survival::cox(), data = dd,
            dry_run = "frame")
  av <- fr$aterm_values$time
  cat(sprintf(paste0("basis: Kbhaz brms %d frmtmb %d  max|dZbhaz|/",
                     "max|Zbhaz| %.3e  max|dZcbhaz|/max|Zcbhaz| %.3e\n"),
              sd$Kbhaz, ncol(av$Zbhaz),
              max(abs(av$Zbhaz - strip(sd$Zbhaz))) / max(abs(sd$Zbhaz)),
              max(abs(av$Zcbhaz - strip(sd$Zcbhaz))) /
                max(abs(sd$Zcbhaz))))

  # the two frmtmb spellings
  f_bhaz <- frm(frmtmb::bf(m$fo), family = frmtmb.survival::cox(),
                data = dd)
  f_cox <- frm(frmtmb::bf(time | cens(cens) ~ x), family = m$frm_cox,
               data = dd)
  cat("frmtmb logLik, bhaz() spelling",
      format(as.numeric(logLik(f_bhaz)), digits = 17),
      " cox() spelling", format(as.numeric(logLik(f_cox)), digits = 17),
      "\n identical logLik", identical(logLik(f_bhaz), logLik(f_cox)),
      " identical estimates", identical(f_bhaz$estimates, f_cox$estimates),
      " identical vcov", identical(vcov(f_bhaz), vcov(f_cox)), "\n")

  # brms
  t0 <- proc.time()[[3]]
  bfit <- brms::brm(brms::bf(m$fo), data = dd, family = brms::cox(),
                    chains = 2, iter = 1000, seed = 1, refresh = 0,
                    backend = "rstan", silent = 2)
  cat("brms fit in", round(proc.time()[[3]] - t0), "s\n")
  ll <- brms::log_lik(bfit)
  brms_ll <- rowSums(ll)
  dr <- posterior::as_draws_df(bfit)
  K <- sd$Kbhaz
  S <- as.matrix(dr[, paste0("sbhaz[", seq_len(K), "]")])
  eta <- outer(dr$b_Intercept, rep(1, n)) + outer(dr$b_x, dd$x)
  fam <- stats::family(f_bhaz)
  aterms <- av
  aterms$sbhaz_chart <- "softmax"
  ev <- dd$cens == 0
  frm_ll <- vapply(seq_len(nrow(S)), function(d) {
    raw <- log(S[d, -1L]) - log(S[d, 1L])
    mu <- exp(eta[d, ])
    dens <- fam$lpdf(dd$time, list(mu = mu), aterms,
                     list(sbhaz_raw = raw))
    surv <- fam$lccdf(dd$time, list(mu = mu), aterms,
                      list(sbhaz_raw = raw))
    sum(dens[ev]) + sum(surv[!ev])
  }, 0)
  rel <- abs(frm_ll - brms_ll) / abs(brms_ll)
  cat(sprintf(paste0("draws %d  max |frmtmb - brms| / |brms| log lik ",
                     "%.3e  median %.3e\n"),
              length(brms_ll), max(rel), stats::median(rel)))
  mle <- as.numeric(logLik(f_bhaz))
  cat(sprintf(paste0("frmtmb MLE logLik %.10f  max over brms draws ",
                     "%.10f  draws above the MLE %d\n"),
              mle, max(brms_ll), sum(brms_ll > mle)))
  fx <- fixef(f_bhaz)
  cat(sprintf(paste0("b_x: frmtmb MLE %.4f (SE %.4f)  brms posterior ",
                     "mean %.4f (sd %.4f)\n"),
              fx["x", 1], fx["x", 2], mean(dr$b_x), stats::sd(dr$b_x)))
  cat("Rhat max", format(max(brms::rhat(bfit), na.rm = TRUE), digits = 4),
      "\n")
}
