# Reviewer: frmtmb ML fits against brms 2.23.0 posterior means on the
# same censored or truncated data, seven families. A discrete right or
# interval lower bound is handed to brms less one (frmtmb's inclusive
# convention; NEWS says so), and also handed unchanged, to show the
# size of the convention difference. Output: dev/cdfs-rev-out/brmsfit.txt
.libPaths(c("C:/Users/adf44/source/r/wt-cdfs-rev-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages({
  library(brms)
  library(frmtmb)
})
options(frmtmb.notices = FALSE, mc.cores = 1)
bf <- frmtmb::bf
set.seed(20261021)
n <- 300
x <- rnorm(n)
cases <- list()
# student: left below q10, right above q85, 10% interval of width 0.6
y <- 1 + 0.5 * x + 0.8 * rt(n, 5)
d <- data.frame(x = x, y = y, cc = 0, y2 = y + 0.6)
d$cc[y < quantile(y, 0.1)] <- -1
d$cc[y > quantile(y, 0.85)] <- 1
d$cc[sample(which(d$cc == 0), 30)] <- 2
cases$student <- list(d = d, f = frmtmb::student(), b = brms::student())
# Gamma: right above q80, left below q10
y <- rgamma(n, 3, 3 / exp(0.5 + 0.3 * x))
d <- data.frame(x = x, y = y, cc = 0, y2 = y)
d$cc[y > quantile(y, 0.8)] <- 1
d$y[d$cc == 1] <- quantile(y, 0.8)
d$cc[y < quantile(y, 0.1)] <- -1
d$y[d$cc == -1] <- quantile(y, 0.1)
cases$Gamma <- list(d = d, f = stats::Gamma(link = "log"), b = Gamma(link = "log"))
# von_mises: left below -1.2, right above 1.4
y <- 2 * atan(0.4 * x) + rnorm(n, 0, 0.5)
y <- pmin(pmax(y, -pi + 0.01), pi - 0.01)
d <- data.frame(x = x, y = y, cc = 0, y2 = y)
d$cc[y < -1.2] <- -1; d$y[d$cc == -1] <- -1.2
d$cc[y > 1.4] <- 1; d$y[d$cc == 1] <- 1.4
cases$von_mises <- list(d = d, f = frmtmb::von_mises(), b = brms::von_mises())
# hurdle_lognormal: right above q85
y <- ifelse(runif(n) < 0.25, 0, rlnorm(n, 0.3 + 0.4 * x, 0.6))
d <- data.frame(x = x, y = y, cc = 0, y2 = y)
q85 <- quantile(y[y > 0], 0.85)
d$cc[y > q85] <- 1; d$y[d$cc == 1] <- q85
cases$hurdle_lognormal <- list(d = d, f = frmtmb::hurdle_lognormal(),
                               b = brms::hurdle_lognormal())
# negbinomial: right-censored "8 or more", left-censored "1 or fewer"
y <- rnbinom(n, mu = exp(1.2 + 0.4 * x), size = 2.5)
d <- data.frame(x = x, y = pmin(pmax(y, 1), 8), cc = 0)
d$cc[y >= 8] <- 1
d$cc[y <= 1] <- -1
d$y2 <- d$y
cases$negbinomial <- list(d = d, f = frmtmb::negbinomial(), b = brms::negbinomial(),
                          discrete = TRUE)
# zero_inflated_poisson: right-censored "6 or more"
y <- rpois(n, exp(1 + 0.3 * x)) * rbinom(n, 1, 0.75)
d <- data.frame(x = x, y = pmin(y, 6), cc = as.numeric(y >= 6))
d$y2 <- d$y
cases$zero_inflated_poisson <- list(d = d, f = frmtmb::zero_inflated_poisson(),
                                    b = brms::zero_inflated_poisson(), discrete = TRUE)
# beta_binomial: interval-censored into bins [0,2], [3,5], [6,8], [9,12]
# on 40% of rows
y <- rbinom(n, 12, rbeta(n, 3 * plogis(0.4 * x) * 2, 3 * (1 - plogis(0.4 * x)) * 2))
d <- data.frame(x = x, y = y, tr = 12, cc = 0)
lo <- c(0, 3, 6, 9)[findInterval(y, c(0, 3, 6, 9))]
hi <- c(2, 5, 8, 12)[findInterval(y, c(0, 3, 6, 9))]
ib <- sample(n, 120)
d$cc[ib] <- 2
d$y2 <- d$y
d$y[ib] <- lo[ib]
d$y2[ib] <- hi[ib]
cases$beta_binomial <- list(d = d, f = frmtmb::beta_binomial(), b = brms::beta_binomial(),
                            discrete = TRUE, trials = TRUE)
# a truncated fit: Gamma truncated to [0.5, 6]
y <- rgamma(3 * n, 3, 3 / exp(0.5 + 0.3 * rep(x, 3)))
xx <- rep(x, 3)
keep <- y >= 0.5 & y <= 6
d <- data.frame(x = xx[keep], y = y[keep])[1:n, ]
# fix round 1: exgaussian, right above q85, at sigma / beta about 2 and
# about 25 (the regime RTMBdist::dexgauss() approximated)
for (ub in c(2, 25)) {
  y <- 1 + 0.4 * x + rnorm(n, 0, 0.5) + rexp(n, ub / 0.5)
  d <- data.frame(x = x, y = y, cc = 0, y2 = y)
  q85 <- quantile(y, 0.85)
  d$cc[y > q85] <- 1; d$y[d$cc == 1] <- q85
  cases[[paste0("exgaussian_u", ub)]] <- list(d = d, f = frmtmb::exgaussian(),
                                              b = brms::exgaussian())
}
# fix round 1: poisson, interval bins inside trunc(lb = 1, ub = 20)
y <- rpois(3 * n, exp(1.3 + 0.3 * rep(x, 3)))
xx <- rep(x, 3)
keep <- which(y >= 1 & y <= 20)[1:n]
d <- data.frame(x = xx[keep], y = y[keep], cc = 0)
d$y2 <- d$y
ib <- sample(n, 100)
d$cc[ib] <- 2
d$y[ib] <- pmax(d$y2[ib] - 2, 1)
d$cc[d$cc == 2 & d$y2 <= d$y] <- 0
cases$poisson_trunc_interval <- list(d = d, f = poisson(), b = poisson(),
                                     discrete = TRUE, tr2 = TRUE)
cases$Gamma_trunc <- list(d = d, f = stats::Gamma(link = "log"), b = Gamma(link = "log"),
                          trunc = TRUE)

summ_frm <- function(fit) {
  s <- summary(fit)
  rbind(s$fixed[, c("Estimate", "Est.Error")], s$spec_pars[, c("Estimate", "Est.Error")])
}
summ_brm <- function(bfit) {
  s <- summary(bfit)
  rbind(s$fixed[, c("Estimate", "Est.Error")], s$spec_pars[, c("Estimate", "Est.Error")])
}
sel <- strsplit(Sys.getenv("CASES", paste(names(cases), collapse = ",")), ",")[[1]]
for (nm in sel) {
  cs <- cases[[nm]]
  d <- cs$d
  lhs <- if (isTRUE(cs$tr2)) "y | cens(cc, y2) + trunc(lb = 1, ub = 20)" else
    if (isTRUE(cs$trunc)) "y | trunc(lb = 0.5, ub = 6)" else
    if (isTRUE(cs$trials)) "y | trials(tr) + cens(cc, y2)" else "y | cens(cc, y2)"
  fo <- as.formula(paste(lhs, "~ x"))
  fit <- withCallingHandlers(frm(fo, family = cs$f, data = d),
                             warning = function(w) {
                               cat("frm WARN:", conditionMessage(w), "\n")
                               invokeRestart("muffleWarning")
                             })
  sf <- summ_frm(fit)
  db <- d
  if (isTRUE(cs$discrete) && Sys.getenv("NOSHIFT") != "1") {
    sh <- db$cc %in% c(1, 2)
    db$y[sh] <- db$y[sh] - 1
    # brms takes no response below 0: an interval [0, k2] is the left
    # censored row Y <= k2 there
    z0 <- db$cc == 2 & db$y < 0
    db$cc[z0] <- -1
    db$y[z0] <- db$y2[z0]
  }
  t0 <- proc.time()[["elapsed"]]
  bfit <- brm(fo, data = db, family = cs$b, chains = 2, iter = 2000, warmup = 1000,
              seed = 7, refresh = 0, silent = 2, backend = "rstan")
  sb <- summ_brm(bfit)
  cat(sprintf("\n== %s (frmtmb convergence %d, logLik %.4f; brms %.0f s, max Rhat %.3f)\n", nm,
              fit$opt$convergence, as.numeric(logLik(fit)),
              proc.time()[["elapsed"]] - t0, max(rhat(bfit), na.rm = TRUE)))
  rn <- intersect(rownames(sf), rownames(sb))
  for (r in rn) {
    cat(sprintf("  %-12s frmtmb %10.5f (SE %8.5f) | brms %10.5f (SD %8.5f) | z %6.2f\n", r,
                sf[r, 1], sf[r, 2], sb[r, 1], sb[r, 2], (sf[r, 1] - sb[r, 1]) / sb[r, 2]))
  }
  missing <- setdiff(union(rownames(sf), rownames(sb)), rn)
  if (length(missing)) cat("  unmatched rows:", missing, "\n")
  sc <- strsplit(stancode(bfit), "\n")[[1]]
  cat("  brms censored / truncated rows:\n")
  cat(paste0("    ", trimws(grep("lccdf|lcdf|log_diff_exp", sc, value = TRUE))), sep = "\n")
  if (isTRUE(cs$discrete) && Sys.getenv("NOSHIFT") != "1") {
    # the same data handed to brms unshifted: the convention difference
    bfit2 <- update(bfit, newdata = d, refresh = 0, silent = 2, seed = 7)
    sb2 <- summ_brm(bfit2)
    cat("  brms on the UNSHIFTED bounds:\n")
    for (r in rn) {
      cat(sprintf("  %-12s frmtmb %10.5f | brms %10.5f (SD %8.5f) | z %6.2f\n", r,
                  sf[r, 1], sb2[r, 1], sb2[r, 2], (sf[r, 1] - sb2[r, 1]) / sb2[r, 2]))
    }
  }
}
