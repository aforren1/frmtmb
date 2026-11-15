# Lane predraw2, fix round (review FIX 1): the draw check on GLMs with one
# coefficient per level of a factor (`y ~ 0 + f`) and few rows per level,
# truth known. Built from dev/reviews/predraw2-review-d.R. Per fit: the
# check's median excess, its ratio to the median drop the law implies
# (`rel`), the number of test draws and parameters, the verdict; where it
# fires, the coefficients' 90% coverage and the predictive 90% coverage
# (exact, against the true law) on this build and with draw_reweight()
# mocked to return its input (0.82.0's draws: no side law is built here),
# and the number of distinct parameter draws of 1000.
# Usage: Rscript dev/predraw2-kdesign.R <lib> <out.tsv> <seeds a:b>
#        [design ...]
args <- commandArgs(TRUE)
.libPaths(c(args[1],
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("lib:", dirname(find.package("frmtmb")), "\n")
out <- args[2]
sq <- as.integer(strsplit(args[3], ":")[[1]])
seeds <- seq(sq[1], sq[2])
real <- ns$draw_reweight
ident <- function(object, dw) dw$D
designs <- list(
  pois_K100_m5_mu1 = list(K = 100, m = 5, mu = 1, fam = "poisson"),
  pois_K200_m5_mu1 = list(K = 200, m = 5, mu = 1, fam = "poisson"),
  pois_K100_m5_mu1.5 = list(K = 100, m = 5, mu = 1.5, fam = "poisson"),
  pois_K200_m5_mu2 = list(K = 200, m = 5, mu = 2, fam = "poisson"),
  pois_K500_m5_mu2 = list(K = 500, m = 5, mu = 2, fam = "poisson"),
  pois_K200_m3_mu1.5 = list(K = 200, m = 3, mu = 1.5, fam = "poisson"),
  pois_K50_m5_mu2 = list(K = 50, m = 5, mu = 2, fam = "poisson"),
  pois_K100_m10_mu1.5 = list(K = 100, m = 10, mu = 1.5, fam = "poisson"),
  bern_K60_m20_p03 = list(K = 60, m = 20, mu = 0.3, fam = "bernoulli"),
  bern_K100_m10_p03 = list(K = 100, m = 10, mu = 0.3, fam = "bernoulli"),
  gaus_K10_m3 = list(K = 10, m = 3, mu = 1, fam = "gaussian"),
  gaus_K50_m3 = list(K = 50, m = 3, mu = 1, fam = "gaussian"),
  gaus_K200_m3 = list(K = 200, m = 3, mu = 1, fam = "gaussian"),
  gaus_K500_m3 = list(K = 500, m = 3, mu = 1, fam = "gaussian"),
  gamma_K10_m3 = list(K = 10, m = 3, mu = 1, fam = "gamma"),
  gamma_K50_m3 = list(K = 50, m = 3, mu = 1, fam = "gamma"),
  gamma_K200_m3 = list(K = 200, m = 3, mu = 1, fam = "gamma"),
  gamma_K500_m3 = list(K = 500, m = 3, mu = 1, fam = "gamma"),
  # sparser still: where each coordinate's skew alone is large
  gamma_K20_m2 = list(K = 20, m = 2, mu = 1, fam = "gamma"),
  gamma_K50_m2 = list(K = 50, m = 2, mu = 1, fam = "gamma"),
  pois_K20_m2_mu1 = list(K = 20, m = 2, mu = 1, fam = "poisson"),
  pois_K50_m3_mu1 = list(K = 50, m = 3, mu = 1, fam = "poisson"),
  bern_K20_m10_p03 = list(K = 20, m = 10, mu = 0.3, fam = "bernoulli"))
if (length(args) > 3) designs <- designs[args[-(1:3)]]
rows <- list()
for (nm in names(designs)) {
  ds <- designs[[nm]]
  for (sd in seeds) {
    set.seed(sd)
    d <- data.frame(f = factor(rep(seq_len(ds$K), each = ds$m)))
    if (ds$fam == "poisson") {
      tl <- log(ds$mu) + stats::rnorm(ds$K, 0, 0.3)
      lam <- exp(tl)[d$f]
      d$y <- stats::rpois(nrow(d), lam)
      fam <- poisson()
      cover <- function(lo, hi) {
        lo <- ceiling(lo - 1e-9)
        hi <- floor(hi + 1e-9)
        stats::ppois(hi, lam) - stats::ppois(lo - 1, lam)
      }
    } else if (ds$fam == "bernoulli") {
      tl <- stats::qlogis(ds$mu) + stats::rnorm(ds$K, 0, 0.3)
      pr <- stats::plogis(tl)[d$f]
      d$y <- stats::rbinom(nrow(d), 1, pr)
      fam <- bernoulli()
      cover <- function(lo, hi) {
        (lo <= 0 & hi >= 0) * (1 - pr) + (lo <= 1 & hi >= 1) * pr
      }
    } else if (ds$fam == "gaussian") {
      tl <- ds$mu + stats::rnorm(ds$K, 0, 1)
      mu <- tl[d$f]
      d$y <- stats::rnorm(nrow(d), mu, 1)
      fam <- gaussian()
      cover <- function(lo, hi) stats::pnorm(hi, mu, 1) - stats::pnorm(lo, mu, 1)
    } else {
      tl <- log(ds$mu) + stats::rnorm(ds$K, 0, 0.5)
      mu <- exp(tl)[d$f]
      d$y <- stats::rgamma(nrow(d), 2, 2 / mu)
      fam <- Gamma(link = "log")
      cover <- function(lo, hi) {
        stats::pgamma(hi, 2, 2 / mu) - stats::pgamma(lo, 2, 2 / mu)
      }
    }
    fit <- tryCatch(suppressWarnings(suppressMessages(
      frm(bf(y ~ 0 + f), family = fam, data = d))), error = function(e) NULL)
    if (is.null(fit)) {
      cat(nm, sd, "fit failed\n")
      next
    }
    said <- character(0)
    run <- function() {
      set.seed(1)
      withCallingHandlers(
        predict(fit, ndraws = 1000, summary = FALSE),
        warning = function(w) {
          said <<- c(said, class(w)[1])
          invokeRestart("muffleWarning")
        }, message = function(m) {
          said <<- c(said, class(m)[1])
          invokeRestart("muffleMessage")
        })
    }
    t1 <- proc.time()[[3]]
    pl <- run()
    t2 <- proc.time()[[3]]
    chk <- fit$cache[["draw_check"]]
    r <- data.frame(design = nm, seed = sd, np = length(fit$opt$par),
                    fire = isTRUE(chk$fire), excess = chk$excess %||% NA,
                    rel = chk$rel %||% NA, n = chk$n %||% NA,
                    fallback = any(grepl("could not be recovered", said)),
                    uniq = NA, coef_lane = NA, coef_base = NA,
                    pred_lane = NA, pred_base = NA, t_first = t2 - t1,
                    said = paste(unique(said), collapse = ","))
    if (isTRUE(chk$fire)) {
      assignInNamespace("draw_reweight", ident, ns = "frmtmb")
      pb <- run()
      assignInNamespace("draw_reweight", real, ns = "frmtmb")
      pdraw <- function() {
        set.seed(1)
        dr <- suppressWarnings(suppressMessages(
          ns$predict_par_drawer(fit, TRUE, 1000L)))
        sapply(1:1000, function(s) attr(dr, "delta")(s)) + fit$opt$par
      }
      Dl <- pdraw()
      assignInNamespace("draw_reweight", ident, ns = "frmtmb")
      Db <- pdraw()
      assignInNamespace("draw_reweight", real, ns = "frmtmb")
      kk <- seq_len(ds$K)
      pcov <- function(D) {
        lo <- apply(D[kk, , drop = FALSE], 1, stats::quantile, 0.05)
        hi <- apply(D[kk, , drop = FALSE], 1, stats::quantile, 0.95)
        mean(lo <= tl & tl <= hi)
      }
      cov <- function(p) {
        lo <- apply(p, 2, stats::quantile, 0.05)
        hi <- apply(p, 2, stats::quantile, 0.95)
        mean(cover(lo, hi))
      }
      r$uniq <- length(unique(Dl[1, ]))
      r$coef_lane <- pcov(Dl)
      r$coef_base <- pcov(Db)
      r$pred_lane <- cov(pl)
      r$pred_base <- cov(pb)
    }
    rows[[length(rows) + 1L]] <- r
    cat(sprintf(paste("%-20s seed %2d np %3d fire %-5s excess %7.3f rel",
                      "%7.4f n %2d fb %-5s | uniq %s coef %s/%s pred %s/%s |",
                      "t %.1f\n"), nm, sd, r$np, r$fire, r$excess, r$rel,
                r$n, r$fallback, r$uniq, format(r$coef_lane, digits = 3),
                format(r$coef_base, digits = 3),
                format(r$pred_lane, digits = 4),
                format(r$pred_base, digits = 4), r$t_first))
  }
}
utils::write.table(do.call(rbind, rows), out, sep = "\t", row.names = FALSE,
                   quote = FALSE)
