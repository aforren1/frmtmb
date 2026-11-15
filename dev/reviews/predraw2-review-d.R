# Review of lane predraw2: the draw check on ordinary GLMs with many
# fixed-effect levels and few counts per level. Each level's likelihood is
# a little skewed; the paired excess adds those up over the levels. Truth
# known, so the exact coverage of predict()'s 90% interval on the training
# rows is read, lane against the same call with draw_reweight() mocked out
# (which is 0.82.0's draw where no side law is built).
# Usage: Rscript predraw2-review-d.R [design ...]
local({
.libPaths(c("C:/Users/adf44/source/r/wt-predraw2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("lib:", dirname(find.package("frmtmb")), "\n")
real <- ns$draw_reweight
ident <- function(object, dw) dw$D
designs <- list(
  pois_K200_m5_mu1 = list(K = 200, m = 5, mu = 1, fam = "poisson"),
  pois_K100_m5_mu1.5 = list(K = 100, m = 5, mu = 1.5, fam = "poisson"),
  pois_K200_m5_mu2 = list(K = 200, m = 5, mu = 2, fam = "poisson"),
  pois_K100_m10_mu1.5 = list(K = 100, m = 10, mu = 1.5, fam = "poisson"),
  pois_K100_m5_mu1 = list(K = 100, m = 5, mu = 1, fam = "poisson"),
  pois_K50_m5_mu2 = list(K = 50, m = 5, mu = 2, fam = "poisson"),
  pois_K30_m10_mu3 = list(K = 30, m = 10, mu = 3, fam = "poisson"),
  bern_K60_m20_p03 = list(K = 60, m = 20, mu = 0.3, fam = "bernoulli")
)
a <- commandArgs(trailingOnly = TRUE)
if (length(a)) designs <- designs[a]
seeds <- as.integer(strsplit(Sys.getenv("REVD_SEEDS", "1:5"), ":")[[1]])
seeds <- seq(seeds[1], seeds[2])
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
      cover <- function(lo, hi) stats::ppois(hi, lam) - stats::ppois(lo - 1, lam)
    } else {
      tl <- stats::qlogis(ds$mu) + stats::rnorm(ds$K, 0, 0.3)
      pr <- stats::plogis(tl)[d$f]
      d$y <- stats::rbinom(nrow(d), 1, pr)
      fam <- bernoulli()
      cover <- function(lo, hi) {
        (lo <= 0 & hi >= 0) * (1 - pr) + (lo <= 1 & hi >= 1) * pr
      }
    }
    t0 <- proc.time()[[3]]
    fit <- tryCatch(suppressWarnings(suppressMessages(
      frm(bf(y ~ 0 + f), family = fam, data = d))), error = function(e) NULL)
    if (is.null(fit)) { cat(nm, sd, "fit failed\n"); next }
    tf <- proc.time()[[3]] - t0
    said <- character(0)
    run <- function() {
      set.seed(1)
      withCallingHandlers(
        predict(fit, ndraws = 1000, probs = c(0.05, 0.95), summary = FALSE),
        warning = function(w) {
          said <<- c(said, class(w)[1]); invokeRestart("muffleWarning")
        }, message = function(m) {
          said <<- c(said, class(m)[1]); invokeRestart("muffleMessage")
        })
    }
    t1 <- proc.time()[[3]]
    pl <- run()
    t2 <- proc.time()[[3]]
    chk <- fit$cache[["draw_check"]]
    assignInNamespace("draw_reweight", ident, ns = "frmtmb")
    pb <- run()
    assignInNamespace("draw_reweight", real, ns = "frmtmb")
    # the parameter draws behind each: their unique count, and how often
    # the 90% interval of each level's coefficient holds the true value
    pdraw <- function() {
      set.seed(1)
      dr <- suppressWarnings(suppressMessages(
        ns$predict_par_drawer(fit, TRUE, 1000L)))
      if (is.null(attr(dr, "delta"))) return(NULL)
      sapply(1:1000, function(s) attr(dr, "delta")(s)) + fit$opt$par
    }
    Dl <- pdraw()
    assignInNamespace("draw_reweight", ident, ns = "frmtmb")
    Db <- pdraw()
    assignInNamespace("draw_reweight", real, ns = "frmtmb")
    if (is.null(Dl)) { cat(nm, sd, "plug-in fallback (a level with no counts?)
"); next }
    nuniq <- length(unique(Dl[1, ]))
    pcov <- function(D) {
      lo <- apply(D, 1, stats::quantile, 0.05)
      hi <- apply(D, 1, stats::quantile, 0.95)
      mean(lo <= tl & tl <= hi)
    }
    cov <- function(p) {
      lo <- apply(p, 2, stats::quantile, 0.05, type = 7)
      hi <- apply(p, 2, stats::quantile, 0.95, type = 7)
      mean(cover(ceiling(lo - 1e-9), floor(hi + 1e-9)))
    }
    # the spread of the level means across draws, against the base draw's
    sdl <- mean(apply(pl, 2, stats::sd)); sdb <- mean(apply(pb, 2, stats::sd))
    cat(sprintf("%-18s seed %2d np %3d fire %-5s excess %7.3f n %2d | unique draws %4d | coef cov90 lane %.3f base %.3f | pred cov90 lane %.4f base %.4f | mean row sd lane %.3f base %.3f | said %s | t fit %.1f first predict %.1f\n",
                nm, sd, length(fit$opt$par), format(chk$fire), chk$excess,
                chk$n, nuniq, pcov(Dl), pcov(Db), cov(pl), cov(pb), sdl, sdb,
                paste(unique(said), collapse = ","), tf, t2 - t1))
  }
}
})
