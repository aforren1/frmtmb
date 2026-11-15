# Review of lane predraw2: the draw check's excess on every fit of
# test-predict-draw-check.R and test-predict-draw-edge.R that reaches it,
# over 21 roundings of each data set (x times 1 + 2e-15 k, k = 0..20),
# with the thresholds' distances. One R process.
# Usage: Rscript predraw2-review-b.R <label> [kmax]
local({
.libPaths(c("C:/Users/adf44/source/r/wt-predraw2-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
a <- commandArgs(trailingOnly = TRUE)
label <- a[1]
kmax <- if (length(a) > 1) as.integer(a[2]) else 20L
cat("lib:", dirname(find.package("frmtmb")), "label", label,
    "coretype", Sys.getenv("OPENBLAS_CORETYPE"), "\n")
print(La_library())
q <- function(expr) suppressWarnings(suppressMessages(expr))

zib <- function(zi, seed = 7) {
  set.seed(seed)
  n <- 200
  d <- data.frame(x = stats::rnorm(n), N = sample(5:15, n, TRUE))
  d$y <- stats::rbinom(n, d$N, stats::plogis(0.3 + 0.5 * d$x))
  if (zi > 0) d$y[stats::runif(n) < zi] <- 0L
  d
}
cases <- list(
  C1_mixpois = function() {
    set.seed(101)
    x <- stats::rnorm(300)
    invisible(stats::runif(300))
    d <- data.frame(x = x, y = stats::rpois(300, exp(1 + 0.3 * x)))
    list(d = d, f = function(d) frm(bf(y ~ x), family = mixture(poisson(), poisson()), data = d))
  },
  C2_glm = function() {
    set.seed(3)
    d <- data.frame(x = stats::rnorm(100))
    d$y <- stats::rpois(100, exp(0.5 + 0.4 * d$x))
    list(d = d, f = function(d) frm(bf(y ~ x), family = poisson(), data = d))
  },
  E1_zib0 = function() list(d = zib(0), f = function(d)
    frm(bf(y | trials(N) ~ x), family = zero_inflated_binomial(), data = d)),
  E2_zib03 = function() list(d = zib(0.3), f = function(d)
    frm(bf(y | trials(N) ~ x), family = zero_inflated_binomial(), data = d)),
  E3_zib002 = function() list(d = zib(0.02), f = function(d)
    frm(bf(y | trials(N) ~ x), family = zero_inflated_binomial(), data = d)),
  E4_zibx = function() list(d = zib(0), f = function(d)
    frm(bf(y | trials(N) ~ x, zi ~ x), family = zero_inflated_binomial(), data = d)),
  E5_hux = function() {
    set.seed(12)
    dh <- data.frame(x = stats::rnorm(300))
    dh$y <- stats::rpois(300, exp(2.5 + 0.2 * dh$x))
    list(d = dh, f = function(d) frm(bf(y ~ x, hu ~ x), family = hurdle_poisson(), data = d))
  },
  E6_zip = function() {
    set.seed(3)
    d <- data.frame(x = stats::rnorm(200))
    d$y <- stats::rpois(200, exp(1 + 0.3 * d$x))
    list(d = d, f = function(d) frm(bf(y ~ x), family = zero_inflated_poisson(), data = d))
  },
  E7_gsd = function() {
    set.seed(5)
    d <- data.frame(g = factor(rep(1:6, each = 20)), x = stats::rnorm(120))
    d$y <- 1 + 0.5 * d$x + stats::rnorm(6, 0, 0.1)[d$g] + stats::rnorm(120)
    list(d = d, f = function(d) frm(bf(y ~ x + (1 | g)), data = d))
  },
  E8_sesig = function() {
    set.seed(7)
    d <- data.frame(x = stats::rnorm(200), s = stats::runif(200, 0.2, 1))
    d$y <- 1 + d$x + stats::rnorm(200, 0, d$s)
    list(d = d, f = function(d) frm(bf(y | se(s, sigma = TRUE) ~ x), family = gaussian(), data = d))
  },
  E9_slope_h = function() {
    set.seed(3)
    G <- 15; m <- 8
    d <- data.frame(g = factor(rep(seq_len(G), each = m)),
                    h = factor(rep(1:5, times = G * m / 5)),
                    x = stats::rnorm(G * m))
    ug <- stats::rnorm(G, 0, 0.8); sg <- stats::rnorm(G, 0, 0.4)
    uh <- stats::rnorm(5, 0, 0.6)
    d$y <- 1 + 0.5 * d$x + ug[d$g] + sg[d$g] * d$x + uh[d$h] + stats::rnorm(G * m)
    list(d = d, f = function(d) frm(bf(y ~ x + (1 + x | g) + (1 | h)), data = d))
  },
  E14_sleep = function() {
    set.seed(11)
    d <- data.frame(g = factor(rep(1:18, each = 10)), x = rep(0:9, 18))
    d$y <- 250 + 10 * d$x + stats::rnorm(18, 0, 25)[d$g] +
      stats::rnorm(18, 0, 6)[d$g] * d$x + stats::rnorm(180, 0, 25)
    # x is a count here: perturb y
    list(d = d, f = function(d) frm(bf(y ~ x + (x | g)), data = d), pert = "y")
  },
  E16_ordmix = function() {
    set.seed(7)
    d <- data.frame(x = stats::rnorm(400))
    k <- stats::rbinom(400, 1, 0.3)
    d$yo <- cut(d$x + ifelse(k == 1, 2, 0) + stats::rlogis(400),
                c(-Inf, -1, 0.5, 2, Inf), labels = FALSE)
    list(d = d, f = function(d) frm(bf(yo ~ x), family = mixture(cumulative(), cumulative()), data = d))
  }
)
for (lam in c(0.5, 0.6, 0.7, 1)) local({
  lam <- lam
  cases[[paste0("E10_lam", lam)]] <<- function() {
    set.seed(5)
    dd <- data.frame(g = factor(rep(1:20, each = 10)), x = stats::rnorm(200))
    a <- stats::rnorm(20, 0, 0.7)[dd$g]
    b <- stats::rnorm(20, 0, 0.05)[dd$g] * dd$x
    e <- stats::rnorm(200)
    dd$y <- 1 + 0.5 * dd$x + a + lam * b + e
    list(d = dd, f = function(d) frm(bf(y ~ x + (1 + x | g)), data = d))
  }
})
sel <- Sys.getenv("REVB_CASES")
if (nzchar(sel)) cases <- cases[strsplit(sel, ",")[[1]]]
for (nm in names(cases)) {
  cs <- cases[[nm]]()
  ex <- rep(NA_real_, kmax + 1L); nn <- rep(NA_integer_, kmax + 1L)
  rr <- rep(NA_character_, kmax + 1L)
  for (k in 0:kmax) {
    d <- cs$d
    col <- cs$pert %||% "x"
    d[[col]] <- d[[col]] * (1 + 2e-15 * k)
    fit <- tryCatch(q(cs$f(d)), error = function(e) NULL)
    if (is.null(fit)) next
    set.seed(1)
    dr <- tryCatch(q(ns$predict_par_drawer(fit, TRUE, 2L)),
                   error = function(e) NULL)
    chk <- fit$cache[["draw_check"]]
    if (!is.null(chk)) { ex[k + 1L] <- chk$excess; nn[k + 1L] <- chk$n }
    law <- fit$cache[["draw_edge_law"]]
    if (is.list(law)) {
      steep <- function(v) length(v) > 1L
      r <- unlist(lapply(c(law$lo, law$hi), function(v) if (steep(v)) v[2] / 2))
      if (length(r)) rr[k + 1L] <- paste(format(sort(r), digits = 3), collapse = ",")
    }
  }
  thr <- c(0.5, 1, 2)
  dist <- if (all(is.na(ex))) NA else min(abs(outer(ex[!is.na(ex)], thr, "-")))
  cat(sprintf("%-12s excess min %9.4f max %9.4f | n %s | fire %s | min dist to {0.5,1,2} %.4f | r of steep sides (k=0) %s\n",
              nm, min(ex, na.rm = TRUE), max(ex, na.rm = TRUE),
              paste(unique(nn), collapse = "/"),
              paste(unique(ex > 1), collapse = "/"), dist, rr[1]))
  cat("   ", nm, "per k:", format(ex, digits = 6), "\n")
}
})
