# Short check of lane quadseam2's fix round, F2: latres_quadratic() with a
# step of at least 1 on the link scale. Attack: families whose support
# on the link scale is bounded so that a step of 1 leaves it (Beta with
# identity and log links, von_mises identity near pi, Gamma inverse),
# and logit Beta with its mean near 0 and near 1. Expected: the gate
# answers FALSE, so the warning decision is 0.82.0's, and no warning
# escapes from inside the gate. Also lognormal at y x 1e-8 and 1e8, both
# spellings of an unread mi() response: gate TRUE, no latres warning.
# Per fit: gate (lane only), count of frmtmb_latres_laplace, other
# warning classes, logLik.
#   Rscript quadseam2-short-latres.R <lib|base>
a <- commandArgs(TRUE)
if (a[1] == "lib") a[1] <- "C:/Users/adf44/source/r/wt-quadseam2-lib"
.libPaths(c(if (a[1] != "base") a[1], "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
ns <- asNamespace("frmtmb")
has_q <- exists("latres_quadratic", ns)
# warnings raised while the gate itself runs, outside frm()
gate <- function(f) {
  if (!has_q) return("-")
  bk <- ns$frame_latres(f$frame)[[1L]]
  gw <- character()
  g <- withCallingHandlers(
    tryCatch(ns$latres_quadratic(f, bk), error = function(e) "ERR"),
    warning = function(x) {
      gw <<- c(gw, class(x)[1]); invokeRestart("muffleWarning")
    })
  paste0(g, if (length(gw)) paste0("[gatewarn:", paste(gw, collapse = ","),
                                   "]"))
}
one <- function(lbl, ...) {
  w <- character()
  f <- tryCatch(withCallingHandlers(suppressMessages(frm(...)),
    warning = function(x) {
      w <<- c(w, class(x)[1]); invokeRestart("muffleWarning")
    }), error = function(e) e)
  if (inherits(f, "error")) {
    cat(sprintf("ROW %-10s ERROR %s\n", lbl,
                substr(gsub("\\s+", " ", conditionMessage(f)), 1, 110)))
    return(invisible())
  }
  cat(sprintf("ROW %-10s ll %.6f latres %d gate %s other %s\n", lbl,
              as.numeric(logLik(f)), sum(w == "frmtmb_latres_laplace"),
              gate(f), paste(setdiff(w, "frmtmb_latres_laplace"),
                             collapse = ",")))
}
gen <- function(seed, G = 25L, K = 6L, sd = 0.5) {
  set.seed(seed)
  d <- expand.grid(t = seq_len(K), g = factor(seq_len(G)))
  d <- d[order(d$g, d$t), ]
  d$x <- rnorm(nrow(d))
  d$e <- as.vector(replicate(G, arima.sim(list(ar = 0.5), K, sd = sd)))
  d
}
rb <- function(mu, phi = 30) {
  pmin(pmax(rbeta(length(mu), mu * phi, (1 - mu) * phi), 1e-6), 1 - 1e-6)
}
fo <- function(v) bf(as.formula(paste(v, "~ x + ar(t, g, cov = TRUE)")))
for (seed in 1:3) {
  d <- gen(900 + seed)
  # logit Beta, mean near 0 and near 1
  d$blo <- rb(plogis(-3.5 + 0.2 * d$x + d$e))
  d$bhi <- rb(plogis(3.5 + 0.2 * d$x + d$e))
  # identity Beta at a mean of 0.5 (a step of 1 leaves (0, 1) both ways)
  d$bid <- rb(0.5 + 0.05 * d$x + 0.1 * d$e)
  # log Beta at a mean of 0.9 (a step of +1 leaves (0, 1))
  d$blg <- rb(0.9 * exp(0.02 * d$x + 0.05 * d$e - 0.05))
  # von_mises: tan_half link, and identity near pi
  mu_vm <- 2 * atan(0.3 * d$x + d$e)
  d$vm <- (mu_vm + rnorm(nrow(d), 0, 0.3) + pi) %% (2 * pi) - pi
  d$vmi <- pmin(pmax(2.6 + 0.05 * d$x + 0.1 * d$e +
                       rnorm(nrow(d), 0, 0.1), -pi + 1e-6), pi - 1e-6)
  # Gamma, inverse link at eta = 2 (test-latres.R's case)
  d$gi <- rgamma(nrow(d), 50, 50 * (2 + 0.1 * d$x + 0.3 * d$e / 0.5))
  one(paste0("BETAlo/", seed), fo("blo"), family = Beta(), data = d)
  one(paste0("BETAhi/", seed), fo("bhi"), family = Beta(), data = d)
  one(paste0("BETAid/", seed), fo("bid"), family = Beta("identity"),
      data = d)
  one(paste0("BETAlog/", seed), fo("blg"), family = Beta("log"), data = d)
  one(paste0("VMtan/", seed), fo("vm"), family = von_mises(), data = d)
  one(paste0("VMid/", seed), fo("vmi"), family = von_mises("identity"),
      data = d)
  one(paste0("GAMinv/", seed), fo("gi"), family = Gamma(), data = d)
}
ln_data <- function(seed) {
  set.seed(seed)
  d <- expand.grid(t = seq_len(6L), g = factor(seq_len(10L)))
  d <- d[order(d$g, d$t), ]
  d$x <- rnorm(nrow(d))
  e <- as.vector(replicate(10L, arima.sim(list(ar = 0.6), 6L, sd = 0.5)))
  d$y <- exp(0.6 + 0.4 * d$x + e + rnorm(nrow(d), 0, 0.3))
  d$y[sort(sample(nrow(d), 8))] <- NA
  d
}
for (sc in c(1e-8, 1e8)) for (seed in 40401:40403) {
  d <- ln_data(seed)
  d$y <- d$y * sc
  one(sprintf("LNmi/%g/%d", sc, seed),
      bf(y | mi() ~ x + ar(t, g, cov = TRUE)), family = lognormal(), data = d)
  dw <- d; dw$w <- ifelse(is.na(d$y), 0, 1); dw$y[is.na(dw$y)] <- sc
  one(sprintf("LNw0/%g/%d", sc, seed),
      bf(y | weights(w) ~ x + ar(t, g, cov = TRUE)), family = lognormal(),
      data = dw)
}
cat("END\n")
