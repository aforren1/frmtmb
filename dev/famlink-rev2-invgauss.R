## Recheck (rounds 1 and 1b), priority 4: does inverse.gaussian's new
## default link, 1/mu^2, converge across designs beyond the worker's
## four? Lane only. For each design, 20 replicate data sets (seeds
## 20260916 + r): the default-link fit and, as a control, the same data
## on link = "log". Recorded: status, optimizer code, max |gradient| at
## the optimum, pdHess, nonfinite_trials. For fixed-effect designs the
## default-link coefficients are also compared with
## stats::glm(family = inverse.gaussian()), whose canonical link is the
## same 1/mu^2 (the relative difference is reported; it is a reference
## implementation, not a tolerance).
## Usage: Rscript dev/famlink-rev2-invgauss.R
ARM <- "lane"
source("dev/famlink-rev-common.R")
rig <- function(n, mu, lambda) {
  # Michael, Schucany and Haas (1976)
  nu <- rnorm(n); y <- nu^2
  x <- mu + mu^2 * y / (2 * lambda) - mu / (2 * lambda) * sqrt(4 * mu * lambda * y + mu^2 * y^2)
  z <- runif(n)
  ifelse(z <= mu / (mu + x), x, mu^2 / x)
}
designs <- list(
  fixed_small_mean = function() { n <- 200; d <- data.frame(x = rnorm(n)); d$y <- rig(n, exp(-0.5 + 0.3 * d$x), 2); list(y ~ x, d, TRUE) },
  fixed_large_mean = function() { n <- 200; d <- data.frame(x = rnorm(n)); d$y <- rig(n, exp(3 + 0.4 * d$x), 50); list(y ~ x, d, TRUE) },
  fixed_steep = function() { n <- 300; d <- data.frame(x = runif(n, -2, 2)); d$y <- rig(n, exp(0.5 + 1.2 * d$x), 5); list(y ~ x, d, TRUE) },
  factor_pred = function() { n <- 240; d <- data.frame(f = factor(rep(c("a", "b", "c", "d"), 60)), x = rnorm(n)); d$y <- rig(n, exp(c(0, 0.5, 1, 1.5)[d$f] + 0.2 * d$x), 3); list(y ~ f + x, d, TRUE) },
  random_intercept = function() { n <- 300; d <- data.frame(x = rnorm(n), g = factor(rep(1:30, 10))); d$y <- rig(n, exp(1 + 0.3 * d$x + rnorm(30, 0, 0.4)[d$g]), 4); list(y ~ x + (1 | g), d, FALSE) },
  random_slope = function() { n <- 400; d <- data.frame(x = rnorm(n), g = factor(rep(1:40, 10))); d$y <- rig(n, exp(1 + (0.3 + rnorm(40, 0, 0.2)[d$g]) * d$x + rnorm(40, 0, 0.3)[d$g]), 4); list(y ~ x + (x | g), d, FALSE) },
  distributional_shape = function() { n <- 300; d <- data.frame(x = rnorm(n)); d$y <- rig(n, exp(1 + 0.3 * d$x), exp(1 + 0.5 * d$x)); list(bf(y ~ x, shape ~ x), d, FALSE) },
  weights = function() { n <- 200; d <- data.frame(x = rnorm(n), w = runif(n, 0.5, 2)); d$y <- rig(n, exp(0.8 + 0.3 * d$x), 3); list(y | weights(w) ~ x, d, FALSE) },
  censored = function() { n <- 250; d <- data.frame(x = rnorm(n)); y <- rig(n, exp(0.8 + 0.3 * d$x), 3); d$cen <- ifelse(y > 4, "right", "none"); d$y <- pmin(y, 4); list(y | cens(cen) ~ x, d, FALSE) },
  truncated = function() { n <- 250; d <- data.frame(x = rnorm(n)); y <- rig(3 * n, exp(0.8 + 0.3 * rep(d$x, 3)), 3); keep <- which(y > 0.3)[seq_len(n)]; d$y <- y[keep]; d$x <- rep(d$x, 3)[keep]; list(y | trunc(lb = 0.3) ~ x, d, FALSE) },
  intercept_only = function() { n <- 60; d <- data.frame(y = rig(60, 2, 4)); list(y ~ 1, d, TRUE) },
  small_n = function() { n <- 25; d <- data.frame(x = rnorm(n)); d$y <- rig(n, exp(0.5 + 0.5 * d$x), 3); list(y ~ x, d, TRUE) }
)
rows <- list()
for (nm in names(designs)) {
  for (r in 1:20) {
    set.seed(20260916 + r)
    spec <- designs[[nm]]()
    for (lk in c("default", "log")) {
      fam <- if (lk == "default") "inverse.gaussian" else brmsfamily("inverse.gaussian", "log")
      w <- character()
      fit <- withCallingHandlers(tryCatch(frm(spec[[1]], data = spec[[2]], family = fam, se = TRUE), error = function(e) e),
        warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") },
        message = function(m) invokeRestart("muffleMessage"))
      row <- data.frame(design = nm, rep = r, link = lk, status = if (inherits(fit, "error")) "error" else "fit",
                        conv = NA_integer_, grad = NA_real_, pdHess = NA, nonfinite = NA_integer_,
                        warnings = length(w), glm_reldiff = NA_real_,
                        msg = if (inherits(fit, "error")) substr(conditionMessage(fit), 1, 80) else substr(paste(w, collapse = " | "), 1, 80))
      if (!inherits(fit, "error")) {
        row$conv <- fit$opt$convergence
        row$grad <- tryCatch(max(abs(fit$obj$gr(fit$opt$par))), error = function(e) NA_real_)
        row$pdHess <- isTRUE(fit$cache$sdr$pdHess)
        row$nonfinite <- fit$opt$nonfinite_trials %||% NA_integer_
        if (lk == "default" && isTRUE(spec[[3]])) {
          g <- tryCatch(suppressWarnings(stats::glm(spec[[1]], data = spec[[2]], family = stats::inverse.gaussian())), error = function(e) NULL)
          if (!is.null(g) && g$converged) {
            b <- fit$estimates$beta
            row$glm_reldiff <- max(abs(b - stats::coef(g)) / pmax(abs(stats::coef(g)), 1e-8))
          }
        }
      }
      rows[[length(rows) + 1L]] <- row
    }
  }
}
tab <- do.call(rbind, rows)
saveRDS(tab, "dev/famlink-rev2-invgauss.rds")
agg <- do.call(rbind, lapply(split(tab, list(tab$design, tab$link), drop = TRUE), function(s) data.frame(
  design = s$design[1], link = s$link[1], n = nrow(s), fitted = sum(s$status == "fit"),
  conv0 = sum(s$conv == 0, na.rm = TRUE), pdHess = sum(s$pdHess, na.rm = TRUE),
  grad_max = signif(max(s$grad, na.rm = TRUE), 3), with_warnings = sum(s$warnings > 0),
  nonfinite_med = stats::median(s$nonfinite, na.rm = TRUE),
  glm_reldiff_max = signif(suppressWarnings(max(s$glm_reldiff, na.rm = TRUE)), 3))))
print(agg[order(agg$design, agg$link), ], row.names = FALSE)
bad <- tab[tab$status != "fit" | (!is.na(tab$conv) & tab$conv != 0) | tab$warnings > 0, ]
if (nrow(bad)) { cat("\nnon-clean replicates:\n"); print(bad[, c("design", "rep", "link", "status", "conv", "grad", "pdHess", "msg")], row.names = FALSE) }
