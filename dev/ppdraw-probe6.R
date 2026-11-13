# The flat side of draw_edge_scan(): for every parameter of the sweep's
# fits (dev/ppdraw-sweep-gen.R, seed 101), the backlog fits and the core
# suite's one interior hold (test-re-formula-partial.R), the smaller of
# the two log-likelihood changes at +-2 draw-law standard deviations.
# An estimate AT the end of its range changes by rounding there; one
# merely near it changes by more. Prints one line per parameter whose
# smaller change is below 0.5.
# Usage: Rscript dev/ppdraw-probe6.R <lib>
args <- commandArgs(TRUE)
.libPaths(c(args[1], "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
n <- 300L
source("dev/ppdraw-sweep-gen.R")
flat <- function(fit, label) {
  ds <- frmtmb:::fit_draw_space(fit)
  if (is.null(ds$V) || !all(is.finite(ds$V))) {
    cat(sprintf("%-40s no usable covariance\n", label))
    return(invisible())
  }
  p <- fit$opt$par
  off <- length(ds$map$names) - length(p)
  sd <- sqrt(diag(ds$V))
  f0 <- fit$obj$fn(p)
  for (j in seq_along(p)) {
    s <- sd[off + j]
    if (!is.finite(s) || s <= 0) next
    dr <- vapply(c(-1, 1), function(k) {
      q <- p
      q[j] <- p[j] + 2 * k * s
      tryCatch(suppressWarnings(fit$obj$fn(q)) - f0,
               error = function(e) NA_real_)
    }, 0)
    m <- suppressWarnings(min(abs(dr[is.finite(dr)])))
    if (is.finite(m) && m < 0.5) {
      cat(sprintf("%-40s %-42s est %9.3f sd %9.3g flat side %.3g\n",
                  label, rownames(ds$V)[off + j], p[j], s, m))
    }
  }
}
q <- function(expr) suppressMessages(suppressWarnings(expr))
for (fam in fams) for (case in c("interior", "boundary")) {
  s <- spec(fam)
  fit <- q(tryCatch(frm(s[[1]], family = s[[2]], data = gen(fam, case)),
                    error = function(e) NULL))
  if (!is.null(fit)) flat(fit, paste(fam, case))
}
# 3269 and its two neighbors
mk <- function(zi) {
  set.seed(7)
  d <- data.frame(x = stats::rnorm(200), N = sample(5:15, 200, TRUE))
  d$y <- stats::rbinom(200, d$N, stats::plogis(0.3 + 0.5 * d$x))
  if (zi > 0) d$y[stats::runif(200) < zi] <- 0L
  d
}
for (zi in c(0, 0.02, 0.3)) {
  flat(q(frm(bf(y | trials(N) ~ x), family = zero_inflated_binomial(),
             data = mk(zi))), paste("3269 zib zi", zi))
}
# 3277
set.seed(7)
d <- data.frame(x = stats::rnorm(400))
k <- stats::rbinom(400, 1, 0.3)
d$yo <- cut(d$x + ifelse(k == 1, 2, 0) + stats::rlogis(400),
            c(-Inf, -1, 0.5, 2, Inf), labels = FALSE)
flat(q(frm(bf(yo ~ x), family = mixture(cumulative(), cumulative()),
           data = d)), "3277 ordmix")
# se(s, sigma = TRUE)
set.seed(7)
dg <- data.frame(x = rnorm(200), s = runif(200, 0.2, 1))
dg$y <- 1 + dg$x + rnorm(200, 0, dg$s)
flat(q(frm(bf(y | se(s, sigma = TRUE) ~ x), family = gaussian(), data = dg)),
     "se(sigma = TRUE)")
# test-re-formula-partial.R
set.seed(3)
G <- 15
m <- 8
dp <- data.frame(g = factor(rep(seq_len(G), each = m)),
                 h = factor(rep(1:5, times = G * m / 5)),
                 x = stats::rnorm(G * m))
ug <- stats::rnorm(G, 0, 0.8)
sg <- stats::rnorm(G, 0, 0.4)
uh <- stats::rnorm(5, 0, 0.6)
dp$y <- 1 + 0.5 * dp$x + ug[dp$g] + sg[dp$g] * dp$x + uh[dp$h] +
  stats::rnorm(G * m)
flat(q(frm(bf(y ~ x + (1 + x | g) + (1 | h)), data = dp)), "re-formula-partial")
# random slopes with a small true sd, which lme4 does not call singular
for (sd_s in c(0.05, 0.1, 0.2)) for (seed in 1:5) {
  set.seed(seed)
  dd <- data.frame(g = factor(rep(1:20, each = 10)), x = stats::rnorm(200))
  dd$y <- 1 + 0.5 * dd$x + stats::rnorm(20, 0, 0.7)[dd$g] +
    stats::rnorm(20, 0, sd_s)[dd$g] * dd$x + stats::rnorm(200)
  flat(q(frm(bf(y ~ x + (1 + x | g)), data = dd)),
       sprintf("slope sd %.2f seed %d", sd_s, seed))
}
cat("done\n")
