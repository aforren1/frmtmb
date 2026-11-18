# Short check of lane drawdir, check 5: predict()'s draws, seeded, on the
# lane build and on rellib-r22 (0.83.0), for fits the per-direction
# reading does not fire on and for a mixture the summed check fires on
# (0.83.0's reweighting path). Each arm saves its estimates and draws;
# the second arm run compares. Run pinned to one core type
# (dev/rtmb-pitfalls.md item 21).
# Usage: Rscript dev/reviews/drawdir-short-bit.R lane|base
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-drawdir-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
log <- file(sprintf("dev/reviews/drawdir-short-bit-%s.txt", arm), "w")
sink(log)
sink(log, type = "message")
suppressMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
q <- function(expr) suppressMessages(suppressWarnings(expr))
cls <- function(expr) {
  out <- character(0)
  val <- withCallingHandlers(expr,
    warning = function(w) {
      out <<- c(out, class(w)[1])
      invokeRestart("muffleWarning")
    }, message = function(m) {
      out <<- c(out, class(m)[1])
      invokeRestart("muffleMessage")
    })
  list(val = val, cls = unique(out))
}
fits <- list(
  grp = function() {
    set.seed(1)
    d <- data.frame(g = factor(rep(1:30, each = 4)), x = stats::rnorm(120))
    d$y <- 1 + 0.5 * d$x + stats::rnorm(30, 0, 0.3)[d$g] + stats::rnorm(120)
    q(frm(bf(y ~ x + (1 | g)), data = d))
  },
  pois = function() {
    set.seed(3)
    d <- data.frame(x = stats::rnorm(200))
    d$y <- stats::rpois(200, exp(1 + 0.3 * d$x))
    q(frm(bf(y ~ x), family = poisson(), data = d))
  },
  sleep = function() {
    q(frm(bf(Reaction ~ Days + (Days | Subject)), data = lme4::sleepstudy))
  },
  slope = function() {
    set.seed(2003)
    d <- data.frame(g = factor(rep(1:20, each = 8)), x = stats::rnorm(160))
    d$y <- 1 + 0.5 * d$x + stats::rnorm(20, 0, 0.8)[d$g] +
      stats::rnorm(20, 0, 0.5)[d$g] * d$x + stats::rnorm(160)
    q(frm(bf(y ~ x + (1 + x | g)), data = d))
  },
  mix = function() {
    set.seed(101)
    x <- stats::rnorm(300)
    invisible(stats::runif(300))
    d <- data.frame(x = x, y = stats::rpois(300, exp(1 + 0.3 * x)))
    q(frm(bf(y ~ x), family = mixture(poisson(), poisson()), data = d))
  })
out <- list()
for (nm in names(fits)) {
  fit <- fits[[nm]]()
  set.seed(1)
  p1 <- cls(predict(fit, ndraws = 200, summary = FALSE))
  nd <- fit$data[1:3, , drop = FALSE]
  set.seed(7)
  p2 <- cls(predict(fit, newdata = nd, ndraws = 100, summary = FALSE))
  chk <- fit$cache[["draw_check"]]
  out[[nm]] <- list(par = fit$opt$par, p1 = p1$val, p2 = p2$val,
                    cls = union(p1$cls, p2$cls))
  cat(sprintf("%-6s fire %s excess %.4g dir %s warn_only %s | said %s\n",
              nm, chk$fire, chk$excess, format(chk$dir),
              format(chk$warn_only), paste(out[[nm]]$cls, collapse = ",")))
}
saveRDS(out, sprintf("dev/reviews/drawdir-short-bit-%s.rds", arm))
other <- sprintf("dev/reviews/drawdir-short-bit-%s.rds",
                 if (arm == "lane") "base" else "lane")
if (file.exists(other)) {
  o <- readRDS(other)
  for (nm in names(out)) {
    cat(sprintf("%-6s par identical %s | in-sample draws identical %s | newdata draws identical %s | max |diff| %.3g\n",
                nm, identical(out[[nm]]$par, o[[nm]]$par),
                identical(out[[nm]]$p1, o[[nm]]$p1),
                identical(out[[nm]]$p2, o[[nm]]$p2),
                max(abs(out[[nm]]$p1 - o[[nm]]$p1))))
  }
}
cat("DONE\n")
sink(type = "message")
sink()
