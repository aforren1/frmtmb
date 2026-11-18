# Short check of lane sparsese: vcov(full = TRUE) of the random-slope design
# at seed 2 (sparsese-short-slope.R), whose correlation is lost at its
# boundary. Is the lost row NaN, as summary() says?
#   Rscript sparsese-short-vcovneg.R <lane|base>
arm <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-sparsese-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib", dirname(find.package("frmtmb")), "\n")
set.seed(2)
ng <- 30L
n <- 300L
d <- data.frame(x = stats::runif(n, 1, 2),
                g = factor(rep(seq_len(ng), length.out = n)),
                h = factor(sample(rep(seq_len(12), length.out = n))),
                f = factor(rep(c("a", "b"), each = 1, length.out = n)))
for (k in 1:7) d[[paste0("x", k)]] <- d$x^k
ug <- stats::rnorm(ng, 0, 0.5)
sg <- stats::rnorm(ng, 0, 0.4)
uh <- stats::rnorm(12, 0, 0.3)
eta <- sin(3 * d$x) + ug[d$g] + 0.3 * (d$f == "b")
d$y <- eta + sg[d$g] * d$x + stats::rnorm(n, 0, 0.3)
fit <- withCallingHandlers(
  frm(y ~ x1 + x2 + x3 + x4 + x5 + (1 + x1 | g), data = d, REML = TRUE),
  warning = function(w) {
    cat("W:", substr(conditionMessage(w), 1, 150), "\n")
    invokeRestart("muffleWarning")
  }, message = function(m) {
    cat("M:", substr(conditionMessage(m), 1, 150), "\n")
    invokeRestart("muffleMessage")
  })
s <- frmtmb:::sdr_of(fit)
cat("lost:", paste(names(s$se_lost), s$se_lost), " opt par:",
    format(fit$opt$par, digits = 6), "\n")
print(suppressWarnings(vcov(fit, full = TRUE)))
print(suppressWarnings(summary(fit))$random)
