# Review of lane gpgap, item 2: gr = FALSE with duplicated x, through the lane's
# translator (a row-by-row map), both programs. The hand map of handmap.txt
# matched by coordinate and so cannot tell two rows at one x apart. Seed 7301.
.libPaths(c("C:/Users/adf44/source/r/wt-gpgap-lib", "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(testthat)})
H <- new.env(parent = asNamespace("frmtmb"))
for (h in c("helper-brms.R", "helper-stan.R")) sys.source(file.path("tests/testthat", h), envir = H)
set.seed(7301)
n <- 40
d <- data.frame(x = sample(seq(0, 8, by = 0.25), n, TRUE), z = sample(seq(0, 8, by = 0.5), n, TRUE))
d$y <- 2 * sin(0.8 * d$x) + 0.5 * cos(d$z) + rnorm(n, 0, 0.3)
fit <- suppressWarnings(frm(frmtmb::bf(y ~ gp(x, gr = FALSE)), family = gaussian(), data = d))
fx <- function(code) gsub("cov[n, n] += 1e-12;", "cov[n, n] += 1e-6 * square(sdgp);", code, fixed = TRUE)
out <- NULL
with_reporter("silent", out <- H$brms_lp_check(brms::bf(y ~ gp(x, gr = FALSE)), gaussian(), d, fit,
                                               joint = TRUE, code_fix = fx))
cat(sprintf("GRDUP frm's line in brms's program: const %.3e grad %.2e tol %.1e dup rows %d\n",
            out$measured_const, out$max_grad, 1e-6 * abs(out$ours), sum(duplicated(d$x))))
