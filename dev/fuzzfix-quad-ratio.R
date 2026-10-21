# Print the numbers test-quadrature-nested.R asserts on, and a by-split
# fit with three blocks (depth 1 per row) under quadrature.
#   Rscript dev/fuzzfix-quad-ratio.R <lib|base>
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r10"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (a[1] == "base") c(base, ul) else c(a[1], base, ul))
suppressMessages({library(frmtmb); library(testthat)})
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
env <- testthat::test_env("frmtmb")
src <- readLines("tests/testthat/test-quadrature-nested.R")
eval(parse(text = src[1:grep("^test_that", src)[1] - 1]), env)
d <- env$quad_nest_data()
fo <- bf(y | weights(w) ~ 1 + x + (1 | ga/gb)) + bernoulli()
fq <- suppressMessages(frm(fo, data = d, quadrature = TRUE))
fl <- suppressMessages(frm(fo, data = d))
s <- vapply(fq$frame$re_blocks, function(bk) {
  th <- fq$estimates$theta[bk$theta_idx]
  sqrt(frmtmb:::covstruct_registry[[bk$covstruct]]$vcov(th, bk)[1, 1])
}, 1)
names(s) <- vapply(fq$frame$re_blocks, `[[`, "", "group_name")
ex <- env$quad_nest_exact(d, unname(fq$estimates$beta), s[["ga"]], s[["ga:gb"]])
cat(sprintf("sds ga %.6g ga:gb %.6g\n", s[["ga"]], s[["ga:gb"]]))
cat(sprintf("logLik quad %.10f exact %.10f laplace-at-quad %.10f\n",
            as.numeric(logLik(fq)), ex, -fl$obj$fn(fq$opt$par)))
eq <- abs(as.numeric(logLik(fq)) - ex); el <- abs(-fl$obj$fn(fq$opt$par) - ex)
cat(sprintf("err_q %.3g err_l %.3g ratio %.3g\n", eq, el, eq / el))
d$gby <- factor(c("a", "b", "c")[(as.integer(d$ga) %% 3) + 1])
fb <- tryCatch(suppressMessages(suppressWarnings(
  frm(bf(y ~ x + (1 | gr(ga, by = gby))) + bernoulli(), data = d,
      quadrature = TRUE))), error = function(e) e)
cat("by-split blocks:", if (inherits(fb, "error")) conditionMessage(fb)
    else sprintf("%d, left random %d, logLik %.6f",
                 length(fb$frame$re_blocks), length(fb$obj$env$random),
                 as.numeric(logLik(fb))), "\n")
