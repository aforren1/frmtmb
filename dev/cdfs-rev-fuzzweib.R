# Reviewer: the weibull + trunc + diagbar fuzz spec (seed 20446531).
# Is its NaN gradient the truncation normalizer log(1 - F(lb)) on the
# probability scale? Fit it as the fuzz harness does, then again with
# row_lpdf() patched in this session to divide by exp(lccdf(lb)) when
# there is no upper bound (the log survivor every CDF family now has).
.libPaths(c("C:/Users/adf44/source/r/wt-cdfs-rev-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(testthat); library(frmtmb)})
options(frmtmb.notices = FALSE)
wt <- "C:/Users/adf44/source/r/frmtmb-wt-cdfs"
testthat::source_test_helpers(file.path(wt, "tests/testthat"), env = globalenv())
plan <- fuzz_plan(seed = 20260901L, size = 300L)
run_spec <- function(seed, lab) {
  sp <- as.list(plan[plan$seed == seed, , drop = FALSE][1, ])
  d <- fuzz_data(sp)
  txt <- fuzz_call_text(sp)
  r <- tryCatch({
    fit <- suppressWarnings(eval(parse(text = txt)))
    sprintf("fits: logLik %.6f conv %d", as.numeric(logLik(fit)), fit$opt$convergence)
  }, error = function(e) paste("ERROR:", substr(conditionMessage(e), 1, 90)))
  cat(sprintf("%-10s seed %d %s\n   %s\n", lab, seed, fuzz_spec_key(sp), r))
  invisible(list(sp = sp, d = d, txt = txt))
}
x <- run_spec(20446531, "as is")
cat("   call:", gsub("\n", " ", x$txt), "\n")
cat("   trunc_lb", fuzz_families$weibull$trunc_lb, " min y", min(x$d$y), "\n")
ns <- asNamespace("frmtmb")
f <- ns$row_lpdf
b <- paste(deparse(f, width.cutoff = 500L), collapse = "\n")
old <- "ll <- ll - log(Fub - Flb)"
stopifnot(grepl(old, b, fixed = TRUE))
b <- sub(old, paste0("ll <- ll - (if (is.null(av[[\"trunc_ub\"]]) && !is.null(lb) && ",
                     "has_lccdf(fam)) fam_lccdf(fam, lb, dpv, av, extra) else ",
                     "log(Fub - Flb))"), b, fixed = TRUE)
g <- eval(parse(text = b))
environment(g) <- ns
utils::assignInNamespace("row_lpdf", g, "frmtmb")
run_spec(20446531, "patched")
cat("== the other two specs, as is (no trunc in either)\n")
run_spec(20502220, "as is")
