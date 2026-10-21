# Reviewer: which fuzz random-effect structures give every row its own
# effect (so the weights refusal can fire), and what the 0.6 scaling of
# fuzz_design() does to sum(w) / n.
#   Rscript dev/fuzzfix-rev-fuzzw.R <lib>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
env <- new.env()
sys.source("tests/testthat/helper-fuzz.R", envir = env)

# 1. sum(w) / n under the generator, by design size, 20000 seeds each
for (rn in 0:3) {
  seeds <- rn + 4L * (0:19999)
  r <- vapply(seeds, function(s) {
    set.seed(s)
    sp <- list(seed = s, re = "ar1")
    d <- env$fuzz_design(sp)
    sum(d$w) / nrow(d)
  }, 1)
  set.seed(rn)
  n <- nrow(env$fuzz_design(list(seed = rn, re = "ar1")))
  cat(sprintf("design %d n %d: ar1 sum(w)/n min %.4f max %.4f (bound 1.2)\n",
              rn, n, min(r), max(r)))
}

# 2. ownership: does any re structure besides ar1 give every row its own
# effect? With weights of 1.5 the refusal fires exactly when it does.
for (rn in 0:3) {
  for (re in names(env$fuzz_re)) {
    set.seed(rn)
    d <- env$fuzz_design(list(seed = rn, re = "none"))
    d$w15 <- 1.5
    d$y <- rnorm(nrow(d))
    term <- env$fuzz_re[[re]]$term
    fo <- stats::as.formula(paste("y | weights(w15) ~ x",
                                  if (!is.null(term)) paste("+", term)))
    r <- tryCatch(suppressWarnings(suppressMessages(
      frm(bf(fo), family = gaussian(), data = d,
          control = frmtmb_control(check_se = "ignore")))),
      error = function(e) e)
    cat(sprintf("design %d re %-8s %s\n", rn, re,
                if (inherits(r, "error") &&
                    grepl("has no maximum", conditionMessage(r))) "REFUSED"
                else if (inherits(r, "error")) paste("error:",
                  substr(conditionMessage(r), 1, 60))
                else "fits"))
  }
}
