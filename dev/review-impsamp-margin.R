# Item 10: how near is each rewritten test to a branch rounding could decide?
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-common.R")
set_libs("after")
suppressPackageStartupMessages(library(frmtmb))
set.seed(601)
ng <- 80L
per <- 4L
imp_scalar <- local({
  src <- readLines("C:/Users/adf44/source/r/frmtmb-wt-impsamp/tests/testthat/test-importance.R")
  i <- grep("^imp_scalar_data <- function", src)
  j <- i + which(src[(i + 1L):length(src)] == "}")[1L]
  eval(parse(text = src[i:j]))
})
d <- imp_scalar()
for (R in c(1L, 5L)) {
  f <- suppressWarnings(frm(bf(y ~ x + (1 | g)) + bernoulli(), data = d, importance = 200L,
                            control = frmtmb_control(importance_rounds = R)))
  cat("rounds cap", R, "moves", format(f$importance$moves, digits = 4), "grads",
      format(f$importance$grads, digits = 4), "capped", f$importance$capped,
      "grad_tol", f$control$grad_tol, "\n")
}
# the guarded refusal: its condition on the design the test uses, with the
# wrapper moving nothing, how far below the tolerance?
f <- suppressWarnings(frm(bf(y ~ x + (1 | g)) + bernoulli(), data = d, importance = 16L))
cat("16 draws: value", format(f$opt$objective, digits = 8), "start", format(f$importance$start_value, digits = 8),
    "rise", format(f$opt$objective - f$importance$start_value, digits = 3), "tol",
    format(frmtmb:::imp_worse_tol(f$importance$mcse, f$importance$mcse), digits = 3), "\n")
# textbook model at the test's draw counts
txt <- readLines("C:/Users/adf44/source/r/frmtmb-wt-impsamp/tests/testthat/test-importance.R")
i <- grep("^imp_textbook_data <- function", txt)
j <- i + which(txt[(i + 1L):length(txt)] == "}")[1L]
imp_textbook <- eval(parse(text = txt[i:j]))
dt <- imp_textbook(124L)
dt$y <- rbinom(nrow(dt), 10L, plogis(dt$eta))
for (k in c(200L, 512L, 1000L)) {
  f <- suppressWarnings(frm(bf(y | trials(ntr) ~ x + (1 | g)), family = binomial(), data = dt,
                            importance = k))
  cat("textbook N", k, "rounds", f$importance$rounds, "moves",
      format(f$importance$moves, digits = 3), "grads", format(f$importance$grads, digits = 3), "\n")
}
# the ESS test's design, across seeds near the pinned one: does the pin rest on
# seed 5? (expect_lt(min ess, floor) and expect_gt(median ess, floor))
for (s in 1:12) {
  set.seed(s)
  g <- factor(rep(seq_len(30L), each = 2L))
  du <- data.frame(y = rbinom(60L, 1, 0.4), x = rnorm(60L), g = g)
  f <- suppressWarnings(tryCatch(frm(bf(y ~ x + (x | g)) + bernoulli(), data = du, importance = 500L),
                                 error = function(e) e))
  if (inherits(f, "error")) {
    cat("seed", s, "ERROR", substr(conditionMessage(f), 1, 100), "\n")
  } else {
    cat("seed", s, "min ess", format(min(f$importance$ess), digits = 3), "median",
        format(median(f$importance$ess), digits = 3), "\n")
  }
}
