# Lane priors2: which conditions two ordinal-mixture fits of the suite
# raise, on 0.82.0 and on the lane build, now that the step-start
# finding is a message (test-ordinal-mixture.R "simulate() draws ...",
# test-nonfinite-gradient.R "the three-component mixture ...").
# PRIORS2_LIB unset = lane library; PRIORS2_LIB=none = rellib-r21 only.
lib <- Sys.getenv("PRIORS2_LIB", "C:/Users/adf44/source/r/wt-priors2-lib")
if (identical(lib, "none")) lib <- ""
.libPaths(c(if (nzchar(lib)) lib, "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb", dirname(find.package("frmtmb")), "\n")
show <- function(label, expr) {
  cat("\n===", label, "\n")
  f <- withCallingHandlers(expr,
    warning = function(w) {
      cat("WARN [", paste(class(w)[1:2], collapse = ","), "]:",
          substr(gsub("\n", " ", conditionMessage(w)), 1, 150), "\n")
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      cat("MSG [", paste(class(m)[1:2], collapse = ","), "]:",
          substr(gsub("\n", " ", conditionMessage(m)), 1, 150), "\n")
      invokeRestart("muffleMessage")
    })
  cat("logLik", format(as.numeric(logLik(f)), digits = 12), "code",
      f$opt$convergence, "\n")
}
# test-ordinal-mixture.R's omx_data()
omx_data <- function(seed, n = 400) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), g = factor(sample(c("a", "b"), n, TRUE)))
  cls <- rbinom(n, 1, 0.4)
  lat <- ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) + rlogis(n)
  d$y <- as.integer(cut(lat, c(-Inf, -1.5, 0, 1.5, Inf)))
  d
}
src <- readLines("tests/testthat/test-ordinal-mixture.R")
i <- grep("^omx_data <- function", src)
if (length(i)) {
  j <- i + which(src[(i + 1):length(src)] == "}")[1L]
  eval(parse(text = src[i:j]))
}
d <- omx_data(20261011, n = 300)
show("ordinal-mixture simulate fit",
     frm(bf(y | thres(gr = g) ~ x), family = mixture(cumulative(), sratio()),
         data = d))
set.seed(20261005 + 31)
n <- 400
x <- rnorm(n)
rnorm(n)
sample(c("a", "b"), n, TRUE)
rbinom(n, 1, 0.4)
rlogis(n)
runif(n)
runif(n, 0.5, 2)
cls3 <- sample(1:3, n, TRUE, prob = c(0.3, 0.3, 0.4))
lat <- c(1.5, -0.8, 0.3)[cls3] * x + c(1.5, -1.5, 0)[cls3] + rlogis(n)
d3 <- data.frame(y = as.integer(cut(lat, c(-Inf, -1.5, 0, 1.5, Inf))),
                 x = x)
show("three-component mixture",
     frm(bf(y ~ x), family = mixture(cumulative("probit"), sratio("cloglog"),
                                     acat()), data = d3))
