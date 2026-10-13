# Reviewer of lane seflags, claim 1: the "0.5 * largest projection" cap.
# Where the noise bound tau is meaningless, flat_tol is the cap, and a
# prediction with a cosine below it against an EXACTLY flat direction is
# reported finite. (1) ported fixture 1 (an exact alias); (2) the lane's
# alias construction plus a near-collinear pair that shrinks the gap.
#   Rscript dev/seflags-rev-pred2.R <lane|base>
source("C:/Users/adf44/source/r/frmtmb-wt-seflags/dev/seflags-rev-common.R")
source("C:/Users/adf44/source/r/frmtmb-wt-seflags/dev/seflags-rev-fns.R")
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true")
suppressMessages(library(testthat))

alias_null <- function(fit, d, xname) {
  X <- fit$frame$linpreds[[1]]$X
  cn <- colnames(X)
  ifx <- grep(paste0("^s[(]", xname, "[)][.]?fx"), cn)
  cf <- coef(stats::lm(as.numeric(X[, ifx[1]]) ~ d[[xname]]))
  on <- ns$outer_par_names(fit)
  nex <- numeric(length(on))
  nex[on == "(Intercept)"] <- cf[1]
  nex[on == xname] <- cf[2]
  nex[on == paste0("s(", xname, ").fx1")] <- -1
  cat("  exact alias: fx1 =", format(cf, digits = 6), "\n")
  nex
}

## (1) fixture 1
env <- new.env(parent = ns)
root <- "C:/Users/adf44/source/r/frmtmb-wt-seflags/tests/testthat"
for (h in list.files(root, "^helper-.*[.]R$", full.names = TRUE)) {
  sys.source(h, envir = env)
}
withr::local_seed(3L)
fit <- suppressWarnings(suppressMessages(env$brms_fixture(1)))
cat("\n== fixture 1 | code", fit$opt$convergence, "\n")
d1 <- fit$data
nex <- alias_null(fit, d1, "Age")
cos_sweep(fit, nex, "fixture 1")
nd <- d1[1:5, ]
p <- cap(fitted(fit, newdata = nd, re_formula = NA))
cat("  fitted(newdata, re_formula = NA): warnings", length(p$w), "se",
    paste(signif(p$v[, 2], 4), collapse = " "), "\n")

## (2) alias + near-collinear pair + (1 | g)
for (s in c(1e-3, 1e-4, 1e-5)) {
  set.seed(1)
  n <- 200
  x <- stats::rnorm(n)
  x <- x - mean(x) + 0.1
  z1 <- stats::rnorm(n)
  d <- data.frame(x = x, z1 = z1, z2 = z1 + stats::rnorm(n, 0, s),
                  g = factor(rep(1:10, n / 10)))
  d$y <- 1 + 0.5 * d$x + sin(2 * d$x) + 0.3 * d$z1 + 0.3 * d$z2 +
    stats::rnorm(10, 0, 0.5)[d$g] + stats::rnorm(n, 0, 0.3)
  r <- cap(frm(y ~ x + s(x) + z1 + z2 + (1 | g), data = d))
  fit <- r$v
  cat("\n== alias + z2 = z1 + N(0,", s, ") | code", fit$opt$convergence,
      "| fit warnings", length(r$w), "\n")
  nex <- alias_null(fit, d, "x")
  cos_sweep(fit, nex, "alias+collinear")
  # the identified pair against lm() on the same fixed design with g fixed
  se <- se_outer(fit)
  cat("  SE z1, z2:", format(se[c("z1", "z2")], digits = 5), "\n")
}
