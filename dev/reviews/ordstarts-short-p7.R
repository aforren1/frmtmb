# Reviewer short check of lane ordstarts (build fix round 1), item 4: the
# weights' reference at K = 2 and K = 3. For each swapped start, component
# k must get component sigma(k)'s mixing weight: max |pi'_k - pi_sigma(k)|
# over components and assignments (intercept-only weights).
#   Rscript ordstarts-short-p7.R > ordstarts-short-p7.txt
.libPaths(c("C:/Users/adf44/source/r/wt-ordstarts-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
q <- function(e) suppressWarnings(suppressMessages(e))
pis <- function(par, nm, K, ref) {
  L <- numeric(K)
  for (k in setdiff(seq_len(K), ref)) {
    L[k] <- par[nm == paste0("theta", k, "_(Intercept)")]
  }
  exp(L) / sum(exp(L))
}
check <- function(label, f, K) {
  nm <- ns$outer_par_map(f)$names
  ref <- f$spec$responses$y$family$mix$ref %||% K
  sw <- ns$mixture_ord_swap_starts(f, "y")
  p0 <- pis(f$opt$par, nm, K, ref)
  err <- vapply(sw, function(s) {
    p1 <- pis(s$par, nm, K, ref)
    max(abs(p1 - p0[s$sigma]))
  }, 0)
  cat(sprintf("%-34s ref %d | weights %s = %s | starts %d | max |pi'_k - pi_sigma(k)| %.2g\n",
              label, ref, paste(grep("^theta", nm, value = TRUE), collapse = ","),
              paste(round(f$opt$par[grep("^theta", nm)], 4), collapse = ","),
              length(sw), max(err)))
}
mk <- function(seed, n) {
  set.seed(seed)
  d <- data.frame(x = rnorm(n), z = rnorm(n))
  cls <- rbinom(n, 1, 0.4)
  lat <- ifelse(cls == 1, 1.5 * d$x + 1, -0.8 * d$x - 1) + rlogis(n)
  d$y <- 1L + (lat > -1 + 0.1 * d$z) + (lat > 0) + (lat > 1 - 0.1 * d$z)
  d
}
d2 <- mk(13, 1000)
for (form in list(bf(y ~ x + z), bf(y ~ x + z, theta1 ~ 1),
                  bf(y ~ x + z, theta2 ~ 1))) {
  f <- q(frm(form, family = mixture(cumulative(), sratio()), data = d2,
             control = frmtmb_control(starts = 0)))
  check(paste("K=2", paste(names(form$pforms), collapse = ",")), f, 2)
}
set.seed(3)
n <- 900
d <- data.frame(x = rnorm(n), z = rnorm(n))
cls <- sample(1:3, n, TRUE)
lat <- c(1.5, -0.8, 0)[cls] * d$x + c(1, -1, 0)[cls] + rlogis(n)
d$y <- 1L + (lat > -1) + (lat > 0) + (lat > 1)
for (form in list(bf(y ~ x), bf(y ~ x, theta1 ~ 1, theta2 ~ 1),
                  bf(y ~ x, theta2 ~ 1, theta3 ~ 1),
                  bf(y ~ x, theta1 ~ 1, theta3 ~ 1))) {
  f <- q(frm(form, family = mixture(cumulative(), sratio(), acat()),
             data = d, control = frmtmb_control(starts = 0)))
  check(paste("K=3", paste(names(form$pforms), collapse = ",")), f, 3)
}
cat("DONE p7\n")
