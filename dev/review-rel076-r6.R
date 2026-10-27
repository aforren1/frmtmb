# Review, short check of R6 (ac_bound_theta()'s snap at a computed end):
# the snap is autocor_bound_chart()'s; a bound genuinely past an end is
# still refused, with the ends printed to 15 digits; ub = 1 works under
# frm() and frm_sample() for both sar() types; an lb a few ulps inside
# 1 / min eigenvalue is no bound, an lb well inside is a bound.
#   Rscript dev/review-rel076-r6.R > dev/review-rel076-out/r6.txt
.libPaths(c("C:/Users/adf44/source/r/rellib-r14",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages({
  library(frmtmb); library(frmtmb.sample); library(frmtmb.spatial)
})
ns <- asNamespace("frmtmb")
src <- function(f) paste(deparse(get(f, ns)), collapse = "\n")
cat("snap in autocor_bound_chart():",
    grepl("64 * .Machine$double.eps * max(1, abs(b))",
          src("autocor_bound_chart"), fixed = TRUE),
    "| in ac_bound_theta():",
    grepl("64 * .Machine$double.eps * max(1, abs(e))",
          src("ac_bound_theta"), fixed = TRUE), "\n")
n <- 20
A <- matrix(0, n, n)
for (i in 1:n) {
  j <- i %% n + 1; A[i, j] <- A[j, i] <- 1
  j <- (i + 1) %% n + 1; A[i, j] <- A[j, i] <- 1
}
A[1, 10] <- A[10, 1] <- 1
W <- A / rowSums(A)
set.seed(3)
d <- data.frame(x = stats::rnorm(n))
d$y <- as.vector(solve(diag(n) - 0.4 * W, 1 + 0.5 * d$x +
                         stats::rnorm(n, 0, 0.7)))
d2 <- list(W = W)
try1 <- function(lab, expr) {
  r <- tryCatch({v <- expr; if (inherits(v, "frmtmb_fit"))
    paste("FIT, logLik", format(as.numeric(logLik(v)), digits = 12)) else
      "RAN"}, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(sprintf("  %-44s %s\n", lab, r))
}
eps <- .Machine$double.eps
for (ty in c("lag", "error")) {
  fo <- stats::as.formula(paste0("y ~ x + sar(W, type = '", ty, "')"))
  cl <- paste0(ty, "sar")
  dp <- frmtmb::default_prior(fo, data = d, data2 = d2)
  lo <- dp$lb[dp$class == cl]; hi <- dp$ub[dp$class == cl]
  cat("\n==", cl, ": table ends", format(lo, digits = 17),
      format(hi, digits = 17), "\n")
  f0 <- frm(fo, data = d, data2 = d2)
  try1("no bound", f0)
  sp <- function(lb = NA, ub = NA) set_prior("", class = cl, lb = lb, ub = ub)
  try1("ub = 1", frm(fo, data = d, data2 = d2, prior = sp(ub = 1)))
  try1("lb = table lb, ub = 1", frm(fo, data = d, data2 = d2,
                                     prior = sp(lb = lo, ub = 1)))
  try1("ub = 1 + 1e-9", frm(fo, data = d, data2 = d2,
                             prior = sp(ub = 1 + 1e-9)))
  try1("ub = hi + 65 ulps", frm(fo, data = d, data2 = d2,
                                 prior = sp(ub = hi + 65 * eps)))
  try1("ub = hi + 60 ulps (snapped)", frm(fo, data = d, data2 = d2,
                                           prior = sp(ub = hi + 60 * eps)))
  try1("lb = lo - 1e-9", frm(fo, data = d, data2 = d2,
                              prior = sp(lb = lo - 1e-9)))
  try1("lb = lo + 8 ulps (inside, snapped)",
       frm(fo, data = d, data2 = d2, prior = sp(lb = lo + 8 * eps * abs(lo))))
  try1("lb = lo - 8 ulps (outside, snapped)",
       frm(fo, data = d, data2 = d2, prior = sp(lb = lo - 8 * eps * abs(lo))))
  # a real bound well inside binds: lb = 0.6 above the ML estimate
  try1("lb = 0.6 (a real bound)", frm(fo, data = d, data2 = d2,
                                       prior = sp(lb = 0.6)))
  ds <- tryCatch(suppressWarnings(suppressMessages(frm_sample(
    fo, data = d, data2 = d2, prior = sp(lb = 0, ub = 1), chains = 1,
    iter = 300, seed = 2, refresh = 0))), error = function(e) e)
  cat("  frm_sample(lb = 0, ub = 1):", if (inherits(ds, "error"))
    paste("ERROR:", conditionMessage(ds)) else {
      m <- as.matrix(ds$draws)
      paste("draws from", format(min(m[, cl])), "to", format(max(m[, cl])))
    }, "\n")
  ds2 <- tryCatch(suppressWarnings(suppressMessages(frm_sample(
    fo, data = d, data2 = d2, prior = sp(ub = 1 + 1e-9), chains = 1,
    iter = 300, seed = 2, refresh = 0))), error = function(e) e)
  cat("  frm_sample(ub = 1 + 1e-9):", if (inherits(ds2, "error"))
    paste("ERROR:", conditionMessage(ds2)) else "RAN", "\n")
}
# the other maps: ar's ends -1 and 1, a hair past them
d3 <- expand.grid(t = 1:5, g = 1:10); set.seed(1); d3$y <- rnorm(50)
try1("ar cov: ub = 1 + 1e-12", frm(y ~ ar(t, g, cov = TRUE), data = d3,
     prior = set_prior("", class = "ar", ub = 1 + 1e-12)))
try1("ar cov: lb = -1 - 8 ulps", frm(y ~ ar(t, g, cov = TRUE), data = d3,
     prior = set_prior("", class = "ar", lb = -1 - 8 * eps)))
