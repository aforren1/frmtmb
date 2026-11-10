# Review of lane code8, item 2 (cost): the fit, the outer Hessian and the
# SE analysis on large random-effect fits, before and after.
#   RV_ARM=before|after Rscript dev/review-code8-cost.R <label>
lab <- commandArgs(TRUE)[1]
source("dev/review-code8-common.R")
mk <- function(n, ng, seed) {
  set.seed(seed)
  d <- data.frame(g = factor(sample(seq_len(ng), n, TRUE)),
                  x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n, 5, 2))
  u0 <- rnorm(ng, 0, 0.5); u1 <- rnorm(ng, 0, 0.3)
  d$y <- 1 + (0.5 + u1[d$g]) * d$x1 - 0.3 * d$x2 + 0.1 * d$x3 +
    u0[d$g] + rnorm(n)
  d$yp <- rpois(n, exp(0.2 + 0.3 * d$x1 + u0[d$g]))
  d$yn <- 2 + 3 * exp(-0.4 * (d$x3 + u0[d$g])) + rnorm(n, 0, 0.5)
  d
}
cases <- list(
  list("lmm 1e4", 1e4, function(d) frm(y ~ x1 + x2 + x3 + (1 + x1 | g), data = d)),
  list("lmm 1e5", 1e5, function(d) frm(y ~ x1 + x2 + x3 + (1 + x1 | g), data = d)),
  list("glmm 1e5", 1e5, function(d) frm(yp ~ x1 + x2 + (1 | g), data = d,
                                          family = poisson())),
  list("nl 2e4", 2e4, function(d) frm(bf(yn ~ c0 + A * exp(-k * x3),
                                         c0 ~ 1, A ~ 1 + (1 | g), k ~ 1,
                                         nl = TRUE), data = d,
                                      prior = c(prior(normal(2, 5), nlpar = "c0"),
                                                prior(normal(3, 5), nlpar = "A"),
                                                prior(normal(0.4, 1), nlpar = "k")))),
  list("lmm 1e5 p12", 1e5, function(d) {
    for (j in 1:8) d[[paste0("z", j)]] <- rnorm(nrow(d))
    frm(y ~ x1 + x2 + x3 + z1 + z2 + z3 + z4 + z5 + z6 + z7 + z8 + (1 | g),
        data = d)
  })
)
out <- list()
for (cs in cases) {
  d <- mk(cs[[2]], 200, 11)
  ng <- 0L
  t0 <- proc.time()[["elapsed"]]
  r <- withCallingHandlers(cs[[3]](d), warning = function(w) {
    cat("  W:", substr(conditionMessage(w), 1, 100), "\n")
    invokeRestart("muffleWarning")
  }, message = function(m) invokeRestart("muffleMessage"))
  tfit <- proc.time()[["elapsed"]] - t0
  f <- r
  p <- f$opt$par
  # the Hessian as the SE check builds it, counted
  gr0 <- f$obj$gr
  cnt <- 0L
  gcount <- function(x) { cnt <<- cnt + 1L; gr0(x) }
  eps <- ns$outer_fd_steps(f$frame, names(p))
  t1 <- proc.time()[["elapsed"]]
  h <- ns$fd_hessian(p, gcount, eps = eps)
  th <- proc.time()[["elapsed"]] - t1
  d0 <- diag(h$H)
  nredo <- sum(0.5 * d0 * eps^2 > 0.01)
  # the analysis from scratch
  f$obj$env$.frmtmb_outer_h <- NULL
  for (nm in c("se_fit_an", "hessian_fixed", "se_analysis", "sdr")) {
    if (exists(nm, f$cache)) rm(list = nm, envir = f$cache)
  }
  t2 <- proc.time()[["elapsed"]]
  an <- quiet(ns$se_fit_analysis0(f, TRUE))
  ta <- proc.time()[["elapsed"]] - t2
  s <- se_of(f)
  out[[length(out) + 1L]] <- data.frame(case = cs[[1]], np = length(p),
    t_fit = tfit, t_hess = th, grads = cnt, redo = nredo, t_an = ta,
    tier = an$tier, lost = length(s$lost))
  cat(cs[[1]], "SE:", paste(sprintf("%s=%.10g", names(s$se), s$se),
                             collapse = " "), "\n")
  print(out[[length(out)]])
}
X <- do.call(rbind, out)
print(X)
write.table(X, sprintf("dev/review-code8-out/cost-%s.tsv", lab), sep = "\t",
            quote = FALSE, row.names = FALSE)
cat("DONE\n")
