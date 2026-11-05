# Item 3, part 1: global state after a fit; TapeConfig restored; compress flag;
# factor failure on the tape; tied rows; huge length scale; duplicates; by, iso,
# multi-d; numeric factor against tape factor. Seeds below.
source("C:/Users/adf44/source/r/frmtmb-wt-ncgp/dev/review-ncgp-lib.R")
rev_libs("lane")
ns <- asNamespace("frmtmb")
cat("== global state\n")
library(RTMB)
tc0 <- RTMB::TapeConfig()
print(unlist(tc0))
cmp0 <- TMB::config(DLL = "RTMB")$tmbad.sparse_hessian_compress
set.seed(1)
d <- data.frame(x = runif(30, 0, 5)); d$y <- sin(d$x) + rnorm(30, sd = 0.2)
# node count of a non-gp fit before and after a gp fit in one process
dg <- data.frame(y = rnorm(200), z = rnorm(200), g = factor(rep(1:20, 10)))
nodes <- function(f) {
  c(F = nrow(RTMB::GetTape(f$obj, "ADFun")$data.frame()),
    H = nrow(RTMB::GetTape(f$obj, "ADHess")$data.frame()))
}
f0 <- frm(y ~ z + (1 | g), dg)
n0 <- try(nodes(f0))
fg <- frm(y ~ gp(x), d)
print(unlist(RTMB::TapeConfig()))
cat("compress before", cmp0, "after", TMB::config(DLL = "RTMB")$tmbad.sparse_hessian_compress, "\n")
cat("TapeConfig restored:", identical(unlist(tc0), unlist(RTMB::TapeConfig())), "\n")
f1 <- frm(y ~ z + (1 | g), dg)
n1 <- try(nodes(f1))
print(rbind(n0, n1))
cat("non-gp ll identical:", identical(logLik(f0), logLik(f1)), "\n")
# fit error inside the scope: is state restored?
e <- try(frm(y ~ gp(x), d, start = list(theta = c(1, 2, 3))), silent = TRUE)
cat("after an error in a fit: ", paste(unlist(RTMB::TapeConfig()), collapse = ","), " compress ",
    TMB::config(DLL = "RTMB")$tmbad.sparse_hessian_compress, "\n")

cat("== factor failure on the tape\n")
bk <- Filter(function(k) identical(k$covstruct, "gp"), fg$frame$re_blocks)[[1]]
ob <- fg$obj
p <- fg$opt$par
for (lrho in log(c(1e2, 1e4, 1e6, 1e7, 1e8, 1e9, 1e10, 1e12, 1e15, 1e30))) {
  q <- p; q[which(names(q) == "theta")[2]] <- lrho
  v <- try(ob$fn(q), silent = TRUE)
  g <- try(ob$gr(q), silent = TRUE)
  cat(sprintf("rho %.0e fn %s gr %s\n", exp(lrho), format(v), paste(format(g, digits = 3), collapse = " ")))
}
v <- ob$fn(p); cat("back at the optimum, delta fn:", v - fg$opt$objective, "\n")
for (lrho in log(c(1e-10, 1e-30, 1e-100, 1e-300))) {
  q <- p; q[which(names(q) == "theta")[2]] <- lrho
  v <- try(ob$fn(q), silent = TRUE)
  cat(sprintf("small rho %.0e fn %s\n", exp(lrho), format(v)))
}
for (lsd in c(-800, -400, -100, 100, 400, 800)) {
  q <- p; q[which(names(q) == "theta")[1]] <- lsd
  v <- try(ob$fn(q), silent = TRUE); g <- try(ob$gr(q), silent = TRUE)
  cat(sprintf("log sd %g fn %s gr %s\n", lsd, format(v), paste(format(g, digits = 3), collapse = " ")))
}
cat("== tape built at a failing or huge start\n")
for (lrho in log(c(1e6, 1e10, 1e15, 1e30))) {
  r <- rev_fit(y ~ gp(x), d, start = list(theta = c(0, lrho)))
  if (inherits(r$fit, "error")) {
    cat(sprintf("start rho %.0e ERROR %s\n", exp(lrho), substr(conditionMessage(r$fit), 1, 160)))
  } else {
    g <- rev_gp_read(r$fit)
    cat(sprintf("start rho %.0e ll %.5f rho %.4g sd %.4g conv %d nwarn %d\n", exp(lrho), g$ll, g$rho, g$sd, g$conv, length(r$warn)))
  }
}
cat("== numeric factor against tape factor: field and objective, hard designs\n")
for (cs in list(list(n = 60, s = 1e-3, x = "unif"), list(n = 200, s = 1e-4, x = "unif"),
                list(n = 200, s = 1e-2, x = "grid"), list(n = 300, s = 1, x = "tight"))) {
  set.seed(3)
  x <- switch(cs$x, unif = sort(runif(cs$n, 0, 10)), grid = seq(0, 10, length.out = cs$n),
              tight = sort(runif(cs$n, 0, 0.01)))
  y <- sin(2 * x) + cs$s * rnorm(cs$n)
  r <- rev_fit(y ~ gp(x), data.frame(y, x))
  f <- r$fit
  if (inherits(f, "error")) { cat("ERROR", conditionMessage(f), "\n"); next }
  e <- f$obj$env
  # the tape's joint at the mode, against the objective on numbers
  nll <- ns$build_objective(f$frame)
  num <- nll(f$estimates)
  tp <- e$f(e$last.par.best)
  bk <- Filter(function(k) identical(k$covstruct, "gp"), f$frame$re_blocks)[[1]]
  th <- f$estimates$theta[bk$theta_idx]
  Lnum <- ns$gp_factor_ordered(th, bk)
  rtmbL <- local({
    K <- ns$gp_corr(th, bk)[bk$gp_order, bk$gp_order]
    tape <- RTMB::MakeTape(function(p) { M <- RTMB::matrix(p, nrow(K), nrow(K)); t(chol(M)) }, as.vector(K))
    matrix(tape(as.vector(K)), nrow(K), nrow(K))
  })
  cat(sprintf("n %d s %g %s: joint tape %.10f numbers %.10f diff %.2e | L max abs diff %.2e | field diff %.2e | conv %d nwarn %d cond %.1e\n",
              cs$n, cs$s, cs$x, tp, num, tp - num, max(abs(Lnum - rtmbL)),
              max(abs((Lnum - rtmbL) %*% f$estimates$b[bk$b_idx])) * exp(th[1]),
              f$opt$convergence, length(r$warn), kappa(ns$gp_corr(th, bk), exact = TRUE)))
}
cat("== duplicates, ties, by, iso, multi-d at low noise\n")
set.seed(5)
u <- sort(runif(25, 0, 6))
xr <- rep(u, each = 4)
dd <- data.frame(x = xr, z = rep(runif(25, 0, 6), each = 4), f = factor(rep(c("a", "b"), 50)),
                 w = rep(c(1, 2), 50))
dd$y <- sin(dd$x) + cos(dd$z) * (dd$f == "a") + 1e-3 * rnorm(100)
for (spec in c("y ~ gp(x, gr = FALSE)", "y ~ gp(x)", "y ~ gp(x, by = f)", "y ~ gp(x, by = f, cmc = FALSE)",
               "y ~ gp(x, by = w)", "y ~ gp(x, z)", "y ~ gp(x, z, iso = FALSE)", "y ~ gp(x, z, gr = FALSE)",
               "y ~ gp(x) + gp(z)", "y ~ gp(x, scale = FALSE)")) {
  t0 <- Sys.time()
  r <- rev_fit(as.formula(spec), dd)
  if (inherits(r$fit, "error")) {
    cat(sprintf("%-32s ERROR %s\n", spec, substr(conditionMessage(r$fit), 1, 200)))
    next
  }
  f <- r$fit
  ft <- try(fitted(f), silent = TRUE)
  cat(sprintf("%-32s ll %.4f conv %d nwarn %d sdgp %s lscale %s fittedOK %s maxse %.3g sec %.1f\n", spec, as.numeric(logLik(f)),
              f$opt$convergence, length(r$warn),
              paste(signif(exp(f$estimates$theta[vapply(f$frame$re_blocks, function(k) k$theta_idx[1], 1L)]), 4), collapse = ","),
              paste(signif(exp(f$estimates$theta[unlist(lapply(f$frame$re_blocks, function(k) k$theta_idx[-1]))]), 4), collapse = ","),
              !inherits(ft, "try-error"), if (inherits(ft, "try-error")) NA else max(ft[, 2]),
              as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  if (length(r$warn)) cat("   W:", substr(r$warn[1], 1, 140), "\n")
  # gp_order is a permutation of each block
  for (bk in Filter(function(k) identical(k$covstruct, "gp"), f$frame$re_blocks)) {
    if (!identical(sort(bk$gp_order), seq_len(length(bk$b_idx)))) cat("   gp_order NOT a permutation\n")
  }
}
cat("DONE\n")
