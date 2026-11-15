# Lane gpfield: what the relinearization reads at a fit's optimum, on the
# lane build. Args: tags(comma)
a <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-gpfield/dev/gpfield-lib.R")
gf_libs("lane")
source("C:/Users/adf44/source/r/frmtmb-wt-gpfield/dev/gpfield-rows.R")
ns <- asNamespace("frmtmb")
for (tag in strsplit(a[1], ",")[[1]]) {
  r <- gf_row(tag)
  f <- gf_fit(r$fml, data = r$data, family = r$family)$fit
  par <- f$opt$par
  cat("DBG", tag, "ll", as.numeric(logLik(f)), "theta",
      f$estimates$theta, "\n")
  m <- tryCatch(ns$gp_mean_relin(f$obj, f$frame, par),
                error = function(e) {cat("ERR", conditionMessage(e), "\n"); NULL})
  if (!is.null(m)) cat("DBG mean relin idx", m$idx, "theta", m$theta, "\n")
  h <- tryCatch(ns$gp_hsgp_relin(f$obj, f$frame, par),
                error = function(e) {cat("ERR", conditionMessage(e), "\n"); NULL})
  if (!is.null(h)) {
    cat("DBG hsgp relin start theta", h[names(h) == "theta"], "fn",
        f$obj$fn(h), "vs", f$obj$fn(par), "\n")
    o <- nlminb(h, f$obj$fn, f$obj$gr)
    cat("DBG descent from it", -o$objective, "conv", o$convergence,
        "theta", o$par[names(o$par) == "theta"], "\n")
  }
  # objective along the length scale at the fit's other values
  th <- which(names(par) == "theta")
  for (ls in log(c(0.01, 0.03, 0.05, 0.08, 0.1, 0.2, 0.3, 0.5, 1))) {
    for (sd in log(c(0.03, 0.05, 0.1, 0.3, 1))) {
      p <- par; p[th[1]] <- sd; p[th[2]] <- ls
      v <- tryCatch(-f$obj$fn(p), error = function(e) NA)
      cat(sprintf("GRID ls %.3g sd %.3g ll %.4f\n", exp(ls), exp(sd), v))
    }
  }
}
