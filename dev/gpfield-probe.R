# Lane gpfield: a row's default fit against 15 theta starts (sd 0.3, 1.5, 6
# by length scale 0.1, 0.3, 1, 3, 8 in x's units), with the length scale over
# the median nearest-neighbour gap, the inner Hessian's smallest eigenvalue,
# and, if ndraw > 0, the marginal by importance sampling (gf_is()) at the
# default, at the best start and at the best "regular" start (code 0, no
# warning, every exact length scale at or above half the median gap).
# Args: arm tags(comma) ndraw
a <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-gpfield/dev/gpfield-lib.R")
gf_libs(a[1])
source("C:/Users/adf44/source/r/frmtmb-wt-gpfield/dev/gpfield-rows.R")
tags <- strsplit(a[2], ",")[[1]]
ndraw <- if (length(a) > 2) as.integer(a[3]) else 0L
desc <- function(f, x) {
  th <- f$estimates$theta
  out <- character(0)
  rat <- NA_real_
  for (bk in f$frame$re_blocks) {
    if (!bk$covstruct %in% c("gp", "hsgp")) next
    t2 <- th[bk$theta_idx]
    ls <- exp(t2[-1])
    if (bk$covstruct == "hsgp") ls <- ls * bk$gp_dmax
    r <- min(ls) / gf_mednn(x)
    rat <- min(rat, r, na.rm = TRUE)
    out <- c(out, sprintf("%s sd %.4g ls %.4g ls/mednn %.3g", bk$covstruct,
                          exp(t2[1]), min(ls), r))
  }
  list(txt = paste(out, collapse = " ; "), rat = rat)
}
for (tag in tags) {
  r <- gf_row(tag)
  t0 <- proc.time()[[3]]
  q <- gf_fit(r$fml, data = r$data, family = r$family)
  wall <- proc.time()[[3]] - t0
  f <- q$fit
  if (inherits(f, "error")) {
    cat(sprintf("PROBE %s ERROR %s\n", tag, conditionMessage(f)))
    next
  }
  dd <- desc(f, r$data$x)
  ev <- tryCatch(gf_inner_ev(f$obj, f$opt$par), error = function(e) NA)
  cat(sprintf("PROBE %s default ll %.6f conv %d nwarn %d wall %.2f inner_ev %.3g | %s | %s\n",
              tag, as.numeric(logLik(f)), f$opt$convergence, length(q$warn),
              wall, ev, dd$txt,
              paste(unique(sub(":.*", "", q$warn)), collapse = ";")))
  th0 <- f$estimates$theta
  best <- list(ll = -Inf); bestr <- list(ll = -Inf)
  for (rho0 in c(0.1, 0.3, 1, 3, 8)) for (sd0 in c(0.3, 1.5, 6)) {
    th <- th0
    for (bk in f$frame$re_blocks) {
      if (!bk$covstruct %in% c("gp", "hsgp")) next
      th[bk$theta_idx[1]] <- log(sd0)
      th[bk$theta_idx[-1]] <- log(if (bk$covstruct == "hsgp") rho0 / bk$gp_dmax
                                  else rho0)
    }
    q2 <- gf_fit(r$fml, data = r$data, family = r$family,
                 start = list(theta = unname(th)))
    if (inherits(q2$fit, "error")) next
    v <- as.numeric(logLik(q2$fit))
    if (!is.finite(v)) next
    d2 <- desc(q2$fit, r$data$x)
    reg <- q2$fit$opt$convergence == 0 && !length(q2$warn) &&
      isTRUE(d2$rat >= 0.5)
    cat(sprintf("PROBES %s rho0 %g sd0 %g ll %.6f conv %d nwarn %d reg %s | %s\n",
                tag, rho0, sd0, v, q2$fit$opt$convergence, length(q2$warn),
                reg, d2$txt))
    if (v > best$ll) best <- list(ll = v, fit = q2$fit, txt = d2$txt)
    if (reg && v > bestr$ll) bestr <- list(ll = v, fit = q2$fit, txt = d2$txt)
  }
  cat(sprintf("PROBEREF %s default %.6f best %.6f (%.4f) | best regular %.6f (%.4f) | %s\n",
              tag, as.numeric(logLik(f)), best$ll, as.numeric(logLik(f)) - best$ll,
              bestr$ll, as.numeric(logLik(f)) - bestr$ll, bestr$txt %||% ""))
  if (ndraw > 0) {
    for (wh in c("default", "best", "bestreg")) {
      ff <- switch(wh, default = f, best = best$fit, bestreg = bestr$fit)
      if (is.null(ff)) next
      v <- tryCatch(gf_is(ff$obj, ff$opt$par, ndraw), error = function(e) NULL)
      if (is.null(v)) next
      cat(sprintf(paste0("PROBEIS %s %-8s laplace %.4f is %.4f (halves %.4f ",
                         "%.4f) ess %.0f of %d ev_min %.3g\n"),
                  tag, wh, v[["laplace"]], v[["is"]], v[["h1"]], v[["h2"]],
                  v[["ess"]], ndraw, v[["ev_min"]]))
    }
  }
  flush.console()
}
cat("DONE\n")
