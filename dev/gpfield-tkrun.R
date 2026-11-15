# Lane gpfield: gf_tk() (dev/gpfield-tk.R) against the importance-sampling
# marginal (gf_is(), set.seed(9001)) at each row's fit. A tag may end in
# "@rho0/sd0": the fit from that theta start (every gp block). Args: arm
# tags(comma) ndraw
a <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-gpfield/dev/gpfield-lib.R")
gf_libs(a[1])
source("C:/Users/adf44/source/r/frmtmb-wt-gpfield/dev/gpfield-rows.R")
source("C:/Users/adf44/source/r/frmtmb-wt-gpfield/dev/gpfield-tk.R")
tags <- strsplit(a[2], ",")[[1]]
ndraw <- if (length(a) > 2) as.integer(a[3]) else 20000L
for (tg in tags) {
  p <- strsplit(tg, "@")[[1]]
  r <- gf_row(p[1])
  st <- NULL
  if (length(p) > 1) {
    rs <- as.numeric(strsplit(p[2], "/")[[1]])
    f0 <- gf_fit(r$fml, data = r$data, family = r$family)$fit
    th <- f0$estimates$theta
    for (bk in f0$frame$re_blocks) {
      if (!bk$covstruct %in% c("gp", "hsgp")) next
      th[bk$theta_idx[1]] <- log(rs[2])
      th[bk$theta_idx[-1]] <- log(rs[1])
    }
    st <- list(theta = unname(th))
  }
  q <- if (is.null(st)) gf_fit(r$fml, data = r$data, family = r$family) else
    gf_fit(r$fml, data = r$data, family = r$family, start = st)
  f <- q$fit
  if (inherits(f, "error")) {
    cat(sprintf("TK %s ERROR %s\n", tg, conditionMessage(f)))
    next
  }
  t0 <- proc.time()[[3]]
  tkp <- tryCatch(asNamespace("frmtmb")$gp_laplace_tk(f), error = function(e) NA)
  tk <- tryCatch(gf_tk(f), error = function(e) {
    cat("TKERR", conditionMessage(e), "\n"); c(tk = NA, t4 = NA, vmax = NA)
  })
  ttk <- proc.time()[[3]] - t0
  v <- if (ndraw > 0) tryCatch(gf_is(f$obj, f$opt$par, ndraw),
                               error = function(e) NULL)
  th <- f$estimates$theta
  cat(sprintf(paste0("TK %s ll %.4f conv %d nwarn %d theta %s | pkg tk %.4f | tk %.4f t4 %.4f ",
                     "vmax %.3g (%.2f s) | is %s ess %s is-lap %s ev_min %s\n"),
              tg, as.numeric(logLik(f)), f$opt$convergence, length(q$warn),
              paste(sprintf("%.3f", th), collapse = ","), as.numeric(tkp %||% NA), tk[["tk"]],
              tk[["t4"]], tk[["vmax"]], ttk,
              if (is.null(v)) "NA" else sprintf("%.4f", v[["is"]]),
              if (is.null(v)) "NA" else sprintf("%.0f", v[["ess"]]),
              if (is.null(v)) "NA" else sprintf("%.4f", v[["is"]] - v[["laplace"]]),
              if (is.null(v)) "NA" else sprintf("%.3g", v[["ev_min"]])))
  flush.console()
}
cat("DONE\n")
