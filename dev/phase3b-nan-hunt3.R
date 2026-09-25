# NaN hunt, contaminant density: rows above and below the non-decision
# time, the unreachable variant a stated max_ndt selects, wild values.
.libPaths(c("C:/Users/adf44/source/r/phase3b-lib2",
            "C:/Users/adf44/source/r/phase3b-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(pkgload::load_all("extensions/frmtmb.eam", quiet = TRUE,
                                   export_all = TRUE))
y <- c(0.12, 0.2, 0.26, 0.3, 0.5, 1, 3, 4.9)
up <- c(0, 1, 0, 1, 1, 0, 1, 0)
cr <- c(0.12, 4.9)
delta <- 1e-9 * 0.12
base_unreach <- function(y, dp, at) {
  ddm_lpdf_both(ddm_floor(y - dp$ndt - delta, delta), dp$mu, dp$bs,
                dp$bias, at$dec)
}
base_plain <- function(y, dp, at) {
  ddm_lpdf_both(y - dp$ndt, dp$mu, dp$bs, dp$bias, at$dec)
}
for (nm in c("unreach", "plain")) {
  lp <- ddm_cont_lpdf(if (nm == "unreach") base_unreach else base_plain, cr)
  f <- function(p) {
    dp <- list(mu = p[1], bs = exp(p[2]), ndt = p[3], bias = 0.5,
               lambda = 1 / (1 + exp(-p[4])), .eta_lambda = p[4])
    sum(lp(y, dp, list(dec = up)))
  }
  tp <- RTMB::MakeTape(f, c(0, 0, 0.1, -3))
  grid <- expand.grid(v = c(-30, -3, 0, 3, 30), la = c(-5, -2, 0, 1, 3, 6),
                      t0 = if (nm == "unreach") c(0.01, 0.2, 0.29, 0.45)
                           else c(0.01, 0.1, 0.119),
                      el = c(-40, -20, -3, 0, 3))
  bad <- 0; ex <- NULL; big <- 0
  for (i in seq_len(nrow(grid))) {
    p <- as.numeric(grid[i, ])
    g <- tp$jacobian(p)
    if (!all(is.finite(g)) || !is.finite(tp(p))) {
      bad <- bad + 1
      if (is.null(ex)) ex <- c(p, tp(p), g)
    }
    big <- max(big, abs(g[is.finite(g)]))
  }
  cat(nm, ": non-finite at", bad, "of", nrow(grid), "; largest |grad|",
      format(big, digits = 3), "\n")
  if (!is.null(ex)) cat("  first:", format(ex, digits = 4), "\n")
}
