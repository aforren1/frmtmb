# Reviewer: which censoring code drives the AD / FD mismatch of the
# skew_normal (alpha = -10) and hurdle_gamma (shape e^4) objectives in
# cdfs-rev-grad.R: the same data with one code at a time switched off.
.libPaths(c("C:/Users/adf44/source/r/wt-cdfs-rev-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
options(frmtmb.notices = FALSE)
src <- readLines("C:/Users/adf44/source/r/frmtmb-wt-cdfs/dev/cdfs-rev-grad.R")
cut <- grep("^pick <- ", src)
for (cs in list(list("skew_normal", 2, "alpha", -10), list("hurdle_gamma", 1, "shape", 4),
                list("von_mises", 3, "kappa", 4))) {
  Sys.setenv(GROUP = cs[[2]])
  eval(parse(text = src[1:(cut - 1)]))
  nm <- cs[[1]]
  fm <- fams[[nm]]
  dpn <- setdiff(frmtmb:::as_frmtmb_family(fm$fam)$dpars, "mu")
  fo <- do.call(bf, c(list(y | cens(cc, y2) ~ x),
                      lapply(dpn, function(dp) as.formula(paste(dp, "~ 1")))))
  for (keep in c("all", "-1", "1", "2")) {
    d <- fm$d
    if (keep != "all") d$cc[d$cc != as.numeric(keep)] <- 0
    ob <- frm(fo, family = fm$fam, data = d, dry_run = "objective")$obj
    p <- ob$par
    jd <- which(names(p) == "betad")
    p[jd[match(cs[[3]], dpn)]] <- cs[[4]]
    g <- as.numeric(ob$gr(p))
    fd <- numDeriv::grad(ob$fn, p, method = "Richardson")
    cat(sprintf("%-14s codes kept %-4s (n %3d): max rel |AD - FD| %.2e  fn %.8g\n", nm, keep,
                sum(d$cc != 0), max(abs(g - fd) / pmax(1, abs(fd))), ob$fn(p)))
  }
  # the interval rows: F(y2) - F(y) on the probability scale vs the
  # same mass from the log pairs
  d <- fm$d
  ii <- which(d$cc == 2)
  f <- frmtmb:::as_frmtmb_family(fm$fam)
  ob <- frm(fo, family = fm$fam, data = d, dry_run = "objective")$obj
  p <- ob$par
  jd <- which(names(p) == "betad")
  p[jd[match(cs[[3]], dpn)]] <- cs[[4]]
  mu_eta <- p[1] + p[2] * d$x
  dp <- list(mu = f$links$mu$linkinv(mu_eta))
  for (j in seq_along(dpn)) dp[[dpn[j]]] <- f$links[[dpn[j]]]$linkinv(p[jd[j]])
  F1 <- f$lcdf(d$y[ii], lapply(dp, function(v) if (length(v) > 1) v[ii] else v), list())
  F2 <- f$lcdf(d$y2[ii], lapply(dp, function(v) if (length(v) > 1) v[ii] else v), list())
  S1 <- f$lccdf(d$y[ii], lapply(dp, function(v) if (length(v) > 1) v[ii] else v), list())
  S2 <- f$lccdf(d$y2[ii], lapply(dp, function(v) if (length(v) > 1) v[ii] else v), list())
  viaF <- log(F2 - F1)
  viaS <- S1 + log(-expm1(S2 - S1))
  k <- which.max(abs(viaF - viaS))
  cat(sprintf("  interval rows: worst |log(F2 - F1) - log(S1 - S2)| %.3g at y %.4g y2 %.4g (F1 %.17g, viaF %.6g viaS %.6g)\n",
              abs(viaF - viaS)[k], d$y[ii][k], d$y2[ii][k], F1[k], viaF[k], viaS[k]))
}
