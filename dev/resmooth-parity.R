# Lane wt-resmooth. Design-space parity with brms for a smooth indexed
# by a grouping factor, at re_formula = NA.
#
# brms's coefficients cannot be pushed into a frmtmb fit (different
# basis rotation), so the comparison is the one the task allows instead:
# does frmtmb's re_formula = NA design SPAN brms's re_formula = NA
# answer? The residual of least squares of brms's posterior-mean epred
# on frmtmb's design is reported relative to the norm of that epred
# after centring, once against the design this change produces and once
# against the design 0.64.0 produced (the group-indexed smooth blocks
# removed). A span that holds the answer gives a residual at the level
# of the linear algebra; one that does not cannot.
#   Rscript dev/resmooth-parity.R > dev/resmooth-parity.txt
.libPaths(c("C:/Users/adf44/source/r/wt-resmooth-lib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true")
suppressMessages(library(frmtmb))
cat("frmtmb", as.character(packageVersion("frmtmb")), "from",
    find.package("frmtmb"), "\n")
bd <- readRDS("dev/resmooth-brms-design.rds")
d <- bd$data

forms <- list(re = bf(y ~ s(x) + s(g, bs = "re")),
              fs = bf(y ~ s(x, g, bs = "fs", k = 5)),
              t2 = bf(y ~ t2(x, g, bs = c("cr", "re"))),
              sg = bf(y ~ s(x) + (1 | g)))

# relative residual of lm(target ~ D), against the centred target
resid_rel <- function(D, target) {
  D <- cbind(1, as.matrix(D))
  q <- qr(D, tol = 1e-10)
  r <- as.vector(qr.resid(q, target))
  sqrt(sum(r^2)) / sqrt(sum((target - mean(target))^2))
}

for (nm in names(forms)) {
  fit <- suppressWarnings(frm(forms[[nm]], data = d))
  lp <- fit$frame$linpreds[[1L]]
  bl <- fit$frame$re_blocks
  sm <- which(vapply(bl, function(b) b[["covstruct"]] == "smooth", NA))
  grp <- unique(unlist(lapply(fit$frame$linpreds,
                              frmtmb:::smooth_group_block_ids)))
  # what re_formula = NA spans now, and what it spanned through 0.64.0
  Znow <- if (length(sm)) {
    as.matrix(lp[["Z"]][, unlist(lapply(bl[sm], `[[`, "c_idx")),
                        drop = FALSE])
  }
  keep_old <- setdiff(sm, grp)
  Zold <- if (length(keep_old)) {
    as.matrix(lp[["Z"]][, unlist(lapply(bl[keep_old], `[[`, "c_idx")),
                        drop = FALSE])
  }
  Xf <- as.matrix(lp[["X"]])
  Dnow <- if (is.null(Znow)) Xf else cbind(Xf, Znow)
  Dold <- if (is.null(Zold)) Xf else cbind(Xf, Zold)
  tg <- bd[[nm]]$epred_na
  # brms's own design, for a symmetric statement
  Db <- do.call(cbind, lapply(bd[[nm]]$pieces, as.matrix))
  cat(sprintf(paste0("%-3s frmtmb cols X %d + Zsm %d (0.64.0 Zsm %d) | ",
                     "brms cols %d\n"),
              nm, ncol(Xf), if (is.null(Znow)) 0L else ncol(Znow),
              if (is.null(Zold)) 0L else ncol(Zold), ncol(Db)))
  cat(sprintf(paste0("    resid of brms epred(NA) on frmtmb NA design: ",
                     "now %.3e | 0.64.0 %.3e | on brms's own design ",
                     "%.3e\n"),
              resid_rel(Dnow, tg), resid_rel(Dold, tg), resid_rel(Db, tg)))
  # and the reverse direction: frmtmb's own eta(NA) in brms's span
  eta <- as.vector(frm_linpred(fit, re_formula = NA))
  cat(sprintf("    resid of frmtmb eta(NA) on brms design: %.3e\n",
              resid_rel(Db, eta)))
  cat(sprintf("    brms epred: identical(NA, NULL) %s max|NA-NULL| %.3e\n",
              identical(bd[[nm]]$epred_na, bd[[nm]]$epred_null),
              max(abs(bd[[nm]]$epred_na - bd[[nm]]$epred_null))))
}
