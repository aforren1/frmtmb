# Reviewer: the gradient mismatches flagged by cdfs-rev-grad.R, with
# central differences at three step sizes per component, so a FD
# artifact (kink, curvature) can be told from an AD error. Same data
# and seeds as cdfs-rev-grad.R.
.libPaths(c("C:/Users/adf44/source/r/wt-cdfs-rev-lib",
            "C:/Users/adf44/source/r/rellib-r10",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
options(frmtmb.notices = FALSE)
src <- readLines("C:/Users/adf44/source/r/frmtmb-wt-cdfs/dev/cdfs-rev-grad.R")
cut <- grep("^pick <- ", src)
cases <- list(list("compois", 3, "nu", 4), list("student", 1, "sigma", -10),
              list("hurdle_lognormal", 2, "sigma", -10), list("shifted_lognormal", 1, "sigma", 10),
              list("zero_inflated_beta_binomial", 1, "phi", 10),
              list("compois", 3, "nu", 3))
for (cs in cases) {
  Sys.setenv(GROUP = cs[[2]])
  eval(parse(text = src[1:(cut - 1)]))
  nm <- cs[[1]]
  fm <- fams[[nm]]
  d <- fm$d
  tri <- nm %in% c("binomial", "beta_binomial", "zero_inflated_binomial",
                   "zero_inflated_beta_binomial")
  dpn <- setdiff(frmtmb:::as_frmtmb_family(fm$fam)$dpars, "mu")
  rhs <- paste(c("y |", if (tri) "trials(tr) +", "cens(cc, y2) ~ x"), collapse = " ")
  fo <- do.call(bf, c(list(as.formula(rhs)),
                      lapply(dpn, function(dp) as.formula(paste(dp, "~ 1")))))
  ob <- frm(fo, family = fm$fam, data = d, dry_run = "objective")$obj
  p <- ob$par
  jd <- which(names(p) == "betad")
  p[jd[match(cs[[3]], dpn)]] <- cs[[4]]
  g <- as.numeric(ob$gr(p))
  cat(sprintf("== %s, %s eta = %g: fn %.10g\n", nm, cs[[3]], cs[[4]], ob$fn(p)))
  for (j in seq_along(p)) {
    fds <- vapply(c(1e-3, 1e-5, 1e-7), function(h) {
      e <- replace(numeric(length(p)), j, h)
      (ob$fn(p + e) - ob$fn(p - e)) / (2 * h)
    }, 0)
    cat(sprintf("  par %d: AD %-14.8g FD(1e-3) %-14.8g FD(1e-5) %-14.8g FD(1e-7) %-14.8g\n",
                j, g[j], fds[1], fds[2], fds[3]))
  }
}
