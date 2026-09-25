# Lane mv, 2026-09-25: the models that existed before this lane keep
# their objective BIT FOR BIT. Run once against the base library and
# once against the lane library, then diff the two outputs:
#   FRMTMB_LIB=base Rscript dev/mv-bitwise.R > base.txt
#   FRMTMB_LIB=/opt/rlib/lane-mv Rscript dev/mv-bitwise.R > lane.txt
#   diff base.txt lane.txt
# Values are printed in hexadecimal (%a), so any change in any bit
# shows. The point is the start template plus a fixed perturbation, so
# no optimizer enters.
lib <- Sys.getenv("FRMTMB_LIB", "")
.libPaths(c(if (nzchar(lib) && lib != "base") lib, "/opt/rlib/base",
            "/opt/rlib/deps", "/opt/r/lib/R/library"))
suppressMessages(library(frmtmb))
set.seed(1)
n <- 120
dd <- data.frame(x = rnorm(n), g = factor(rep(1:12, 10)),
                 t = rep(1:10, 12))
dd$y1 <- 1 + dd$x + rnorm(n)
dd$y2 <- -1 + 0.5 * dd$x + rnorm(n)
dd$yt <- dd$x + rt(n, 4)
dd$o <- cut(dd$x + rlogis(n), c(-Inf, -1, 0, 1, Inf), labels = FALSE)
models <- list(
  cumulative = bf(o ~ x + (1 | g)) + cumulative(),
  acat = bf(o ~ x) + acat(),
  mv_gauss = bf(y1 ~ x + (1 | g)) + bf(y2 ~ x + (1 | g)) + gaussian(),
  rescor = bf(y1 ~ x) + bf(y2 ~ x) + set_rescor(TRUE) + gaussian(),
  student_ar = bf(yt ~ x + ar(t, g, cov = TRUE)) + student()
)
for (nm in names(models)) {
  obj <- frm(models[[nm]], data = dd, dry_run = "objective")$obj
  set.seed(2)
  p <- obj$par + rnorm(length(obj$par), 0, 0.1)
  cat(nm, "fn", sprintf("%a", obj$fn(p)), "\n")
  cat(nm, "gr", sprintf("%a", as.numeric(obj$gr(p))), "\n")
}
