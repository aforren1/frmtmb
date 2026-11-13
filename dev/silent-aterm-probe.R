# Item 4675: two registered addition terms that share a key.
# Usage: Rscript dev/silent-aterm-probe.R base|lane
arm <- commandArgs(TRUE)[1]
libs <- c("C:/Users/adf44/source/r/rellib-r20",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (identical(arm, "lane")) {
  libs <- c("C:/Users/adf44/source/r/wt-silent-lib", libs)
}
.libPaths(libs)
suppressPackageStartupMessages(library(frmtmb))
cat("arm", arm, "frmtmb from", find.package("frmtmb"), "\n")
try_reg <- function(...) {
  tryCatch({frmtmb_register_aterm(...); "registered"},
           error = function(e) paste("REFUSED:", conditionMessage(e)))
}
show_spec <- function(f) {
  r <- tryCatch({
    sp <- frm(bf(f) + gaussian(), data = NULL, dry_run = "spec")
    a <- sp$responses[[1]]$aterms
    paste(names(a), vapply(a, deparse1, ""), sep = "=", collapse = "  ")
  }, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat("  ", deparse1(f), "->", r, "\n")
}
local({
  frmtmb_local_registries()
  cat("register x (arity 3):", try_reg("x", arity = 3), "\n")
  cat("register x2:", try_reg("x2"), "\n")
  show_spec(y | x2(z) + x(u, v, w) ~ 1)
  show_spec(y | x(u, v, w) + x2(z) ~ 1)
})
local({
  frmtmb_local_registries()
  cat("register vint1:", try_reg("vint1"), "\n")
  cat("register cens_y (arity 2):", try_reg("cens_y", arity = 2), "\n")
})
