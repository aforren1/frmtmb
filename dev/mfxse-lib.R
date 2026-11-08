# Shared set-up of lane mfxse's scripts.
#   MFX_ARM  "after" (default): the lane's library, then rellib-r18;
#            "before": rellib-r18 alone (frmtmb 0.79.1 at 3cd2a0c8)
mfx_root <- "C:/Users/adf44/source/r/frmtmb-wt-mfxse"
mfx_user <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
mfx_base <- "C:/Users/adf44/source/r/rellib-r18"
mfx_lib <- "C:/Users/adf44/source/r/wt-mfxse-lib"
mfx_arm <- Sys.getenv("MFX_ARM", "after")
if (identical(mfx_arm, "before")) {
  .libPaths(c(mfx_base, mfx_user))
} else {
  .libPaths(c(mfx_lib, mfx_base, mfx_user))
}
Sys.setenv(FRMTMB_STAN_CACHE = file.path(mfx_root, "dev/stan-cache"),
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           NOT_CRAN = "true")
suppressPackageStartupMessages({
  library(frmtmb)
  library(marginaleffects)
})
cat("lib:", find.package("frmtmb"), format(packageVersion("frmtmb")),
    "arm:", mfx_arm, "\n")
quiet <- function(e) suppressWarnings(suppressMessages(e))
