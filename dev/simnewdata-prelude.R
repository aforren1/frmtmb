# Library order for lane wt-simnewdata: the lane library first, the
# round's read-only base build second, the user library last.
# SIMNEWDATA_LIB = "base" drops the lane library, so the same script
# measures the base build.
lane_lib <- Sys.getenv("SIMNEWDATA_LIB",
                       "C:/Users/adf44/source/r/simnewdata-lib")
.libPaths(c(if (!identical(lane_lib, "base")) lane_lib,
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-simnewdata/dev/stan-cache",
           NOT_CRAN = "true")
