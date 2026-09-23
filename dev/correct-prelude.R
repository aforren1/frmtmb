# Library order for lane wt-correct: the lane library first, the round's
# read-only base build second, the user library last.
# CORRECT_LIB lets a second worker on this lane run the same scripts
# against its own private library, which the R library rule requires.
.libPaths(c(Sys.getenv("CORRECT_LIB", "C:/Users/adf44/source/r/correct-lib"),
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE =
             "C:/Users/adf44/source/r/frmtmb-wt-correct/dev/stan-cache",
           NOT_CRAN = "true")
