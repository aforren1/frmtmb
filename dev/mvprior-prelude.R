# Library order for lane wt-mvprior. The LANE arm reads the lane library
# first; the BASE arm reads the round's read-only build of 0.62.0 first.
# The user library is last in both, as a fallback only.
mvprior_arm <- Sys.getenv("MVPRIOR_ARM", "lane")
MVPRIOR_ROOT <- "C:/Users/adf44/source/r/frmtmb-wt-mvprior"
MVPRIOR_LIB <- "C:/Users/adf44/source/r/mvprior-lib"
.libPaths(c(if (identical(mvprior_arm, "lane")) MVPRIOR_LIB,
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = file.path(MVPRIOR_ROOT, "dev/stan-cache"),
           NOT_CRAN = "true")
