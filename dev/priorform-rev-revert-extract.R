# Reviewer, lane wt-priorform: save the BASE build's versions of the
# functions the lane changed, so a lane process can put one back and
# show which test fails without it.
#   Rscript dev/priorform-rev-revert-extract.R
.libPaths(c("C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
ns <- asNamespace("frmtmb")
nms <- c("bf", "lf", "nlf", "parse_linpred", "assemble_frame", "as_priorlist",
         "set_prior", "print.frmtmb_priorlist", "check_dpar_name", "get_prior")
fns <- mget(nms, envir = ns)
saveRDS(fns, "C:/Users/adf44/source/r/frmtmb-wt-priorform/dev/priorform-rev-revert-base.rds")
cat("saved", length(fns), "\n")
