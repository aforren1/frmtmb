setwd("C:/Users/adf44/source/r/frmtmb-wt-sarfcor/extensions/frmtmb.spatial")
rep1 <- function(p, old, new) {
  txt <- paste(readLines(p), collapse = "\n")
  n <- lengths(regmatches(txt, gregexpr(old, txt, fixed = TRUE)))
  if (n != 1L) stop(p, " matches: ", n, " for ", substr(old, 1, 60))
  writeLines(strsplit(sub(old, new, txt, fixed = TRUE), "\n",
                      fixed = TRUE)[[1]], p)
}
rep1("R/zzz.R", "# rows in frm_compat(), through frmtmb_register_compat(). The structures are
# registered FIRST,", "# rows in frm_compat(), through frmtmb_register_compat(). The structures
# are registered FIRST,")
rep1("R/zzz.R", " sigma and nu must be constant, as in brms.\",
    override = TRUE)", " sigma and nu must be constant, as in brms.\")")
rep1("DESCRIPTION", "    rstan,\n", "    rstan,\n    spdep,\n")
rep1("DESCRIPTION", "    fixed data. Each is fitted", "    fixed data. sar() is the simultaneous autoregressive term of 'brms',
    in its lag and error forms, registered through the residual
    correlation seam. Each is fitted")
