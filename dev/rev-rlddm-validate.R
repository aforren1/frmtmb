# REVIEW: input validation on the exported ndt_bound(), which the help
# page invites a caller to use with a hand-assembled aterms list.
.libPaths(c("C:/Users/adf44/source/r/rev-rlddm-lib","C:/Users/adf44/source/r/pinlib","C:/Users/adf44/source/r/rellib-0552","C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb);library(frmtmb.eam)})
p <- function(lab, e) {
  r <- tryCatch({ v <- e; paste("ACCEPTED, ub =", format(v[["ub"]], digits = 9)) },
                error = function(err) paste("REFUSED:", substr(conditionMessage(err), 1, 70)),
                warning = function(w) paste("WARNING:", substr(conditionMessage(w), 1, 50)))
  cat(sprintf("%-38s %s\n", lab, r))
}
p("ndt_bound(c(0.31, 0.42))",        ndt_bound(c(0.31, 0.42)))
p("ndt_bound(numeric(0))",           ndt_bound(numeric(0)))
p("ndt_bound(c(0.3, NA))",           ndt_bound(c(0.3, NA)))
p("ndt_bound(c(-0.3, 0.4))",         ndt_bound(c(-0.3, 0.4)))
p("ndt_bound(c(0, 0.4))",            ndt_bound(c(0, 0.4)))
p("ndt_bound('a')",                  ndt_bound("a"))
p("ndt_bound(rt, max_ndt = 0)",      ndt_bound(c(0.31, 0.42), max_ndt = 0))
p("ndt_bound(rt, max_ndt = -1)",     ndt_bound(c(0.31, 0.42), max_ndt = -1))
p("ndt_bound(rt, max_ndt = c(1,2))", ndt_bound(c(0.31, 0.42), max_ndt = c(1, 2)))
p("ndt_bound(rt, max_ndt = NA)",     ndt_bound(c(0.31, 0.42), max_ndt = NA_real_))