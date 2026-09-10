# PUNCH ROUND 1, B2. Which families ndt_bound_attach() refuses, and
# what the old guard let through.
#
#   Rscript dev/rlddm-scripts/rlddm-guard.R <lib>
#
# No seed: the table is exhaustive over this package's five families.
#
# The first guard was `!is.null(fam[["ndt_raw"]])`, the marker
# ddm_ndt_install() leaves. Only FOUR of the five families carry it.
# This prints the marker beside the outcome so that the correlation the
# old guard rested on, and the family it does not cover, are both
# visible rather than argued.

args <- commandArgs(trailingOnly = TRUE)
lib <- args[[1L]]
.libPaths(c(lib,
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/source/r/rellib-0552",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.eam)
})
cat("lib:", lib, " eam", format(packageVersion("frmtmb.eam")), "\n\n")

rt <- c(0.31, 0.42, 0.55, 0.61, 0.78)
at <- list(ndt_group = ndt_bound_key(c("a", "a", "a", "b", "b")))
bd <- ndt_bound(rt, at, what = "seam")

fams <- list(wiener = wiener(), lba = lba(2), rdm = rdm(2),
             wiener_gng = wiener_gng(), gddm = gddm())
cat(sprintf("%-11s %8s %-24s %s\n", "family", "ndt_raw", "ndt link",
            "attach"))
for (nm in names(fams)) {
  f <- fams[[nm]]
  lk <- f[["links"]][["ndt"]]
  lkn <- if (is.list(lk)) lk[["name"]] else as.character(lk)
  out <- tryCatch({
    ndt_bound_attach(f, bd)
    "ACCEPTED"
  }, error = function(e) "REFUSED by name")
  cat(sprintf("%-11s %8s %-24s %s\n", nm,
              !is.null(f[["ndt_raw"]]), lkn, out))
}

cat("\n-- what the OLD guard would have done to gddm(), by hand --\n")
# the old condition, reproduced rather than described
g <- gddm()
cat("  old guard would refuse it:", !is.null(g[["ndt_raw"]]), "\n")
cat("  ndt link before          :", g[["links"]][["ndt"]][["name"]], "\n")
cat("  ndt start before         :",
    format(g[["init_dpars"]][["ndt"]](rt, at), digits = 6), "\n")
cat("  its density reads ndt_floor:",
    any(grepl("ndt_floor", deparse(g[["lpdf"]]), fixed = TRUE)), "\n")

cat("\n-- a family from ANOTHER package is still accepted --\n")
foreign <- frmtmb::frmtmb_family(
  "notours", dpars = c("mu", "ndt"),
  links = list(mu = "identity", ndt = "log"), primary_dpars = "mu",
  type = "continuous",
  lpdf = function(y, dpars, aterms) {
    stats::dnorm(y - ndt_apply(dpars, aterms, "notours"), dpars[["mu"]],
                 1, log = TRUE)
  },
  init_dpars = list(mu = function(y, aterms) 0,
                    ndt = function(y, aterms) 0.5 * min(y)))
ok <- tryCatch({
  f2 <- ndt_bound_attach(foreign, bd)
  paste("ACCEPTED, ndt link now", f2[["links"]][["ndt"]],
        "and aterm_data adds",
        paste(names(f2[["aterm_data"]](rt, at)), collapse = ", "))
}, error = function(e) paste("REFUSED:", conditionMessage(e)))
cat(" ", ok, "\n")
