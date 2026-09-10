# REVIEW, item 1.0b, attack 3: the seam is an EXPORT, so the case that
# matters is the one where the caller does not do what the prose says.
#
#   Rscript dev/rev-rlddm-seam-abuse.R <lib>
#
# Seed 4242. Three questions.
#
# 1. ndt_bound_attach() refuses wiener(), lba(), rdm() and wiener_gng()
#    by name. What about gddm(), which is in the same package, has an
#    `ndt`, and whose density does NOT multiply a floor out?
# 2. What does a consumer who copies the three documented lines and
#    NOT the refusal get?
# 3. Does the exported bound object survive being hand-edited, which is
#    what a caller holding an API object will eventually do?

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

set.seed(4242)
rt <- c(0.31, 0.42, 0.55, 0.61, 0.78)
g <- c("a", "a", "a", "b", "b")
bd <- ndt_bound(rt, list(ndt_group = ndt_bound_key(g)))

say <- function(lab, e) {
  out <- tryCatch({ e; "ACCEPTED" },
                  error = function(err)
                    paste("REFUSED:", substr(conditionMessage(err), 1, 150)))
  cat(sprintf("%-34s %s\n", lab, out))
}

cat("-- 1. which families the seam refuses to bound --\n")
for (nm in c("wiener", "lba", "rdm", "wiener_gng", "gddm")) {
  f <- tryCatch(switch(nm,
                       wiener = wiener(), lba = lba(2), rdm = rdm(2),
                       wiener_gng = wiener_gng(), gddm = gddm()),
                error = function(e) e)
  if (inherits(f, "error")) {
    cat(sprintf("%-34s (could not construct: %s)\n", nm,
                substr(conditionMessage(f), 1, 60)))
    next
  }
  has_raw <- !is.null(f[["ndt_raw"]])
  has_ndt <- "ndt" %in% names(f[["links"]])
  cat(sprintf("%-16s ndt_raw=%-5s ndt link=%-5s  ", nm, has_raw, has_ndt))
  say("", ndt_bound_attach(f, bd))
}

cat("\n-- 2. the three documented lines, with and without the refusal --\n")
# exactly what R/ndt-seam.R's @section tells a consumer to write
three_lines <- function(dpars, aterms) {
  ndt <- dpars[["ndt"]]
  fl <- aterms[["ndt_floor"]]
  if (!is.null(fl)) ndt <- ndt * fl
  ndt
}
cat("grouped, floor present :",
    format(three_lines(list(ndt = 0.5),
                       list(ndt_group = c(1, 2),
                            ndt_floor = c(0.4, 0.6))), digits = 6), "\n")
cat("grouped, floor MISSING :",
    format(three_lines(list(ndt = 0.5), list(ndt_group = c(1, 2))),
           digits = 6),
    " <- a fraction returned as a time, silently\n")
cat("what rlddm's own copy does instead:",
    tryCatch(format(frmtmb.eam:::ddm_ndt_scaler(
      c(a = 0.4, b = 0.6), "x")(list(ndt_group = c("a", "b")))),
      error = function(e) conditionMessage(e)), "\n")

cat("\n-- 3. a hand-edited bound object --\n")
bad <- bd
bad[["floors"]] <- c(a = 9, b = 9)      # floors above every response
f <- gddm()
say("attach a floors-above-data bound", {
  fam <- tryCatch(ndt_bound_attach(wiener(), bad), error = function(e) e)
  if (inherits(fam, "error")) stop(conditionMessage(fam))
  fam
})
bad2 <- structure(list(ub = 0.31, floors = NULL, pending = FALSE,
                       sizes = NULL, what = "hand"),
                  class = "frmtmb_eam_ndt_bound")
say("a hand-built object with the class", {
  fam <- ln <- NULL
  ndt_bound_attach(structure(list(links = list(ndt = "log"),
                                  init_dpars = list()),
                             class = "frmtmb_family"), bad2)
})
