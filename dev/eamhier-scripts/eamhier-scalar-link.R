# Lane eamhier: on a GLOBAL bound, does predict(dpar = "ndt", type =
# "response") report a time or a fraction of the bound?
#
# It decides whether this lane's own records for the global-bound arms
# are right, and whether tests/testthat/test-scale.R's `eam-unbounded`
# and `eam-sv` global rows record what they say they record: that file
# computes the population non-decision time as
#
#   ndt_frac * mean(if (is.null(bound$floors)) bound$ub else floors)
#
# which is a second multiplication by the bound if the scalar link
# already returns a time.
#
# Run: Rscript --vanilla dev/eamhier-scripts/eamhier-scalar-link.R
source("dev/eamhier-scripts/eamhier-common.R")
eamhier_libs()
suppressMessages({
  library(frmtmb)
  library(frmtmb.eam)
})

d <- eamhier_data(20260910L, "A", ns = 6L, nt = 80L)
key <- d[match(levels(d$s), as.character(d$s)), , drop = FALSE]

probe <- function(tag, fam) {
  fit <- frm(eamhier_form("A", "gl"), family = fam, data = d,
             se = TRUE)
  bd <- frmtmb::single_response(fit)[["family"]][["ndt_bound"]]
  lk <- family(fit)[["links"]][["ndt"]]
  eta <- unname(unlist(fixef(fit))[["ndt.(Intercept)"]])
  p1 <- suppressWarnings(stats::predict(
    fit, newdata = key[1L, , drop = FALSE], dpar = "ndt",
    type = "response", re.form = NA))
  p2 <- suppressWarnings(stats::predict(
    fit, newdata = key[1L, , drop = FALSE], dpar = "ndt",
    type = "response", re.form = NA, se.fit = TRUE))
  cat("\n== ", tag, " ==\n", sep = "")
  cat("  link name          ", if (is.list(lk)) lk$name else lk, "\n")
  cat("  bound ub           ", bd[["ub"]], " floors NULL: ",
      is.null(bd[["floors"]]), "\n", sep = "")
  cat("  eta (population)   ", eta, "\n")
  cat("  linkinv(eta)       ",
      if (is.list(lk)) lk$linkinv(eta) else NA, "\n")
  cat("  predict re.form=NA ", as.numeric(p1)[1L], "\n")
  cat("  predict + se.fit   ", as.numeric(p2$fit)[1L], "\n")
  cat("  ndt_time(re.form=NA)",
      as.numeric(ndt_time(fit, newdata = key[1L, , drop = FALSE],
                          re.form = NA)), "\n")
  cat("  plogis(eta)        ", stats::plogis(eta), "\n")
  cat("  min(rt) in data    ", min(d$rt), "\n")
}

probe("default max_ndt (= min(rt))", wiener())
probe("max_ndt = 0.45, allow_unreachable",
      wiener(max_ndt = 0.45, allow_unreachable = TRUE))
