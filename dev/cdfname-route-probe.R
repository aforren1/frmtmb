# Which check raises the rename refusal on each route of
# tests/testthat/test-cdf-slot-name.R: the frame route must raise from
# assemble_frame(), not from as_frmtmb_family().
#   Rscript dev/cdfname-route-probe.R
.libPaths(c("C:/Users/adf44/source/r/wt-cdfname-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), "\n")
pois_family <- function(...) {
  frmtmb_family(
    "pois_cdfname", dpars = "mu", links = list(mu = "log"),
    type = "discrete", accepts_aterms = c("weights", "cens", "trunc"),
    lpdf = function(y, dpars, aterms) {
      RTMB::dpois(y, dpars[["mu"]], log = TRUE)
    },
    init_dpars = list(mu = function(y, aterms) mean(y)), ...)
}
pois_cdf <- function(q, dpars, aterms) RTMB::ppois(q, dpars[["mu"]])
set.seed(3088)
d <- data.frame(y = stats::rpois(80, 4))
d$cc <- ifelse(d$y > 6, "right", "none")
d$y <- pmin(d$y, 6)
where <- function(label, expr) {
  site <- NA_character_
  # the calling handler inside the exiting one, so it sees the stack
  tryCatch(withCallingHandlers(
    expr,
    error = function(e) {
      if (!inherits(e, "frmtmb_cdf_slot_renamed")) return()
      fns <- vapply(sys.calls(), function(cl) {
        f <- cl[[1L]]
        if (is.name(f)) as.character(f) else ""
      }, "")
      known <- c("refuse_family_dots", "as_frmtmb_family",
                 "assemble_frame", "$.frmtmb_family", "mixture")
      hit <- intersect(rev(fns), known)
      site <<- if (length(hit)) hit[1L] else "?"
    }), error = function(e) NULL)
  cat(sprintf("ROUTE %-28s raised in %s\n", label, site))
}
where("frmtmb_family(lcdf=)", pois_family(lcdf = pois_cdf))
on <- pois_family(); on[["lcdf"]] <- pois_cdf
where("frm(family=)", frm(y | cens(cc) ~ 1, family = on, data = d))
where("bf() + family", bf(y | cens(cc) ~ 1) + on)
where("mixture()", mixture(on, poisson()))
b <- bf(y | cens(cc) ~ 1) + pois_family(cdf = pois_cdf)
b$family[["lcdf"]] <- pois_cdf
where("hand-edited bf", frm(b, data = d))
fz <- pois_family(cdf = pois_cdf, family_finalize = function(fam, y, av) {
  fam[["lcdf"]] <- fam[["cdf"]]
  fam
})
where("family_finalize()", frm(y ~ 1, family = fz, data = d))
where("$lcdf read", poisson_fam <- as_frmtmb_family(poisson())$lcdf)
