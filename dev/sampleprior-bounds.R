# Lane sampleprior: brms's bound on each class that gets a natural-scale
# default against the image of frmtmb's link for that dpar, read off the
# link objects (linkinv at -30 and 30 on the link scale).
#   Rscript dev/sampleprior-bounds.R > dev/sampleprior-log/bounds.txt
.libPaths(c("C:/Users/adf44/source/r/wt-sampleprior-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.eam)
})
ln <- utils::read.delim(
  "C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/sampleprior-log/table-lane.tsv",
  colClasses = "character", na.strings = "NA", quote = "")
k <- ln$variant == "base" & !is.na(ln$brms_lb) &
  grepl("^[a-z]+[0-9]*[|]{5}$", ln$key) & !startsWith(ln$key, "sd|") &
  !startsWith(ln$key, "theta|") & ln$brms != "(flat)" &
  !grepl("^(lncor|sbhaz)[|]", ln$key)
z <- unique(ln[k, c("family", "key", "brms", "brms_lb", "brms_ub")])
z$dpar <- sub("[|]+$", "", z$key)
fam_of <- function(nm) {
  # the fastest response of the table's wiener data (table-lane.tsv)
  if (identical(nm, "wiener")) {
    return(frmtmb.eam::wiener(max_ndt = 0.354641719488436))
  }
  if (identical(nm, "mixture")) {
    return(frmtmb::mixture(stats::gaussian(), stats::gaussian()))
  }
  frmtmb::brmsfamily(nm)
}
img <- vapply(seq_len(nrow(z)), function(i) {
  f <- tryCatch(fam_of(z$family[i]), error = function(e) NULL)
  dp <- z$dpar[i]
  if (identical(z$family[i], "com_poisson") && dp == "shape") dp <- "nu"
  l <- f$links[[dp]]
  if (is.null(l) || !is.function(l$linkinv)) {
    return(sprintf("%s: link set from the data", l$name %||% "?"))
  }
  tryCatch(sprintf("%s: (%.4g, %.4g)", l$name, l$linkinv(-30),
                   l$linkinv(30)), error = function(e) conditionMessage(e))
}, "")
out <- data.frame(family = z$family, dpar = z$dpar, brms = z$brms,
                  brms_bounds = paste0("[", z$brms_lb, ", ", z$brms_ub, "]"),
                  frmtmb_link_image = img)
options(width = 200)
print(out[order(out$dpar, out$family), ], row.names = FALSE, right = FALSE)
