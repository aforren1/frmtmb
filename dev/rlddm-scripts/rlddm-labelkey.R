# The defect filed against item 1.0b: key the per-group floor table on
# the group's LABEL rather than on a digest of it.
#
#   Rscript dev/rlddm-scripts/rlddm-labelkey.R <lib>
#
# Seed 4242. The recorded obstacle is that the addition-term registry's
# coercion "must return numbers for the tape and is stateless". This
# measures the stronger fact underneath it: the label does not merely
# fail to be REMEMBERED, it is destroyed before any family can see it,
# because frmtmb wraps every registered coercion's result in
# as.numeric() at both call sites. So the question is not whether
# frmtmb.eam could hold a label set; it is whether frmtmb could deliver
# a label at all, and today it cannot.
#
# Constructed rather than argued: a term is registered whose coercion
# attaches the labels to its own return value, and the value is then
# read back everywhere a family or a prediction can see it.

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
d <- ddm_simulate(200, mu = 1.2, bs = 1.5, ndt = 0.25)
d$g <- factor(rep(c("a", "b"), length.out = nrow(d)))

# ndt_group() itself, re-registered with a coercion that TRIES to
# carry the label along with the code it must return
frmtmb_register_aterm("ndt_group", arity = 1L, coerce = function(x) {
  structure(ndt_bound_key(x), zz_labels = as.character(x))
})
fit <- frm(bf(rt | dec(upper) + ndt_group(g) ~ 1, bs ~ 1, ndt ~ 1,
              bias = 0.5), family = wiener(), data = d)
av <- fit$frame$aterm_values[["rt"]]
cat("aterm names in the fitted frame :",
    paste(names(av), collapse = ", "), "\n")
cat("attributes on the fitted value  :",
    if (length(attributes(av[["ndt_group"]]))) {
      paste(names(attributes(av[["ndt_group"]])), collapse = ", ")
    } else "NONE", "\n")

nd <- d[1:6, ]
anew <- frmtmb:::aterms_for_newdata(
  frmtmb::single_response(fit), nd)
cat("aterm names on newdata          :",
    paste(names(anew), collapse = ", "), "\n")
cat("attributes on the newdata value :",
    if (length(attributes(anew[["ndt_group"]]))) {
      paste(names(attributes(anew[["ndt_group"]])), collapse = ", ")
    } else "NONE", "\n")
cat("attributes after v[1:3]         :",
    if (length(attributes(av[["ndt_group"]][1:3]))) {
      paste(names(attributes(av[["ndt_group"]][1:3])), collapse = ", ")
    } else "NONE", "\n")

cat("\nthe two call sites that strip it:\n")
src <- function(f, pat) {
  z <- deparse(f)
  z[grepl(pat, z)]
}
cat("  frame.R  :", trimws(src(frmtmb:::assemble_frame,
                               "^\\s*as\\.numeric\\(v\\)$")[1L]), "\n")
cat("  predict.R:",
    trimws(src(frmtmb:::aterms_for_newdata, "as\\.numeric\\(")[1L]), "\n")

# and the residual that stays open, quoted from the source rather than
# re-derived: the code is two rolling hashes packed below 2^49
cat("\nthe code the floor table is keyed on:\n")
cat("  ndt_bound_key('a')        :",
    format(ndt_bound_key("a"), digits = 17), "\n")
cat("  below 2^49                :", ndt_bound_key("zzzzzzzz") < 2^49,
    "\n")
cat("  distinct codes over 10^6 identifiers S000000..S999999: ")
lab <- sprintf("S%06d", 0:999999)
co <- ndt_bound_key(lab)
cat(length(unique(co)), "of", length(co), "\n")
