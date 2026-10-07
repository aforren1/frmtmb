# Summarize dev/setier-singular.R: per design and arm, how many fits
# lme4 calls singular, and what frmtmb told the user.
#   Rscript dev/setier-singular-sum.R <tsv> ...
for (f in commandArgs(TRUE)) {
  r <- read.delim(f, stringsAsFactors = FALSE)
  cat("==", basename(f), "\n")
  r$told <- ifelse(r$se_warn > 0, "se_warning",
                   ifelse(r$boundary_msg > 0, "boundary_msg", "nothing"))
  for (d in unique(r$design)) {
    x <- r[r$design == d, ]
    cat(sprintf(paste0("%-7s n=%3d lme4 singular %3d | lost %3d | ",
                       "boundary msg %3d (lme4 singular %3d) | SE warning ",
                       "%3d (lme4 singular %3d) | other warning %3d | ",
                       "VarCorr warning %3d | finite theta SE > 100 %3d | ",
                       "fit s %.2f\n"),
                d, nrow(x), sum(x$lme4_singular), sum(x$n_lost > 0),
                sum(x$boundary_msg > 0),
                sum(x$boundary_msg > 0 & x$lme4_singular),
                sum(x$se_warn > 0), sum(x$se_warn > 0 & x$lme4_singular),
                sum(x$other_warn > 0), sum(x$varcorr_warn > 0),
                sum(is.finite(x$se_theta_max) & x$se_theta_max > 100),
                median(x$fit_s)))
  }
  cat(sprintf("all: n=%d lme4 singular %d; boundary msg %d; SE warning %d\n",
              nrow(r), sum(r$lme4_singular), sum(r$boundary_msg > 0),
              sum(r$se_warn > 0)))
  r$lost[is.na(r$lost)] <- ""
  lt <- unlist(strsplit(as.character(r$lost[r$lost != ""]), ";"))
  lt <- sub("^theta_[0-9]+", "theta", lt)
  print(table(lt))
}
