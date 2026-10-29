# Review of lane meme: brms 2.23.0 against frmtmb on the extended
# shapes, from dev/review-meme-out/{brms,frm}-ext.rds.
#   Rscript dev/review-meme-cmp-ext.R > dev/review-meme-out/cmp-ext.txt
out <- "C:/Users/adf44/source/r/frmtmb-wt-meme/dev/review-meme-out"
B <- readRDS(file.path(out, "brms-ext.rds"))
F <- readRDS(file.path(out, "frm-ext.rds"))
key <- function(t, with_prior = TRUE, me_only = FALSE) {
  if (is.character(t)) return(paste("ERROR", t))
  t <- as.data.frame(t)
  for (k in c("resp", "dpar", "nlpar", "coef", "group")) {
    if (is.null(t[[k]])) t[[k]] <- ""
  }
  if (me_only) t <- t[t$class %in% c("meanme", "sdme", "corme", "Lme"), ]
  lb <- suppressWarnings(as.numeric(as.character(t$lb)))
  ub <- suppressWarnings(as.numeric(as.character(t$ub)))
  pr <- ifelse(nzchar(t$prior), t$prior, "(flat)")
  paste(if (with_prior) pr, t$class, t$coef, t$group, t$resp, t$dpar,
        ifelse(is.na(lb), "", lb), ifelse(is.na(ub), "", ub), sep = "|")
}
cmp <- function(lab, b, f) {
  same <- identical(b, f)
  cat(sprintf("  %-38s %s\n", lab, if (same) "IDENTICAL" else "DIFFERENT"))
  if (!same) {
    ob <- setdiff(b, f)
    of <- setdiff(f, b)
    if (length(ob)) cat("    brms only:  ", ob, sep = "\n      ")
    if (length(of)) cat("    frmtmb only:", of, sep = "\n      ")
    if (!length(ob) && !length(of)) {
      cat("    same set, other order\n      brms:  ",
          paste(b, collapse = " ; "), "\n      frmtmb:",
          paste(f, collapse = " ; "), "\n")
    }
    cat("\n")
  }
  same
}
drop_b <- "^(Xme|lp__|lprior|Ymi|Yl|zme|Intercept$|Intercept_)"
drop_f <- "^(Xme|lp__|logsdme|thetame|miss_|Intercept$|Intercept_)"
tally <- c()
for (nm in names(B)) {
  cat("====", nm, "\n")
  b <- B[[nm]]
  f <- F[[nm]]
  if (is.character(b$default)) {
    cat("  brms default_prior ERROR:", b$default, "\n")
    cat("  frmtmb default_prior:", if (is.character(f$default)) f$default
        else "a table", "\n")
    cat("  brms standata:", if (is.character(b$standata) &&
                                length(b$standata) == 1) b$standata
        else "ok", "\n")
    next
  }
  if (is.character(b$standata) && length(b$standata) == 1) {
    cat("  brms standata ERROR:", b$standata, "\n")
  }
  if (is.character(f$default) && length(f$default) == 1) {
    cat("  frmtmb default_prior:", f$default, "\n")
  }
  tally[paste(nm, "default sample, me rows")] <-
    cmp("default_prior me rows (sample route)",
        key(b$default, me_only = TRUE), key(f$default_sample, me_only = TRUE))
  tally[paste(nm, "default sample, all rows")] <-
    cmp("default_prior all rows (sample route)", key(b$default),
        key(f$default_sample))
  tally[paste(nm, "default fit, keys")] <-
    cmp("default_prior all keys (fit route)", key(b$default, FALSE),
        key(f$default, FALSE))
  if (!is.null(b$validate)) {
    kb <- paste(key(b$validate), as.data.frame(b$validate)$source)
    kf <- paste(key(f$validate), as.data.frame(f$validate)$source)
    tally[paste(nm, "validate")] <- cmp("validate_prior (with source)", kb, kf)
  }
  if (!is.null(b$variables)) {
    tally[paste(nm, "variables")] <- cmp(
      "variables (draws)", grep(drop_b, b$variables, value = TRUE,
                                invert = TRUE),
      grep(drop_f, f$variables, value = TRUE, invert = TRUE))
    tally[paste(nm, "fit variables")] <- cmp(
      "variables (frm fit vs brms, no prior_)",
      grep(paste0(drop_b, "|^prior_"), b$variables, value = TRUE,
           invert = TRUE), f$fit_variables)
    tally[paste(nm, "prior_draws")] <- cmp("prior_draws names",
                                           b$prior_draws, f$prior_draws)
    tally[paste(nm, "prior_summary me")] <- cmp(
      "prior_summary me rows (brms Lme)",
      key(b$prior_summary, me_only = TRUE),
      key(f$prior_summary, me_only = TRUE))
  } else if (!is.null(b$fit_error)) {
    cat("  brms fit ERROR:", b$fit_error, "\n")
  }
  if (!is.null(f$fit_error)) cat("  frmtmb fit:", f$fit_error, "\n")
  if (!is.null(f$sample_error)) cat("  frmtmb sample:", f$sample_error, "\n")
}
cat(sprintf("\nRESULT identical %d of %d comparisons\n", sum(tally),
            length(tally)))
print(names(tally)[!tally])
