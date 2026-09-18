.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressPackageStartupMessages(library(brms))
cat("brms", format(packageVersion("brms")), "\n")
print(args(brms:::validate_newdata))

# example4 has two grouping factors: subject and (from its formula) more
for (nm in paste0("brmsfit_example", 1:5)) {
  x <- get(nm, envir = asNamespace("brms"))
  gr <- names(x$ranef)
  grv <- unique(x$ranef$group)
  cat("\n", nm, " groups:", paste(grv, collapse = ","), "\n")
  if (length(grv) < 1) next
  nd <- x$data[1:3, , drop = FALSE]
  drop1 <- grv[1]
  if (!drop1 %in% names(nd)) { cat("   not a column; skip\n"); next }
  nd[[drop1]] <- NULL
  # re_formula that drops the block whose column is missing
  keep <- setdiff(grv, drop1)
  rf <- if (length(keep))
    stats::as.formula(paste("~", paste0("(1|", keep, ")",
                                        collapse = " + "))) else NA
  r <- tryCatch(brms:::validate_newdata(nd, x, re_formula = rf,
                                        allow_new_levels = FALSE),
                error = function(e) e)
  cat("   drop", drop1, "with re_formula =", deparse(rf), "->",
      if (inherits(r, "condition"))
        paste("ERR:", substr(conditionMessage(r), 1, 80)) else
        paste("ACCEPTED, cols:", paste(names(r), collapse = ",")), "\n")
  r2 <- tryCatch(brms:::validate_newdata(nd, x, re_formula = NA,
                                         allow_new_levels = FALSE),
                 error = function(e) e)
  cat("   drop", drop1, "with re_formula = NA ->",
      if (inherits(r2, "condition"))
        paste("ERR:", substr(conditionMessage(r2), 1, 80)) else
        "ACCEPTED", "\n")
}
