# Reviewer, claim 1: compare dev/fams2-rev-out/ord-{base,lane}.rds with
# identical(), output by output.
od <- "C:/Users/adf44/source/r/frmtmb-wt-fams2/dev/fams2-rev-out"
b <- readRDS(file.path(od, "ord-base.rds"))
l <- readRDS(file.path(od, "ord-lane.rds"))
stopifnot(identical(names(b), names(l)))
n_same <- 0; n_diff <- 0; n_err <- 0
for (m in names(b)) {
  bm <- b[[m]]; lm <- l[[m]]
  if (!is.null(bm$fit)) {
    fe <- !identical(bm$fit, lm$fit)
    if (fe) { cat("DIFF", m, "fit\n"); str(bm$fit); str(lm$fit) }
    if (inherits(bm$fit$value, "cap_error")) {
      cat(sprintf("%-15s fit ERROR both arms? %s: %s\n", m,
                  !fe, bm$fit$value$msg))
    }
  }
  items <- if (!is.null(bm$post)) bm$post else bm
  litems <- if (!is.null(lm$post)) lm$post else lm
  for (k in names(items)) {
    if (k == "fit") next
    same <- identical(items[[k]], litems[[k]])
    isErr <- inherits(items[[k]]$value, "cap_error")
    if (same) {
      n_same <- n_same + 1
      if (isErr) n_err <- n_err + 1
    } else {
      n_diff <- n_diff + 1
      cat("DIFF", m, k, "\n")
      ae <- all.equal(items[[k]], litems[[k]])
      print(head(ae, 5))
    }
  }
  errs <- names(items)[vapply(items, function(z) {
    inherits(z$value, "cap_error")
  }, TRUE)]
  if (length(errs)) {
    cat(sprintf("%-15s outputs that error in both arms: %s\n", m,
                paste(errs, collapse = ", ")))
  }
}
cat(sprintf("\nidentical outputs %d (of which errors identical in both arms %d); differing %d\n",
            n_same, n_err, n_diff))
