## Pair dev/famlink-rev2-optrace/{base,lane}-<file>.rds record by record.
## An optimizer call is "identical" when par, objective, convergence,
## message, iterations and evaluations are all identical().
## Usage: Rscript dev/famlink-rev2-optrace-cmp.R
dir <- "dev/famlink-rev2-optrace"
tags <- unique(sub("^(base|lane)-", "", sub("[.]rds$", "", list.files(dir))))
flds <- c("par", "objective", "convergence", "message", "iterations", "evaluations")
tot <- c(calls = 0, same = 0, mapped_fits = 0, nonfinite = 0)
rows <- list()
for (t in sort(tags)) {
  fb <- file.path(dir, paste0("base-", t, ".rds")); fl <- file.path(dir, paste0("lane-", t, ".rds"))
  if (!file.exists(fb) || !file.exists(fl)) { cat(t, ": missing an arm\n"); next }
  b <- readRDS(fb); l <- readRDS(fl)
  nb <- length(b$records); nl <- length(l$records); m <- min(nb, nl)
  same <- 0L; first_diff <- NA_integer_
  for (i in seq_len(m)) {
    ok <- all(vapply(flds, function(f) identical(b$records[[i]][[f]], l$records[[i]][[f]]), TRUE))
    if (ok) same <- same + 1L else if (is.na(first_diff)) first_diff <- i
  }
  nf <- vapply(l$records, function(r) as.integer(r$nonfinite_trials %||% 0L), 1L)
  tot <- tot + c(nl, same, sum(nf > 0), sum(nf))
  rows[[t]] <- data.frame(file = t, calls_base = nb, calls_lane = nl, identical = same,
    first_diff = first_diff, lane_calls_nonfinite = sum(nf > 0), nonfinite_total = sum(nf),
    warn_base = b$summary[["warning"]], warn_lane = l$summary[["warning"]],
    fail_err_base = b$summary[["fail"]] + b$summary[["error"]],
    fail_err_lane = l$summary[["fail"]] + l$summary[["error"]])
}
tab <- do.call(rbind, rows)
print(tab, row.names = FALSE)
cat("\nfiles:", nrow(tab), " lane optimizer calls:", tot[["calls"]],
    " identical to base (paired by order, over files with equal call counts):",
    sum(tab$identical[tab$calls_base == tab$calls_lane]), "of",
    sum(tab$calls_lane[tab$calls_base == tab$calls_lane]), "\n")
cat("files with unequal call counts:", sum(tab$calls_base != tab$calls_lane), "\n")
cat("lane calls with nonfinite_trials > 0:", tot[["mapped_fits"]], " total mapped trials:", tot[["nonfinite"]], "\n")
