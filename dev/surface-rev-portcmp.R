# Reviewer: mechanical port spell pass, reviewer rerun vs the lane's
# lane and r6 results, for the vignettes rerun.
for (v in c("brms_missings", "brms_overview", "brms_multilevel")) {
  r <- lapply(c(r6 = "dev/surface-port-out/r6/results-spell",
                lane = "dev/surface-port-out/lane/results-spell",
                rev = "dev/surface-rev-port-out/results-spell"),
              function(d) readRDS(file.path(d, paste0(v, ".rds"))))
  st <- function(x) {
    if (is.data.frame(x)) return(x$status)
    vapply(x, function(e) as.character(e$status %||% e[["status"]] %||% NA), "")
  }
  `%||%` <- function(a, b) if (is.null(a)) b else a
  s <- lapply(r, st)
  cat(v, "rows", length(s$rev), "| OK r6", sum(s$r6 == "OK"), "lane",
      sum(s$lane == "OK"), "rev", sum(s$rev == "OK"),
      "| rev vs lane status differ:", sum(s$rev != s$lane), "\n")
}
