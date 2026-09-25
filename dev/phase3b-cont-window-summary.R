# Summarise dev/phase3b-cont-window.R. Output: dev/phase3b-log/contwin.txt
f <- list.files("dev/phase3b-log/contwin", full.names = TRUE)
rows <- do.call(rbind, lapply(f, function(p) {
  x <- readRDS(p)
  do.call(rbind, lapply(x$res, function(r) {
    data.frame(design = x$design, seed = x$seed, truth = x$truth,
               window = r$window,
               est = if (is.null(r$est)) NA else r$est,
               lwr = if (is.null(r$lwr)) NA else r$lwr,
               upr = if (is.null(r$upr)) NA else r$upr,
               top = if (is.null(r$crange)) NA else r$crange[2],
               conv = if (is.null(r$conv)) NA else r$conv,
               err = if (is.null(r$error)) "" else r$error)
  }))
}))
rows$hit <- rows$lwr < qlogis(rows$truth) & qlogis(rows$truth) < rows$upr
out <- character(0)
for (g in split(rows, list(rows$design, rows$window), drop = TRUE)) {
  out <- c(out, sprintf(
    "%-8s %-9s seeds %2d errors %d  lambda mean %.4f (truth mean %.4f)  covered %d of %d  code 0 on %d  window top mean %.3f",
    g$design[1], g$window[1], nrow(g), sum(nzchar(g$err)),
    mean(plogis(g$est), na.rm = TRUE), mean(g$truth),
    sum(g$hit, na.rm = TRUE), sum(!is.na(g$hit)),
    sum(g$conv == 0, na.rm = TRUE), mean(g$top, na.rm = TRUE)))
}
writeLines(out, "dev/phase3b-log/contwin.txt")
cat(out, sep = "\n")
