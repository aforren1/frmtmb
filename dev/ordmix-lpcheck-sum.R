# Summary of dev/ordmix-lpcheck.R's logs (dev/ordmix-lpcheck-log/), the
# block the findings paste: per case, the largest ulp count over the
# four points, the largest relative difference, brms's gradient at
# frmtmb's optimum and frmtmb's own, and the cases that stopped.
root <- "C:/Users/adf44/source/r/frmtmb-wt-ordmix/dev/ordmix-lpcheck-log"
fs <- sort(list.files(root, pattern = "[.]txt$", full.names = TRUE))
num <- function(x, key) {
  as.numeric(sub(paste0("^.*", key, "=([^ ]+).*$"), "\\1", x))
}
rows <- list()
for (f in fs) {
  case <- sub("[.]txt$", "", basename(f))
  x <- readLines(f, warn = FALSE)
  lp <- grep("^LP ", x, value = TRUE)
  gr <- grep("^GRAD ", x, value = TRUE)
  gf <- grep("^GRADFRM ", x, value = TRUE)
  e <- readLines(sub("[.]txt$", ".err", f), warn = FALSE)
  err <- grep("^Error", e, value = TRUE)
  rows[[case]] <- data.frame(
    case = case, points = length(lp),
    nan_points = sum(grepl("brms=NaN", lp)),
    max_ulp = if (length(lp)) max(num(lp, "ulp"), na.rm = TRUE) else NA,
    max_rel = if (length(lp)) max(num(lp, "reldiff"), na.rm = TRUE) else NA,
    brms_grad = if (length(gr)) sub("^.*\\|=([^ ]+).*$", "\\1", gr) else NA,
    frm_grad = if (length(gf)) sub("^.*\\|=([^ ]+).*$", "\\1", gf) else NA,
    error = if (length(err)) substr(err[1], 1, 60) else "")
}
tab <- do.call(rbind, rows)
print(tab, row.names = FALSE)
ok <- tab$points == 4 & !nzchar(tab$error)
# the single-family cases localize brms's gradient; they test no code
# of this lane
single <- grepl("_1$", tab$case)
hurdle <- grepl("^hu_", tab$case)
mix <- !single & !hurdle
grp <- list(mixtures = mix, hurdle = hurdle, single = single)
cat(sprintf("\ncases run: %d; with all four points: %d\n", nrow(tab),
            sum(ok)))
for (g in names(grp)) {
  k <- ok & grp[[g]]
  cat(sprintf("%-9s cases %d of %d, largest ulp %.1f, largest relative %.3g\n",
              g, sum(k), sum(grp[[g]]), max(tab$max_ulp[k]),
              max(tab$max_rel[k])))
}
cat("points where both densities are NaN:", sum(tab$nan_points), "\n")
both <- vapply(fs, function(f) {
  x <- grep("^LP .*brms=NaN", readLines(f, warn = FALSE), value = TRUE)
  all(grepl("frm=NaN", x))
}, NA)
cat("and frmtmb NaN at each of them too:", all(both), "\n")
