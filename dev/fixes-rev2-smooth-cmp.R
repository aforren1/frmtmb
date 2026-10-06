# Reviewer of lane fixes, re-check: summarize dev/fixes-rev2-smooth.R,
# base against lane, fit by fit.
#   Rscript dev/fixes-rev2-smooth-cmp.R
rd <- function(lib) {
  f <- list.files("dev/fixes-rev2-log", paste0("^smooth-", lib, "-.*rds$"),
                  full.names = TRUE)
  unlist(lapply(f, readRDS), recursive = FALSE)
}
L <- rd("wt-fixes-lib")
B <- rd("rellib-r5")
key <- function(r) paste(r$case, r$seed)
names(L) <- vapply(L, key, "")
names(B) <- vapply(B, key, "")
ks <- intersect(names(L), names(B))
cat("fits per build:", length(L), length(B), " paired:", length(ks), "\n")
cases <- unique(vapply(L, `[[`, "", "case"))
tab <- NULL
odd <- list()
for (cs in cases) {
  kk <- ks[startsWith(ks, paste0(cs, " "))]
  g <- function(X, f) vapply(X[kk], function(r) {
    if (!is.na(r$err)) NA else as.numeric(f(r))
  }, 0)
  err_l <- sum(vapply(L[kk], function(r) !is.na(r$err), NA))
  err_b <- sum(vapply(B[kk], function(r) !is.na(r$err), NA))
  conv_l <- g(L, function(r) r$conv != 0)
  conv_b <- g(B, function(r) r$conv != 0)
  se_l <- g(L, function(r) !r$se_ok)
  se_b <- g(B, function(r) !r$se_ok)
  dll <- g(L, function(r) r$ll) - g(B, function(r) r$ll)
  rel <- dll / abs(g(B, function(r) r$ll))
  # gap to mgcv's ML, each build
  ml_l <- g(L, function(r) r$ll - r$ml)
  ml_b <- g(B, function(r) r$ll - r$ml)
  fd <- vapply(kk, function(k) {
    if (!is.na(L[[k]]$err) || !is.na(B[[k]]$err)) return(NA_real_)
    max(abs(L[[k]]$fitted - B[[k]]$fitted)) / max(abs(B[[k]]$fitted))
  }, 0)
  tl <- g(L, function(r) r$time)
  tb <- g(B, function(r) r$time)
  tab <- rbind(tab, data.frame(
    case = cs, n = length(kk), err = paste(err_l, err_b, sep = "/"),
    conv_bad = paste(sum(conv_l, na.rm = TRUE), sum(conv_b, na.rm = TRUE),
                     sep = "/"),
    se_bad = paste(sum(se_l, na.rm = TRUE), sum(se_b, na.rm = TRUE),
                   sep = "/"),
    dll_min = signif(min(dll, na.rm = TRUE), 3),
    dll_max = signif(max(dll, na.rm = TRUE), 3),
    rel_max = signif(max(abs(rel), na.rm = TRUE), 3),
    ml_gap_l = signif(min(ml_l, na.rm = TRUE), 3),
    ml_gap_b = signif(min(ml_b, na.rm = TRUE), 3),
    fit_rel_max = signif(max(fd, na.rm = TRUE), 3),
    time_ratio = signif(sum(tl, na.rm = TRUE) / sum(tb, na.rm = TRUE), 3)))
  bad <- kk[which((conv_l | se_l | abs(dll) > 1e-4) %in% TRUE)]
  for (k in bad) odd[[k]] <- k
}
print(tab, row.names = FALSE)
cat("\nerr messages (lane):",
    unique(vapply(L, function(r) if (is.na(r$err)) "" else
      substr(r$err, 1, 90), "")), sep = "\n  ")
cat("\nfits with lane conv != 0, a non-finite SE, or |dlogLik| > 1e-4:\n")
for (k in names(odd)) {
  l <- L[[k]]; b <- B[[k]]
  cat(sprintf(paste0("%-8s lane conv %d se_ok %s ll %.6f | base conv %d ",
                     "se_ok %s ll %.6f | mgcv ML %.6f | lane theta %s | ",
                     "base theta %s\n"), k, l$conv, l$se_ok, l$ll, b$conv,
              b$se_ok, b$ll, l$ml,
              paste(round(l$theta, 2), collapse = ","),
              paste(round(b$theta, 2), collapse = ",")))
}
cat("\nbase-only bad fits:\n")
for (k in ks) {
  b <- B[[k]]; l <- L[[k]]
  if (!is.na(b$err) || !is.na(l$err)) next
  if ((b$conv != 0 || !b$se_ok) && !(k %in% names(odd))) {
    cat(sprintf("%-8s base conv %d se_ok %s ll %.6f | lane ll %.6f\n", k,
                b$conv, b$se_ok, b$ll, l$ll))
  }
}
