## Join dev/famlink-rev-designs-{brms,base,lane}.rds: false alarms
## (brms ok, base ok, lane error), new refusals, message counts, and
## bitwise equality of logLik, estimates, SEs and fitted values where
## both frmtmb arms fit.
br <- readRDS("dev/famlink-rev-designs-brms.rds")
ba <- readRDS("dev/famlink-rev-designs-base.rds")
la <- readRDS("dev/famlink-rev-designs-lane.rds")
nms <- names(la)
tab <- data.frame(case = nms,
  brms = vapply(nms, function(k) br[[k]]$status %||% "?", ""),
  base = vapply(nms, function(k) ba[[k]]$status %||% "?", ""),
  lane = vapply(nms, function(k) la[[k]]$status %||% "?", ""),
  msg_brms = vapply(nms, function(k) br[[k]]$n_bern_msg %||% NA_integer_, 1L),
  msg_lane = vapply(nms, function(k) la[[k]]$n_bern_msg %||% NA_integer_, 1L),
  msg_lane_update = vapply(nms, function(k) la[[k]]$n_bern_msg_update %||% NA_integer_, 1L),
  row.names = NULL)
cat("constructions:", nrow(tab), "\n")
cat("false alarms (brms ok, base ok, lane error):",
    sum(tab$brms == "ok" & tab$base == "ok" & tab$lane == "error"), "\n")
cat("lane refuses where brms accepts (any base):",
    sum(tab$brms == "ok" & tab$lane == "error"), "\n")
cat("lane accepts where brms refuses:",
    sum(tab$brms == "error" & tab$lane == "ok"), "\n")
print(tab[tab$brms != tab$lane | tab$base != tab$lane |
          (!is.na(tab$msg_lane) & tab$msg_brms != tab$msg_lane) |
          (!is.na(tab$msg_lane_update) & tab$msg_lane_update > 0), ],
      row.names = FALSE)
both <- nms[tab$base == "ok" & tab$lane == "ok"]
cat("\nfitted in both frmtmb arms:", length(both), "\n")
nb <- 0L
for (k in both) {
  flds <- c("logLik", "est", "se", "fitted")
  same <- vapply(flds, function(f) identical(ba[[k]][[f]], la[[k]][[f]]), TRUE)
  if (!all(same)) {
    nb <- nb + 1L
    cat(sprintf("  %-28s not identical in: %s (max abs diff est %.3g)\n", k,
                paste(flds[!same], collapse = ","),
                max(abs(ba[[k]]$est - la[[k]]$est))))
  }
}
cat("fitted in both and bitwise identical in logLik, estimates, SEs, fitted:",
    length(both) - nb, "of", length(both), "\n")
