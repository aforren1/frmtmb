## Copy of dev/famlink-rev2-trials-cmp.R joining the punch round 2 lane rerun.
## Join dev/famlink-rev2-trials-{brms,base,lane}.rds.
br <- readRDS("dev/famlink-rev2-trials-brms.rds"); ba <- readRDS("dev/famlink-rev2-trials-base.rds"); la <- readRDS("dev/famlink-p2-trials-lane.rds")
k <- names(la)
tab <- data.frame(case = k, brms = sapply(k, function(i) br[[i]]$status), base = sapply(k, function(i) ba[[i]]$status),
  lane = sapply(k, function(i) la[[i]]$status), gp_brms = sapply(k, function(i) br[[i]]$get_prior),
  gp_lane = substr(sapply(k, function(i) la[[i]]$get_prior), 1, 5), row.names = NULL)
cat("constructions:", nrow(tab), "\n")
cat("lane refuses, brms accepts (frm):", sum(tab$brms == "ok" & tab$lane == "error"), "\n")
cat("lane accepts, brms refuses (frm):", sum(tab$brms == "error" & tab$lane == "ok"), " (", paste(tab$case[tab$brms == "error" & tab$lane == "ok"], collapse = ", "), ")\n")
cat("get_prior: lane refuses, brms accepts:", sum(tab$gp_brms == "ok" & tab$gp_lane == "error"), "\n")
print(tab, row.names = FALSE)
cat("\nfitted in both frmtmb arms, logLik and estimates identical; post-fit paths:\n")
for (i in k) if (ba[[i]]$status == "ok" && la[[i]]$status == "ok") {
  cat(sprintf("%-26s identical %-5s | predict(nd, no trials col) base %-10s lane %-10s | simulate lane %-4s | frm_simulate base %-6s lane %s\n", i,
    identical(ba[[i]]$ll, la[[i]]$ll) && identical(ba[[i]]$est, la[[i]]$est),
    substr(ba[[i]]$predict_newdata_no_trials, 1, 10), substr(la[[i]]$predict_newdata_no_trials, 1, 10),
    substr(la[[i]]$simulate, 1, 4), substr(ba[[i]]$frm_simulate %||% "-", 1, 6), substr(la[[i]]$frm_simulate %||% "-", 1, 60)))
}
for (i in k) if (!identical(ba[[i]]$predict_newdata_no_trials, la[[i]]$predict_newdata_no_trials) && la[[i]]$status == "ok") cat(i, "\n  base:", ba[[i]]$predict_newdata_no_trials, "\n  lane:", la[[i]]$predict_newdata_no_trials, "\n")
