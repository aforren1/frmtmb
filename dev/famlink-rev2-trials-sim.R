## simulate() status for the mvbf trial cases, both arms.
ba <- readRDS("dev/famlink-rev2-trials-base.rds"); la <- readRDS("dev/famlink-rev2-trials-lane.rds")
for (i in c("mv_both_trials", "mv_mixed_families")) cat(i, "\n base:", ba[[i]]$simulate, "\n lane:", la[[i]]$simulate, "\n")
