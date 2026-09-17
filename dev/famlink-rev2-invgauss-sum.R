## Summarize dev/famlink-rev2-invgauss.rds: warning kinds per design and link.
tab <- readRDS("dev/famlink-rev2-invgauss.rds")
tab$kind <- ifelse(tab$status != "fit", "error",
  ifelse(grepl("Hessian is not positive", tab$msg), "pdHess warning",
  ifelse(grepl("gradient", tab$msg), "gradient warning",
  ifelse(grepl("converge|false convergence|singular", tab$msg, ignore.case = TRUE), "convergence warning",
  ifelse(tab$warnings > 0, "other warning", "clean")))))
print(table(paste(tab$design, tab$link), tab$kind))
cat("\nlog-link replicates not clean:\n"); print(tab[tab$link == "log" & tab$kind != "clean", c("design", "rep", "status", "conv", "msg")], row.names = FALSE)
cat("\nexample messages, default link:\n"); print(unique(substr(tab$msg[tab$link == "default" & tab$kind != "clean"], 1, 80))[1:8])
cat("\nclean replicates: default", sum(tab$link == "default" & tab$kind == "clean"), "of", sum(tab$link == "default"),
    "; log", sum(tab$link == "log" & tab$kind == "clean"), "of", sum(tab$link == "log"), "\n")
