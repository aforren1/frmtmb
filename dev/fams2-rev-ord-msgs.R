l <- readRDS("C:/Users/adf44/source/r/frmtmb-wt-fams2/dev/fams2-rev-out/ord-lane.rds")
for (m in c("cum_logit_re", "sratio_cs", "sratio_cs_num", "mv", "cum_thres_gr")) {
  it <- l[[m]]$post
  for (k in names(it)) if (inherits(it[[k]]$value, "cap_error"))
    cat(m, k, ":", substr(it[[k]]$value$msg, 1, 160), "\n")
}
# which outputs carry warnings
for (m in names(l)) { it <- l[[m]]$post; for (k in names(it)) if (length(it[[k]]$warnings)) cat("WARN", m, k, substr(it[[k]]$warnings[1], 1, 100), "\n") }
str(l$cum_logit_re$post$dharma$value)
str(l$cum_logit_re$post$ce$value, max.level = 2)
