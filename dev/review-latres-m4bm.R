# Is a condition on the fit better than one on the link? Pool the link-scale forms
# a refusal covers, one row per dataset, by the fit's own margin: the smallest
# distance of eta from the link's edge, in marginal sds of the latent.
dir <- "C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-out/m4b2"
fs <- list.files(dir, "tsv$", full.names = TRUE)
x <- do.call(rbind, lapply(fs, function(f) read.table(f, sep = "\t", header = TRUE, quote = "\"", stringsAsFactors = FALSE)))
link <- c("ginv_b1","ginv_b2","ginv_b3","ginv_s5","gident_b","gid_s5","ig_b1","ig_b2","ig_s5","pid_b","pid_20","pid_8","bid_b","gid_l","bid_l","ginv_l","pid_l")
held <- c("glog_c","glog_s5","glog_w","blog_w","bb_T","zoib_T")
truth <- c(ginv_b1=.1,ginv_b2=.2,ginv_b3=.3,ginv_s5=.3,gident_b=2,gid_s5=2,ig_b1=.3,ig_b2=.5,ig_s5=.5,pid_b=4,pid_20=2,pid_8=1.2,bid_b=.04,gid_l=1,bid_l=.05,ginv_l=.15,pid_l=1,glog_c=.1,glog_s5=.15,glog_w=.05,blog_w=.2,bb_T=.5,zoib_T=.5)
x <- x[x$design %in% c(link, held) & is.finite(x$ll), ]
x$rs <- exp(x$ls - x$x_ls)
x$dar <- x$ar - x$x_ar
x$over <- x$ll - x$x_at_frm
x$lerr <- x$ls - log(truth[x$design])
x$ok <- abs(log(x$rs)) < log(1.25) & abs(x$dar) < .3 & abs(x$over) < 1
x$cov <- abs(x$lerr) <= 1.96 * x$se_ls
L <- x[x$design %in% link, ]
L$band <- cut(L$margin_min, c(-Inf, 1, 2, 3, 4, 5, 6, Inf), right = FALSE)
tab <- do.call(rbind, lapply(split(L, L$band), function(d) data.frame(fits = nrow(d), designs = length(unique(d$design)),
  reads_true = round(mean(d$ok), 3), med_ratio = round(median(d$rs), 3), med_abs_ar = round(median(abs(d$dar)), 3),
  med_gap = round(median(d$over), 3), cov = round(mean(d$cov, na.rm = TRUE), 3), nan_starts = round(mean(d$nan > 0, na.rm = TRUE), 3), wide_span = round(mean(d$span > .01, na.rm = TRUE), 3))))
print(tab)
H <- x[x$design %in% held, ]
cat("held links, the same statistic, all fits:\n")
print(do.call(rbind, lapply(split(H, H$design), function(d) data.frame(fits = nrow(d), reads_true = round(mean(d$ok), 3), med_ratio = round(median(d$rs), 3), med_abs_ar = round(median(abs(d$dar)), 3), med_gap = round(median(d$over), 3), cov = round(mean(d$cov, na.rm = TRUE), 3)))))
cat("per link design, share of fits that read true, and share with margin >= 4:\n")
print(do.call(rbind, lapply(split(L, L$design), function(d) data.frame(fits = nrow(d), reads_true = round(mean(d$ok), 3), m4 = round(mean(d$margin_min >= 4), 3), true_if_m4 = round(mean(d$ok[d$margin_min >= 4]), 3), true_if_lt4 = round(mean(d$ok[d$margin_min < 4]), 3)))))
for (k in c(3, 4, 5)) {
  cat(sprintf("rule margin >= %d: keeps %d fits, of which read true %.3f; stops %d, of which read true %.3f (a refusal on a correct fit)\n", k, sum(L$margin_min >= k), mean(L$ok[L$margin_min >= k]), sum(L$margin_min < k), mean(L$ok[L$margin_min < k])))
}
