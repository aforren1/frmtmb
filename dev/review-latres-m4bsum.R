# Summary of review-latres-m4b.R
dir <- "C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-out/m4b2"
fs <- list.files(dir, "tsv$", full.names = TRUE)
x <- do.call(rbind, lapply(fs, function(f) read.table(f, sep = "\t", header = TRUE, quote = "\"", stringsAsFactors = FALSE)))
tr <- c(ginv_b1=.1,ginv_b2=.2,ginv_b3=.3,gident_b=2,ig_b1=.3,ig_b2=.5,pid_b=4,pid_8=1.2,bid_b=.04,glog_c=.1,gpost_b=.5,wpost_b=.4,epost_b=.4,igpost_b=.4,bepost_b=.1,psqrt_b=.5,bb_post_lo=.5,bb_post_b=.1,zoib_post_lo=.5,zoib_post_b=.1,bb_T=.5,zoib_T=.5,gid_l=1,bid_l=.05,ginv_l=.15,pid_l=1,pid_20=2,glog_w=.05,blog_w=.2,ginv_s5=.3,gid_s5=2,ig_s5=.5,glog_s5=.15)
out <- do.call(rbind, lapply(split(x, x$design), function(d) {
  dn <- d$design[1]
  ok <- is.finite(d$ll) & is.finite(d$x_at_frm) & d$x_at_frm > -1e9
  data.frame(design = dn, n = nrow(d), code1 = sum(d$code != 0, na.rm = TRUE), nw = sum(d$nwarn > 0, na.rm = TRUE),
    s_true = tr[dn], s_frm = median(exp(d$ls), na.rm = TRUE), s_ex = median(exp(d$x_ls), na.rm = TRUE),
    ratio = median(exp(d$ls) / exp(d$x_ls), na.rm = TRUE),
    ar_frm = median(d$ar, na.rm = TRUE), ar_ex = median(d$x_ar, na.rm = TRUE),
    disp_frm = median(exp(d$ldisp), na.rm = TRUE), disp_ex = median(exp(d$x_ldisp), na.rm = TRUE),
    b1_frm = median(d$b1, na.rm = TRUE), b1_ex = median(d$x_b1, na.rm = TRUE),
    gap = median((d$ll - d$x_at_frm)[ok]), gapmax = max(abs(d$ll - d$x_at_frm)[ok]),
    exgap = median((d$x_ll - d$x_at_frm)[ok]),
    span = max(d$span, na.rm = TRUE), nan = sum(d$nan > 0, na.rm = TRUE),
    margin = median(d$margin_min, na.rm = TRUE), margin_lo = min(d$margin_min, na.rm = TRUE),
    cov_ls = mean(abs(d$ls - log(tr[dn])) <= 1.96 * d$se_ls, na.rm = TRUE), row.names = NULL)
}))
print(out[order(out$margin), ], digits = 3, row.names = FALSE, right = TRUE)
write.table(out, "C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-out/m4bsum.tsv", sep = "\t", row.names = FALSE)
