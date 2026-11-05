# Counts quoted in the short check, generated from the benign set.
x <- do.call(rbind, lapply(list.files("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-out/m4b2", "tsv$", full.names = TRUE),
                           function(f) read.table(f, sep = "\t", header = TRUE, quote = "\"", stringsAsFactors = FALSE)))
x <- x[is.finite(x$ll), ]
x$ok <- abs(x$ls - x$x_ls) < log(1.25) & abs(x$ar - x$x_ar) < .3 & abs(x$ll - x$x_at_frm) < 1
print(tapply(x$ok, x$design, sum))
link <- c("ginv_b1","ginv_b2","ginv_b3","ginv_s5","gident_b","gid_s5","ig_b1","ig_b2","ig_s5","pid_b","pid_20","pid_8","bid_b","gid_l","bid_l","ginv_l","pid_l")
L <- x[x$design %in% link, ]
cat("link fits", nrow(L), " read true", sum(L$ok), " below margin 2", sum(L$margin_min < 2), " of them true", sum(L$ok[L$margin_min < 2]), " kept true share", mean(L$ok[L$margin_min >= 2]), "\n")
print(round(tapply(x$margin_min, x$design, median), 1))
s <- x[x$design == "glog_w", ]
cat("glog_w largest span", max(s$span, na.rm = TRUE), "\n")
g <- x[x$design %in% c("bb_T", "zoib_T"), ]
print(table(g$design, g$code != 0))
print(table(x$design[x$design %in% c("gpost_b", "wpost_b", "igpost_b", "epost_b", "bepost_b", "psqrt_b")]))
print(round(tapply(x$margin_min, x$design, min), 2)[c("gpost_b", "wpost_b", "igpost_b", "epost_b", "bepost_b", "psqrt_b")])
print(round(tapply(x$ar, x$design, median), 3)[c("psqrt_b", "bepost_b")])
print(round(tapply(x$x_ar, x$design, median), 3)[c("psqrt_b", "bepost_b")])
print(aggregate(cbind(zoib = 1) ~ design, x[x$design %in% c("bb_post_b", "zoib_post_b", "bb_post_lo", "zoib_post_lo"), ], length))
