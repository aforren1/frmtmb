# Summary of the three random-slope calibrations (M1).
d <- "C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-log/"
r <- do.call(rbind, lapply(c("p1-slopecal.rds", "p1-slopecal2.rds",
                             "p1-slopecal3.rds"),
                           function(f) readRDS(paste0(d, f))))
cat("fits", nrow(r), "; autoscale = TRUE max short", signif(max(r$short_true), 3),
    "; TRUE short by more than 1e-5:", sum(r$short_true > 1e-5), "\n")
a <- aggregate(cbind(n = short_false) ~ s, r, length)
a$false_max <- aggregate(short_false ~ s, r, max)$short_false
a$false_over_1e3 <- aggregate(short_false ~ s, r, function(v) sum(v > 1e-3))$short_false
print(a[order(a$s), ], digits = 4)
