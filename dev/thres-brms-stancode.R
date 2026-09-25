# Lane thres: brms 2.23.0's Stan code and Stan data for thres(), read
# by hand. Nothing is compiled.
#
#   Rscript dev/thres-brms-stancode.R > dev/thres-brms-stancode-log.txt
#
# What the log shows, and where frmtmb follows it:
# - grouped thresholds are one `ordered` (cumulative, sratio) or plain
#   vector per level, merged into merged_Intercept; row n reads
#   merged_Intercept[Jthres[n, 1]:Jthres[n, 2]] (R/thres.R);
# - the design is NOT centered under grouped thresholds (`mu += X * b`,
#   b_Intercept_g = Intercept_g), so frmtmb's threshold prior carries no
#   centering offset there (R/priors.R);
# - cs() with grouped thresholds is refused;
# - the prior rows are class Intercept with group = <level>;
# - an ordered factor with unused levels counts only the levels that
#   occur, per group too.
.libPaths(c("/opt/rlib/base", "/opt/rlib/deps", "/opt/r/lib/R/library"))
suppressMessages(library(brms))
set.seed(1)
d <- data.frame(y = sample(1:4, 60, TRUE), x = rnorm(60),
                g = rep(c("a", "b", "c"), 20))
d$y[d$g == "c"] <- pmin(d$y[d$g == "c"], 3)
d$nth <- ifelse(d$g == "c", 2, 3)
cat(make_stancode(bf(y | thres(gr = g) ~ x), data = d,
                  family = cumulative()))
cat("\n=========== with nth\n")
sd <- make_standata(bf(y | thres(nth, gr = g) ~ x), data = d,
                    family = cumulative())
str(sd[grepl("thres", names(sd), ignore.case = TRUE)])
print(get_prior(bf(y | thres(nth, gr = g) ~ x), data = d,
                family = cumulative()))
cat("\n=========== thres(5), sratio\n")
cat(make_stancode(bf(y | thres(5) ~ x), data = d, family = sratio()))
cat("\n=========== cs() with grouped thresholds\n")
print(tryCatch(make_stancode(bf(y | thres(gr = g) ~ cs(x)), data = d,
                             family = sratio()),
               error = conditionMessage))

cat("\n=========== ordered factor with unused levels\n")
set.seed(1)
d <- data.frame(y = sample(1:3, 40, TRUE), x = rnorm(40),
                g = rep(c("a", "b"), 20))
d$y[d$g == "b"] <- pmin(d$y[d$g == "b"], 2)
d$f <- factor(d$y, levels = 1:4, ordered = TRUE)
d$fm <- factor(c(1, 2, 4)[d$y], levels = 1:4, ordered = TRUE)
s1 <- make_standata(f ~ x, data = d, family = cumulative())
cat("top level unused, ungrouped nthres:", s1$nthres, "\n")
s2 <- make_standata(f | thres(gr = g) ~ x, data = d,
                    family = cumulative())
cat("top level unused, grouped nthres:", s2$nthres, "\n")
s3 <- make_standata(fm ~ x, data = d, family = cumulative())
cat("middle level unused nthres:", s3$nthres, "\n")
