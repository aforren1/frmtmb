# Reviewer check (lane ceplot): fitted(sample_new_levels = "old_levels")
# on y ~ x + (1 + x || g), which frmtmb holds as two blocks of ONE
# grouping factor. brms chooses one seen level per new level of a
# grouping factor (get_new_rdraws() runs once per group), so the new
# group's intercept and slope come from the same seen group. Here: does
# the answer at x = 0 and x = 1 equal some seen group's answer?
#   Rscript dev/ceplot-rev-oldlevels2.R > dev/ceplot-rev-log/oldlevels2.txt
# Data seed 3; calls at seeds 1 to 12.
.libPaths(c("C:/Users/adf44/source/r/wt-ceplot-lib",
            "C:/Users/adf44/source/r/rellib-r4",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(3)
ng <- 12L
d <- data.frame(x = rnorm(360), g = factor(rep(seq_len(ng), 30)))
d$y <- rnorm(360, 1 + rnorm(ng, 0, 2)[d$g] +
               (0.5 + rnorm(ng, 0, 1)[d$g]) * d$x, 0.5)
f <- frm(bf(y ~ x + (1 + x || g)), family = gaussian(), data = d)
cat("blocks:", length(f$frame$re_blocks), "\n")
nd <- data.frame(x = c(0, 1), g = factor(c("new", "new")))
seen <- lapply(seq_len(ng), function(k) {
  unname(fitted(f, newdata = data.frame(x = c(0, 1),
                                        g = factor(k, levels = 1:ng)))[, 1])
})
hit <- vapply(1:12, function(s) {
  set.seed(s)
  v <- unname(fitted(f, newdata = nd, allow_new_levels = TRUE,
                     sample_new_levels = "old_levels")[, 1])
  k <- which(vapply(seen, function(a) identical(a, v), NA))
  if (length(k)) k[1] else NA_integer_
}, 1L)
cat("seed 1..12: the (x = 0, x = 1) answer equals seen group:",
    ifelse(is.na(hit), "none", hit), "\n")
cat("answers that are one seen group's: ", sum(!is.na(hit)), "of 12\n")
