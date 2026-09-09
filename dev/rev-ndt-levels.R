# rev-ndt: a factor grouping is NOT safe on newdata.
#
# ddm_coerce_ndt_group() refuses a character or logical column, because
# it "would be coded on the levels PRESENT in whatever data frame it is
# evaluated in. newdata holding a subset of the groups would then pair a
# row with another group's bound, silently." It codes a factor with
# as.integer(), that is by the level INDEX, on the argument that "a
# data frame carries [the levels] through a subset".
#
# A data frame carries the level SET through a subset. It does not carry
# it through droplevels(), relevel(), factor() called again, or a
# newdata frame the user built themselves. Each of those permutes the
# indices, and the same defect follows.
#
# Seed 808, worktree build.

.libPaths(c(Sys.getenv("REV_LIB",
                       "C:/Users/adf44/source/r/rev-ndt-lib2"),
            "C:/Users/adf44/source/r/reflib-r2",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(frmtmb.eam)})

set.seed(808)
d <- ddm_simulate(900, mu = 1.0, bs = 1.4, ndt = 0.28, bias = 0.5)
d$g <- factor(rep(c("s1", "s2", "s3"), length.out = nrow(d)))
fl <- tapply(d$rt, d$g, min)
cat("group floors:",
    paste(names(fl), sprintf("%.6f", fl), sep = "=", collapse = "  "),
    "\n\n")

fit <- frm(bf(rt | dec(upper) + ndt_group(g) ~ 1, bs ~ 1, ndt ~ 1,
              bias = 0.5), family = wiener(), data = d)
bd <- frmtmb::single_response(fit)[["family"]][["ndt_bound"]]
cat("floors the fit captured, keyed by label code:",
    paste(names(bd$floors), sprintf("%.6f", bd$floors), sep = "=",
          collapse = "  "), "\n\n")

# three rows, one per group, in group order
one <- d[match(levels(d$g), as.character(d$g)), , drop = FALSE]
cat("reference: one row per group, levels as fitted\n")
cat("   ndt_time:", sprintf("%.6f", ndt_time(fit, one)), "\n")
cat("   fraction:",
    sprintf("%.6f", suppressWarnings(
      predict(fit, newdata = one, dpar = "ndt", type = "response"))),
    "\n\n")

# 1. droplevels(): the ordinary way a user takes a subset of subjects
sub <- d[d$g != "s1", , drop = FALSE]
sub <- droplevels(sub)
sub1 <- sub[match(levels(sub$g), as.character(sub$g)), , drop = FALSE]
cat("1. droplevels() after dropping group s1\n")
cat("   levels now:", paste(levels(sub1$g), collapse = " "), "\n")
cat("   ndt_time:", sprintf("%.6f", ndt_time(fit, sub1)), "\n")
cat("   these two rows are s2 and s3; their own floors are",
    sprintf("%.6f", fl[c("s2", "s3")]), "\n")
cat("   fitted mean response:",
    sprintf("%.6f", suppressWarnings(
      predict(fit, newdata = sub1, type = "response"))), "\n")
sub_keep <- d[d$g != "s1", , drop = FALSE]
sub_keep1 <- sub_keep[match(c("s2", "s3"), as.character(sub_keep$g)), ,
                      drop = FALSE]
cat("   WITHOUT droplevels(), the same two rows:\n")
cat("   ndt_time:", sprintf("%.6f", ndt_time(fit, sub_keep1)), "\n")
cat("   fitted mean response:",
    sprintf("%.6f", suppressWarnings(
      predict(fit, newdata = sub_keep1, type = "response"))), "\n\n")

# 2. relevel(): the same rows, one line of reordering
rl <- one
rl$g <- relevel(rl$g, ref = "s3")
cat("2. relevel(ref = 's3'), the same three rows\n")
cat("   levels now:", paste(levels(rl$g), collapse = " "), "\n")
cat("   ndt_time:", sprintf("%.6f", ndt_time(fit, rl)), "\n\n")

# 3. a newdata frame built by hand, which is what a prediction grid is
hand <- data.frame(g = factor(c("s3", "s2", "s1")),
                   upper = one$upper[c(3L, 2L, 1L)],
                   rt = one$rt[c(3L, 2L, 1L)])
cat("3. a hand-built prediction grid, groups s3, s2, s1\n")
cat("   levels now:", paste(levels(hand$g), collapse = " "), "\n")
cat("   ndt_time:", sprintf("%.6f", ndt_time(fit, hand)), "\n")
cat("   their own floors are   ",
    sprintf("%.6f", fl[c("s3", "s2", "s1")]), "\n\n")

# 4. and what a character column would have done, for comparison: it is
#    refused at fit time, so the refusal is doing its job and the
#    factor path is the one left open
cat("4. the same grouping as a character column, at fit time\n")
d$gc <- as.character(d$g)
r <- tryCatch(frm(bf(rt | dec(upper) + ndt_group(gc) ~ 1, bs ~ 1,
                     ndt ~ 1, bias = 0.5), family = wiener(), data = d),
              error = function(e) substr(conditionMessage(e), 1, 70))
cat("  ", if (is.character(r)) paste("REFUSED:", r) else "FITTED", "\n")
