# Merge step 3 (quadseam2 x priors2), R/fit.R: the quadrature sd restart
# replaces `opt` after priors2's check_map_unbounded() has read the first
# one, so the restart's opt is read too before it is kept (quadseam2's
# findings, "After priors2 merges").
f <- "R/fit.R"
x <- readLines(f)
i <- which(x == "      opt <- nonfinite_not_converged(esc$opt, bounds$lower, bounds$upper)")
stopifnot(length(i) == 1)
x <- append(x, c(
  "      # the restart's own stop at an infinite-density end, which the",
  "      # check of the first opt above did not see",
  "      check_map_unbounded(esc$opt, bounds)"), after = i - 1)
con <- file(f, "wb"); writeLines(x, con, sep = "\r\n"); close(con)
cat("done\n")
