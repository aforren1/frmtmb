# Reviewer: build a copy of the lane's run-vignette.R with the
# per-expression set.seed() removed, so the RNG flows through a
# vignette as it did in the 0.34.0 runner. Tests the claim "Before and
# after seeding, the 0.67.0 headline counts were identical."
#
#   Rscript dev/vigport-rev-mknoseed.R
root <- "C:/Users/adf44/source/r/frmtmb-wt-vigport/dev"
x <- readLines(file.path(root, "brms-port/run-vignette.R"), warn = FALSE)
i <- grep("^HERE <- local[(][{]$", x)
stopifnot(length(i) == 1, x[i + 3] == "})")
x <- c(x[seq_len(i - 1)],
       sprintf('HERE <- "%s/brms-port"', root),
       x[(i + 4):length(x)])
s <- grep("set.seed(port_seed(id))", x, fixed = TRUE)
stopifnot(length(s) == 1)
x[s] <- "  # reviewer: no per-expression seed"
dir.create(file.path(root, "vigport-rev-out"), showWarnings = FALSE)
writeLines(x, file.path(root, "vigport-rev-out/run-vignette-noseed.R"))
cat("wrote; removed line", s, "\n")
