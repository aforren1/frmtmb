# Counts of dev/code8-short.R's second run (seeds 1 to 8, xf.tol 1e-3,
# 3e-3, 0.03): the gains of the code-8 stops, short and not.
X <- read.delim("dev/code8-out/short-after-x2.tsv", quote = "")
c8 <- X[grepl("false", X$msg), ]
sh <- c8[!is.na(c8$dll) & c8$dll < -1e-6, ]
ok <- c8[!(c8$case %in% sh$case), ]
cat("fits", nrow(X), "code 8", nrow(c8), "short (> 1e-6)", nrow(sh), "\n")
cat("short: gain", format(range(sh$gain, na.rm = TRUE), digits = 3),
    " smallest five gains", format(head(sort(sh$gain), 5), digits = 3),
    "\n  their dll", format(sh$dll[order(sh$gain)][1:5], digits = 3),
    "\n  warned", sum(sh$w_code > 0), "message", sum(sh$m_atmax > 0), "\n")
cat("not short:", nrow(ok), " gain", format(range(ok$gain, na.rm = TRUE),
                                             digits = 3),
    " message", sum(ok$m_atmax > 0), " warned", sum(ok$w_code > 0), "\n")
cat("short with gain in (1e-6, 1e-4]:", sum(sh$gain > 1e-6 & sh$gain <= 1e-4,
                                             na.rm = TRUE), "\n")
print(sh[sh$gain <= 1e-4, c("case", "gain", "dll", "w_code", "m_atmax")])
