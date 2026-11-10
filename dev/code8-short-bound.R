# Code-8 stops short by more than the bound (2e-5), over both forced
# batteries on the final build: how many keep the warning.
lab <- commandArgs(TRUE)
if (!length(lab)) lab <- c("after2", "after-x2")
X <- rbind(read.delim(sprintf("dev/code8-out/short-%s.tsv", lab[1]), quote = ""),
           read.delim(sprintf("dev/code8-out/short-%s.tsv", lab[2]), quote = ""))
c8 <- X[grepl("false", X$msg), ]
far <- c8[!is.na(c8$dll) & c8$dll < -2e-5, ]
near <- c8[!is.na(c8$dll) & c8$dll >= -2e-5 & c8$dll < -1e-6, ]
cat("code-8 stops", nrow(c8), "; short by more than 2e-5:", nrow(far),
    "warned", sum(far$w_code > 0), "message", sum(far$m_atmax > 0), "\n")
cat("short by 1e-6 to 2e-5:", nrow(near), "message", sum(near$m_atmax > 0),
    "warned", sum(near$w_code > 0), "; their gain/shortfall ratio",
    format(range(near$gain / -near$dll, na.rm = TRUE), digits = 4), "\n")
r <- far$gain / -far$dll
cat("gain/shortfall on the far ones:", format(range(r, na.rm = TRUE),
                                               digits = 4), "\n")
