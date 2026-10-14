# Second look at dev/famuni-out/wiener-density.rds (written by
# dev/famuni-wiener-density.R): where Stan's density is finite but off,
# where it is -Inf above the double underflow, and WienR's worst rows.
# Output: dev/famuni-out/wiener-density2.txt
g <- readRDS(
  "C:/Users/adf44/source/r/frmtmb-wt-famuni/dev/famuni-out/wiener-density.rds")
g$u <- g$td / g$a^2
cols <- c("t", "a", "w", "v", "up", "u", "ref", "eam", "stan", "wienr",
          "rtdists")
cat("-- Stan -Inf where the reference is above -700\n")
print(g[!is.finite(g$stan) & g$ref > -700, cols], digits = 10)
cat("\n-- Stan nonfinite by reference band\n")
br <- cut(g$ref, c(-Inf, -1e4, -1000, -745, -700, Inf))
print(table(band = br, stan_finite = is.finite(g$stan)))
cat("\n-- worst 10 finite Stan rows\n")
e <- abs(g$stan - g$ref); e[!is.finite(e)] <- -1
print(head(g[order(-e), cols], 10), digits = 10)
cat("\n-- Stan error distribution where finite\n")
ok <- is.finite(g$stan)
print(quantile(abs(g$stan - g$ref)[ok], c(0.5, 0.9, 0.99, 0.999, 1)))
cat("rows with Stan error above 1e-8:", sum(abs(g$stan - g$ref)[ok] > 1e-8),
    "of", sum(ok), "\n")
cat("\n-- worst 6 WienR rows\n")
e <- abs(g$wienr - g$ref); e[!is.finite(e)] <- Inf
print(head(g[order(-e), cols], 6), digits = 10)
cat("\n-- eam worst where reference above -700\n")
s <- g$ref > -700
e <- abs(g$eam - g$ref)
cat(sprintf("max abs %.3e  max rel %.3e  n %d\n", max(e[s]),
            max((e / pmax(1, abs(g$ref)))[s]), sum(s)))
cat(sprintf("eam over the whole grid: max abs %.3e at ref %.1f\n",
            max(e), g$ref[which.max(e)]))
