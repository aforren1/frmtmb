## Item 4 of the review, second half: how often would the SHIPPED
## assertion fire on correct data?
##
## dev/coh-falsealarm.R measures the rate of the statistic that was
## REJECTED, and says so in its own header, but the plan's item 2.6 row
## quotes that rate as the shipped test's. This estimates the shipped
## statistic's rate from the 40 paired quotients the lane measured in
## dev/coh-paired.tsv.
.libPaths(c("C:/Users/adf44/source/r/cohrev-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
ln <- readLines("dev/coh-paired.tsv", warn = FALSE)
ln <- ln[nzchar(trimws(ln))]
field <- function(k) {
  as.numeric(vapply(strsplit(ln, "\t", fixed = TRUE), function(p) {
    v <- p[startsWith(p, paste0(k, "="))]
    trimws(sub(paste0("^", k, "="), "", v[1L]))
  }, character(1)))
}
s <- field("seed")
tr <- field("treat")
nu <- field("null")
q <- field("quotient")
cat(sprintf("rows %d, seeds %d to %d, all distinct: %s\n", length(ln),
            min(s), max(s), !anyDuplicated(s)))
cat(sprintf("quotient min %.4f median %.4f max %.4f; at or below 1: %d\n",
            min(q), median(q), max(q), sum(q <= 1)))
cat(sprintf("treat    min %.4f max %.4f\n", min(tr), max(tr)))
cat(sprintf("null     min %.4f max %.4f\n", min(nu), max(nu)))
cat(sprintf("the four pinned seeds 2610:2613: %s\n",
            paste(signif(q[s %in% 2610:2613], 4), collapse = " ")))
cat(sprintf("quotient recomputed from treat/null matches the column: %s\n",
            isTRUE(all.equal(q, tr / nu, tolerance = 1e-6))))

lq <- log(q)
p <- pnorm(0, mean(lq), sd(lq))
cat(sprintf("\nlog quotient: mean %.4f, sd %.4f, so 0 is %.2f sd below the mean\n",
            mean(lq), sd(lq), mean(lq) / sd(lq)))
cat(sprintf("lognormal tail P(quotient <= 1) per seed = %.3g\n", p))
cat(sprintf("the test takes the MIN of 4, so P(fire on correct data) = %.3g, about 1 in %.0f\n",
            1 - (1 - p)^4, 1 / (1 - (1 - p)^4)))
k <- 0
n <- length(q)
z <- qnorm(0.975)
cen <- (k / n + z^2 / (2 * n)) / (1 + z^2 / n)
hw <- z * sqrt((k / n) * (1 - k / n) / n + z^2 / (4 * n^2)) /
  (1 + z^2 / n)
cat(sprintf("model free: %d of %d at or below 1, Wilson 95%% upper bound %.4f per seed,\n",
            k, n, cen + hw))
cat(sprintf("  so a 4-seed min fires at most %.3f of the time at 95%% confidence\n",
            1 - (1 - (cen + hw))^4))
