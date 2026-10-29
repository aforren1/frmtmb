# Short check: recompute the headline, the sensitivity row and the
# per-stratum rates from results.tsv alone, with my own Wilson/Wald.
r <- read.delim("dev/equiv-data/results.tsv", stringsAsFactors = FALSE)
st <- read.delim("dev/equiv-data/strata.tsv")
cp <- read.delim("dev/equiv-frame-out/compat.tsv"); cp <- cp[cp$in_core_only, ]
w <- read.delim("dev/equiv-data/weights-drawn.tsv")
cp$w <- w$weight[match(cp$feature_a, w$feature)] *
  w$weight[match(cp$feature_b, w$feature)]
W <- tapply(cp$w, cp$status, sum); W <- W / sum(W)
cat("W from frame:", round(W, 4), " strata.tsv:",
    round(st$W_share[match(names(W), st$status)], 4), "\n")
stopifnot(nrow(r) == 172, sum(r$draws) == 200)
agree <- r$verdict %in% c("equivalent", "refusal correct") |
  (r$verdict == "compat row wrong" & as.logical(r$agree) %in% TRUE)
cat("agree flag consistent with verdicts:", all(agree == as.logical(r$agree)),
    "\n")
wil <- function(x, n) { z <- 1.959964; p <- x / n; d <- 1 + z^2 / n
  c(p, (p + z^2 / (2 * n) - z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / d,
    (p + z^2 / (2 * n) + z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / d) }
pop <- function(ok) { p <- v <- 0
  for (s in names(W)) { k <- r$status == s; n <- sum(r$draws[k])
    ph <- sum(r$draws[k] * ok[k]) / n
    cat(sprintf("   %-12s %5.1f [%5.1f, %5.1f]\n", s, 100 * ph,
                100 * wil(ph * n, n)[2], 100 * wil(ph * n, n)[3]))
    p <- p + W[[s]] * ph; v <- v + W[[s]]^2 * ph * (1 - ph) / n }
  cat(sprintf("   population %.1f [%.1f, %.1f]\n", 100 * p,
              100 * (p - 1.96 * sqrt(v)), 100 * (p + 1.96 * sqrt(v)))) }
cat("agree:\n"); pop(agree)
weak <- c("p042", "p044", "p047", "p136", "p143")
cat("weak evidence out:\n"); pop(agree & !r$id %in% weak)
cat("no defect:\n"); pop(agree | r$verdict == "divergence by decision")
cat("\nverdict changes vs review expectations:\n")
print(r[r$id %in% c("p029", "p158", "p167", "p169", "p022", "p023", "p094",
                    "p109", "p070", "p087"),
        c("id", "verdict", "agree", "level")], row.names = FALSE)
cat("\nL4-only agree pairs:\n")
print(r[agree & r$top_level == "L4", c("id", "l4_route", "l4_z")],
      row.names = FALSE)
