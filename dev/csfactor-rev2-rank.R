# Is the REFINEMENT refusal right? The reference build cannot answer it,
# because there cs(mr) was a single integer column, so a different model.
# Check the rank condition the refusal rests on, directly, and check the
# accepted neighbours fail it.
#   Rscript dev/csfactor-rev2-rank.R <lib>
args <- commandArgs(TRUE)
.libPaths(c(args[[1L]], "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("== frmtmb", as.character(packageVersion("frmtmb")), "\n")
set.seed(3151)
n <- 600
d <- data.frame(x = rnorm(n), z = rnorm(n))
d$zb <- rbinom(n, 1L, 0.5)
mcode <- sample(0:3, n, TRUE)
d$m <- factor(mcode, levels = 0:3, ordered = TRUE)
d$mc <- factor(c("A", "A", "B", "B")[mcode + 1L])
d$mr <- factor(paste0(mcode, sample(c("p", "q"), n, TRUE)))
m2code <- pmin(3L, pmax(0L, mcode + sample(c(-1L, 0L, 1L), n, TRUE)))
d$m2 <- factor(m2code, levels = 0:3, ordered = TRUE)
d$f <- factor(sample(c("a", "b", "c"), n, TRUE))
eff <- c(0, 0.15, 0.30, 1.30)[mcode + 1L]
eta <- 0.5 * d$x + eff
p1 <- plogis(-0.7 - eta); p2 <- plogis(0.9 - eta)
P <- cbind(p1, p2 - p1, 1 - p2)
d$yo <- apply(P, 1L, function(p) sample.int(3L, 1L, prob = pmax(p, 1e-9)))

rnk <- function(A) {
  s <- sqrt(colSums(A^2)); s[!(s > 0)] <- 1
  qr(sweep(A, 2L, s, "/"))$rank
}
Dmo <- outer(mcode, 1:3, "==") * 1
colnames(Dmo) <- paste0("mom.", 1:3)

# every cs() design the accept/refuse pairs use, built the way the frame
# builds it, then the two rank conditions read off directly
cases <- list(
  "cs(m)"  = ~ m, "cs(mr) refinement" = ~ mr,
  "cs(mc) coarsening" = ~ mc, "cs(m2) correlated" = ~ m2,
  "cs(f) unrelated" = ~ f, "cs(x)" = ~ x
)
cat(sprintf("%-22s %6s %6s %8s %8s  %s\n", "cs design", "ncolZ",
            "rk(1Z)", "rk(1ZD)", "rk(1D)", "verdict"))
for (nmx in names(cases)) {
  Z <- stats::model.matrix(cases[[nmx]], d)
  Z <- Z[, setdiff(colnames(Z), "(Intercept)"), drop = FALSE]
  M <- matrix(1, n, 1L)
  a <- rnk(cbind(M, Z)); b <- rnk(cbind(M, Z, Dmo)); c0 <- rnk(cbind(M, Dmo))
  cat(sprintf("%-22s %6d %6d %8d %8d  %s\n", nmx, ncol(Z), a, b, c0,
              if (b == a && c0 > 1L) "REFUSE" else "accept"))
}
# the interaction basis, z and zb
for (tg in c("z", "zb")) {
  Dm <- Dmo * d[[tg]]
  Z <- stats::model.matrix(~ m, d)
  Z <- Z[, -1L, drop = FALSE]
  M <- matrix(1, n, 1L)
  cat(sprintf("%-22s %6d %6d %8d %8d  %s\n", paste0("cs(m) vs mo(m):", tg),
              ncol(Z), rnk(cbind(M, Z)), rnk(cbind(M, Z, Dm)),
              rnk(cbind(M, Dm)),
              if (rnk(cbind(M, Z, Dm)) == rnk(cbind(M, Z))) "REFUSE"
              else "accept"))
}
cat("\n# the refinement, stated as an identity: every indicator of m is a\n",
    "# sum of indicators of mr, so col(Dmo) is inside col(cbind(1, Zmr)).\n",
    sep = "")
Zr <- stats::model.matrix(~ mr, d)[, -1L, drop = FALSE]
resid_norm <- max(vapply(1:3, function(j) {
  r <- stats::lsfit(cbind(1, Zr), Dmo[, j], intercept = FALSE)$residuals
  sqrt(sum(r^2))
}, 0))
cat("max || residual of a Dmo column on cbind(1, Zmr) || = ",
    sprintf("%.6e", resid_norm), "\n", sep = "")
Zc <- stats::model.matrix(~ mc, d)[, -1L, drop = FALSE]
resid_c <- min(vapply(1:3, function(j) {
  r <- stats::lsfit(cbind(1, Zc), Dmo[, j], intercept = FALSE)$residuals
  sqrt(sum(r^2))
}, 0))
cat("min || residual of a Dmo column on cbind(1, Zmc) || (coarsening) = ",
    sprintf("%.6g", resid_c), "\n", sep = "")
cat("done\n")
