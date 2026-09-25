# Item 3.4: a high-precision reference for the Wiener response-time
# distribution function, both boundaries together, and the sweep that
# sets the blend center.
#
# The reference is the LARGE-time eigenfunction series evaluated in
# Rmpfr at 1000 bits, summed until the next term is below 2^-1000 of the
# running total, so it has no truncation to argue about. It is checked
# against the SMALL-time image series at the same precision, which is a
# separate derivation (method of images against separation of
# variables) and shares no algebra with it. Where the two agree to many
# more digits than a double holds, either is a reference.
#
# Output: dev/phase3b-log/cdf-reference.csv (the grid and the 1000-bit
# values, as 17-digit strings) and dev/phase3b-log/cdf-reference.txt.
#
# Run: Rscript dev/phase3b-cdf-reference.R   (about two minutes)

.libPaths(c("C:/Users/adf44/source/r/phase3b-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(Rmpfr))
PREC <- 1000L
mp <- function(x) mpfr(x, PREC)
PI <- Const("pi", PREC)

# P(T > t) = sum over both boundaries of the remaining mass, large time
ref_surv_large <- function(t, v, a, w) {
  t <- mp(t); v <- mp(v); a <- mp(a); w <- mp(w)
  va2 <- v * v * a * a
  el <- -v * a * w
  eu <- v * a * (1 - w)
  s <- mp(0)
  k <- 1L
  repeat {
    den <- va2 + k * k * PI * PI
    ex <- -den * t / (2 * a * a)
    term <- (k * sin(k * PI * w) / den) *
      (exp(el + ex) + (-1)^(k + 1) * exp(eu + ex))
    s <- s + term
    if (k > 5L && abs(term) < abs(s) * mp(2)^(-1000)) break
    k <- k + 1L
    if (k > 20000L) stop("large series did not converge")
  }
  2 * PI * s
}

pn <- function(x) erfc(-x / sqrt(mp(2))) / 2

# lower-boundary defective CDF, method of images, K images each side
ref_lower_small <- function(t, v, a, w, K = 80L) {
  t <- mp(t); v <- mp(v); a <- mp(a); w <- mp(w)
  st <- sqrt(t)
  out <- mp(0)
  for (j in 0:K) {
    c <- a * (w + 2 * j)
    out <- out + exp(-2 * v * a * (w + j)) * pn((v * t - c) / st) +
      exp(2 * v * a * j) * pn(-(v * t + c) / st)
  }
  for (j in (-K):(-1)) {
    c <- a * (w + 2 * j)
    out <- out - exp(-2 * v * a * (w + j)) * pn(-(v * t - c) / st) -
      exp(2 * v * a * j) * pn((v * t + c) / st)
  }
  out
}
ref_cdf_small <- function(t, v, a, w) {
  ref_lower_small(t, v, a, w) + ref_lower_small(t, -v, a, 1 - w)
}

# the grid test-density.R pins the density on, plus a tail extension
gr_d <- expand.grid(t = c(0.01, 0.05, 0.15, 0.4, 0.8, 1.5, 2.5, 4, 8),
                    a = c(0.4, 0.8, 1.4, 2.2, 3.5),
                    w = c(0.2, 0.35, 0.5, 0.7, 0.85),
                    v = c(-3, -1.5, 0, 1.5, 3))
gr_d$grid <- "density"
gr_x <- expand.grid(t = c(0.02, 0.3, 1, 3, 6, 15),
                    a = c(0.8, 1.4, 2.5, 4),
                    w = c(0.1, 0.9),
                    v = c(-5, -0.5, 0.5, 5))
gr_x$grid <- "tail"
gr <- rbind(gr_d, gr_x)

t0 <- proc.time()
S <- character(nrow(gr)); F <- character(nrow(gr))
agree <- rep(NA_real_, nrow(gr))
for (i in seq_len(nrow(gr))) {
  s <- ref_surv_large(gr$t[i], gr$v[i], gr$a[i], gr$w[i])
  S[i] <- formatMpfr(s, digits = 25)
  F[i] <- formatMpfr(1 - s, digits = 25)
  u <- gr$t[i] / gr$a[i]^2
  # the image series converges fast only for small u; check it where it
  # does, which is the half of the grid where F is small and the large
  # series is doing its hardest subtraction
  if (u < 2) {
    f <- ref_cdf_small(gr$t[i], gr$v[i], gr$a[i], gr$w[i])
    agree[i] <- asNumeric(abs(f - (1 - s)) / (1 - s))
  }
}
el <- (proc.time() - t0)[["elapsed"]]
gr$F_ref <- F
gr$S_ref <- S
dir.create("dev/phase3b-log", showWarnings = FALSE)
utils::write.csv(gr, "dev/phase3b-log/cdf-reference.csv", row.names = FALSE)

sink("dev/phase3b-log/cdf-reference.txt")
cat("1000-bit reference, ", nrow(gr), " rows, ", round(el, 1), " s\n",
    sep = "")
cat("rows where the image series was checked: ", sum(!is.na(agree)), "\n")
cat("worst relative disagreement between the two 1000-bit series: ",
    format(max(agree, na.rm = TRUE), digits = 3), "\n")
cat("smallest F in the grid: ",
    format(min(as.numeric(F)), digits = 3), "\n")
cat("smallest S in the grid: ",
    format(min(as.numeric(S)), digits = 3), "\n")
sink()
cat(readLines("dev/phase3b-log/cdf-reference.txt"), sep = "\n")
