# R3: what the beta-binomial density costs on the tape, by form. Whole fits,
# minimum of repeats, on one build. Candidates are custom families so every
# arm runs through the same frm(). Run as: Rscript dev/nlfix-bbcost.R <reps>
lib <- "C:/Users/adf44/source/r/wt-nlfix-lib2"
.libPaths(c(lib, "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(getNamespaceInfo("frmtmb", "path")), "\n")
ns <- asNamespace("frmtmb")
reps <- as.integer(commandArgs(TRUE)[1])
if (is.na(reps)) reps <- 3L

shift_full <- ns$lgamma_shift_diff

# the finite product for rows with a short shift, the series for the rest
shift_hybrid <- function(M) {
  function(a, s) {
    "[<-" <- RTMB::ADoverload("[<-")
    if (length(a) == 1L) a <- a[rep_len(1L, length(s))]
    long <- which(s > M)
    acc <- rep(0 * a[1L], length(s))
    acc <- a[1L] * 0
    acc <- rep(0, length(s))
    short <- s
    short[long] <- 0
    top <- max(short)
    for (j in seq_len(top) - 1L) {
      i <- which(short > j)
      acc[i] <- acc[i] + log(a[i] + j)
    }
    if (length(long)) acc[long] <- ns$lgamma_shift_diff(a[long], s[long])
    acc
  }
}

fam <- function(tag, shift) {
  frmtmb_family(tag, dpars = c("mu", "phi"),
    links = list(mu = "logit", phi = "log"), type = "discrete",
    required_aterms = "trials",
    lpdf = function(y, dpars, aterms = NULL, extra = NULL) {
      size <- ns$aterm_trials(aterms, tag)
      a <- dpars$mu * dpars$phi
      b <- (1 - dpars$mu) * dpars$phi
      if (is.null(shift)) {
        return(RTMBdist::dbetabinom(y, size, a, b, log = TRUE))
      }
      lgamma(size + 1) - lgamma(y + 1) - lgamma(size - y + 1) +
        shift(a, y) + shift(b, size - y) - shift(dpars$phi, size)
    },
    init_dpars = list(mu = function(y, aterms) 0.5,
                      phi = function(y, aterms) 5))
}

arms <- list(raw = NULL, full = shift_full, h6 = shift_hybrid(6L),
             h12 = shift_hybrid(12L), h24 = shift_hybrid(24L),
             h48 = shift_hybrid(48L))

mk <- function(n, ntr, phi = 8) {
  set.seed(77)
  d <- data.frame(x = rnorm(n), ntr = as.integer(ntr))
  p <- plogis(0.3 + 0.5 * d$x)
  d$y <- rbinom(n, d$ntr, rbeta(n, p * phi, (1 - p) * phi))
  d
}

for (cfg in list(list(n = 3e4, ntr = 10), list(n = 3e4, ntr = 30),
                 list(n = 3e4, ntr = 200))) {
  d <- mk(cfg$n, cfg$ntr)
  cat(sprintf("\n--- n = %g, trials = %g\n", cfg$n, cfg$ntr))
  for (nm in names(arms)) {
    fm <- fam(paste0("bb_", nm), arms[[nm]])
    tt <- numeric(reps)
    f <- NULL
    for (i in seq_len(reps)) {
      gc(FALSE)
      tt[i] <- system.time(
        f <- suppressMessages(suppressWarnings(
          frm(bf(y | trials(ntr) ~ x), family = fm, data = d))))[["elapsed"]]
    }
    p <- f$opt$par
    te <- min(replicate(5L, system.time(for (k in 1:20) f$obj$fn(p))[["elapsed"]])) / 20
    tg <- min(replicate(5L, system.time(for (k in 1:20) f$obj$gr(p))[["elapsed"]])) / 20
    cat(sprintf("%-5s fit %7.2f s  fn %8.5f s  gr %8.5f s  ll %.9f code %d\n",
                nm, min(tt), te, tg, as.numeric(logLik(f)), f$opt$convergence))
  }
}

# second pass: the same finite product, row sums through one sparse product,
# so the tape build makes no per-j copy of the rows
shift_sparse <- function(M) {
  function(a, s) {
    "[<-" <- RTMB::ADoverload("[<-")
    n <- max(length(a), length(s))
    s <- rep_len(s, n)
    if (length(a) == 1L && n > 1L) a <- a[rep_len(1L, n)]
    long <- which(s > M)
    short <- s
    short[long] <- 0
    r <- rep.int(seq_len(n), short)
    acc <- if (length(r)) {
      A <- Matrix::sparseMatrix(i = r, j = seq_along(r), x = 1,
                                dims = c(n, length(r)))
      as.vector(A %*% log(a[r] + (sequence(short) - 1)))
    } else rep(0, n)
    if (length(long)) acc[long] <- ns$lgamma_shift_diff(a[long], s[long])
    acc
  }
}
arms2 <- list(raw = NULL, h24 = shift_hybrid(24L), s12 = shift_sparse(12L),
              s24 = shift_sparse(24L), s48 = shift_sparse(48L),
              s96 = shift_sparse(96L))
for (cfg in list(list(n = 3e4, ntr = 10), list(n = 3e4, ntr = 30),
                 list(n = 3e4, ntr = 200))) {
  d <- mk(cfg$n, cfg$ntr)
  cat(sprintf("\n--- sparse pass: n = %g, trials = %g\n", cfg$n, cfg$ntr))
  for (nm in names(arms2)) {
    fm <- fam(paste0("bb2_", nm), arms2[[nm]])
    tt <- numeric(reps)
    f <- NULL
    for (i in seq_len(reps)) {
      gc(FALSE)
      tt[i] <- system.time(
        f <- suppressMessages(suppressWarnings(
          frm(bf(y | trials(ntr) ~ x), family = fm, data = d))))[["elapsed"]]
    }
    p <- f$opt$par
    tg <- min(replicate(5L, system.time(for (k in 1:20) f$obj$gr(p))[["elapsed"]])) / 20
    cat(sprintf("%-5s fit %7.2f s  gr %8.5f s  ll %.9f code %d\n",
                nm, min(tt), tg, as.numeric(logLik(f)), f$opt$convergence))
  }
}
