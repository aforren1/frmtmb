source("C:/Users/adf44/source/r/frmtmb-wt-predfix/dev/predfix-review2/r2-helpers.R")
suppressPackageStartupMessages(library(frmtmb.latent))
cat("frmtmb.latent from", find.package("frmtmb.latent"), "\n")
# Minor 8 on my own designs: hmm with a covariate at 5e-4 and a
# 40 x 30 panel, and lca(K = 3) at 5e-4. Seed 909 and 910.
q <- function(e) tryCatch(suppressWarnings(suppressMessages(e)),
                          error = function(x) x)
ll <- function(f) if (inherits(f, "error")) NA else as.numeric(logLik(f))
cmp3 <- function(lab, mk, d, col, s = 5e-4) {
  ref <- q(mk(d, FALSE))
  d[[col]] <- d[[col]] * s
  fd <- q(mk(d, NULL)); ff <- q(mk(d, FALSE))
  cat(sprintf("%-34s ref %.6f | default %.6f engaged=%s code %s | FALSE %.6f code %s\n",
              lab, ll(ref), ll(fd),
              if (inherits(fd, "error")) conditionMessage(fd) else !is.null(fd$par_units),
              if (inherits(fd, "error")) NA else fd$opt$convergence,
              ll(ff), if (inherits(ff, "error")) NA else ff$opt$convergence))
  if (!inherits(fd, "error") && !inherits(ref, "error"))
    cat(sprintf("   slope default / (ref * 1/s): %.8f\n",
                fixef(fd)[grep(col, rownames(fixef(fd)))[1], "Estimate"] /
                  (fixef(ref)[grep(col, rownames(fixef(ref)))[1], "Estimate"] / s)))
}
for (seed in 909:910) {
  set.seed(seed)
  G <- matrix(c(0.85, 0.15, 0.2, 0.8), 2, 2, byrow = TRUE)
  dd <- do.call(rbind, lapply(1:40, function(id) {
    s <- integer(30); s[1L] <- sample.int(2, 1L)
    for (t in 2:30) s[t] <- sample.int(2, 1L, prob = G[s[t - 1L], ])
    x <- rnorm(30)
    data.frame(id = id, t = 1:30, x = x,
               y = rnorm(30, c(-1, 2)[s] + 0.8 * x, 0.7))
  }))
  cmp3(paste("hmm K=2 y ~ x, seed", seed),
       function(d, a) frm(bf(y ~ x), family = hmm(K = 2, gaussian(), time = t,
                                                  group = id),
                          data = d, control = frmtmb_control(autoscale = a)),
       dd, "x")
  set.seed(seed)
  n <- 600
  x <- rnorm(n)
  eta2 <- -0.3 + 1.0 * x; eta3 <- 0.2 - 0.9 * x
  pc <- cbind(1, exp(eta2), exp(eta3)); pc <- pc / rowSums(pc)
  cl <- apply(pc, 1, function(p) sample.int(3, 1, prob = p))
  pr <- rbind(c(0.9, 0.85, 0.2, 0.75, 0.8), c(0.15, 0.1, 0.9, 0.2, 0.7),
              c(0.5, 0.9, 0.9, 0.9, 0.1))
  Y <- matrix(0L, n, 5)
  for (j in 1:5) Y[, j] <- 1L + rbinom(n, 1, pr[cl, j])
  dl <- data.frame(x = x); dl$Y <- Y
  cmp3(paste("lca K=3 Y ~ x, seed", seed),
       function(d, a) frm(bf(Y ~ x), family = lca(K = 3), data = d,
                          control = frmtmb_control(autoscale = a)), dl, "x")
}
