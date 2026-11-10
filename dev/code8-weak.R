# Lane code8 (copy of dev/review-code8-weak.R), item 2: weakly curved but identified directions
# (Michaelis-Menten far below K, a biexponential) must keep their SEs,
# and a flat curved ridge (A * B) must lose exactly its parameters.
# Reference per fit: a Richardson Hessian of the objective at the stop
# (numDeriv::jacobian of the gradient); its smallest unit-diagonal
# eigenvalue and the SEs of its inverse.
#   RV_ARM=before|after Rscript dev/review-code8-weak.R <label>
lab <- commandArgs(TRUE)[1]
source("dev/review-code8-common.R")
mkdat <- function(model, seed) {
  set.seed(seed)
  if (model %in% c("mm_re", "mm")) {
    n <- 40
    d <- data.frame(g = factor(rep(1:8, each = 5)), x = runif(n, 0.1, 2))
    u <- rnorm(8, 0, 0.3)[d$g]
    d$y <- 1 + u + 10 * d$x / (5 + d$x) + rnorm(n, 0, 0.1)
  } else if (model %in% c("biexp", "biexp_re")) {
    n <- 40
    d <- data.frame(g = factor(rep(1:8, each = 5)), x = runif(n, 0, 10))
    u <- rnorm(8, 0, 0.3)[d$g]
    d$y <- (5 + if (model == "biexp_re") u else 0) * exp(-1 * d$x) +
      2 * exp(-0.1 * d$x) + rnorm(n, 0, 0.05)
  } else if (model == "prod_re") {
    n <- 60
    d <- data.frame(g = factor(rep(1:10, each = 6)), x = runif(n, 0, 3))
    u <- rnorm(10, 0, 0.5)[d$g]
    d$y <- 1 + u + 3 * exp(-0.7 * d$x) + rnorm(n, 0, 0.2)
  }
  d
}
forms <- list(
  mm = bf(y ~ c0 + Vm * x / (K + x), c0 ~ 1, Vm ~ 1, K ~ 1, nl = TRUE),
  mm_re = bf(y ~ c0 + Vm * x / (K + x), c0 ~ 1 + (1 | g), Vm ~ 1, K ~ 1,
             nl = TRUE),
  biexp = bf(y ~ A1 * exp(-k1 * x) + A2 * exp(-k2 * x), A1 ~ 1, k1 ~ 1,
             A2 ~ 1, k2 ~ 1, nl = TRUE),
  biexp_re = bf(y ~ A1 * exp(-k1 * x) + A2 * exp(-k2 * x), A1 ~ 1 + (1 | g),
                k1 ~ 1, A2 ~ 1, k2 ~ 1, nl = TRUE),
  prod_re = bf(y ~ c0 + A * B * exp(-k * x), c0 ~ 1 + (1 | g), A ~ 1, B ~ 1,
               k ~ 1, nl = TRUE))
starts <- list(
  mm = c(c0 = 1, Vm = 10, K = 5), mm_re = c(c0 = 1, Vm = 10, K = 5),
  biexp = c(A1 = 5, k1 = 1, A2 = 2, k2 = 0.1),
  biexp_re = c(A1 = 5, k1 = 1, A2 = 2, k2 = 0.1),
  prod_re = c(c0 = 1, A = 1.5, B = 2, k = 0.7))
mkprior <- function(st) {
  do.call(c, lapply(names(st), function(nm) {
    set_prior(sprintf("normal(%g, 1000)", st[[nm]]), nlpar = nm)
  }))
}
res <- list()
for (model in names(forms)) for (seed in 1:3) {
  d <- mkdat(model, seed)
  y0 <- d$y
  for (j in c(0, -3, 5, 9, -11, 13, -17, 19)) {
    d$y <- y0 * (1 + 2e-15 * j)
    r <- rfit(forms[[model]], data = d, prior = mkprior(starts[[model]]))
    if (r$err) { cat(model, seed, j, "error", conditionMessage(r$f), "\n"); next }
    f <- r$f
    s <- se_of(f)
    p <- f$opt$par
    J <- tryCatch(numDeriv::jacobian(function(q) drop(f$obj$gr(q)), p,
                                     method = "Richardson"),
                  error = function(e) NULL)
    mev <- NA; rse <- rep(NA, length(p)); ok <- NA
    if (!is.null(J) && all(is.finite(J))) {
      H <- (J + t(J)) / 2
      D <- sqrt(abs(diag(H)))
      ev <- eigen(H / outer(D, D), symmetric = TRUE, only.values = TRUE)$values
      mev <- min(ev) / max(ev)
      Vr <- tryCatch(solve(H), error = function(e) NULL)
      if (!is.null(Vr)) rse <- sqrt(pmax(diag(Vr), NA))
    }
    rr <- rse / s$se
    res[[length(res) + 1L]] <- data.frame(model = model, seed = seed, j = j,
      code = f$opt$convergence, obj = sprintf("%.7f", f$opt$objective),
      ref_minev = signif(mev, 3),
      lost = paste(sort(names(s$lost)), collapse = ","),
      se_ratio_range = paste(signif(range(rr, na.rm = TRUE), 4),
                             collapse = ".."),
      conds = substr(cond_txt(r), 1, 70))
  }
}
X <- do.call(rbind, res)
options(width = 220)
print(X, right = FALSE)
write.table(X, sprintf("dev/code8-out/weak-%s.tsv", lab), sep = "\t",
            quote = FALSE, row.names = FALSE)
cat("DONE\n")
