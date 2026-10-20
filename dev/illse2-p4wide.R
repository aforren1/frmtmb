# Lane illse2, item 4: fits with more than 100 outer parameters. A raw
# polynomial of degree 7 on [1, 2] beside K standard normal covariates,
# gaussian, n rows, so the fit has 8 + K + 1 outer parameters. The
# reference is lm() (its QR), on the ML scale. CAP sets the most outer
# parameters the design's coordinates are taken for (se_reparam_np).
#   ILLSE2_ARM=lane CAP=400 Rscript dev/illse2-p4wide.R <K> <n> <seeds> [re]
source("C:/Users/adf44/source/r/frmtmb-wt-illse2/dev/illse2-common.R")
cap_np <- as.integer(Sys.getenv("CAP", "0"))
if (cap_np > 0L) assignInNamespace("se_reparam_np", cap_np, ns)
a <- commandArgs(TRUE)
K <- as.integer(a[1])
n <- as.numeric(a[2])
seeds <- eval(parse(text = a[3]))
re <- length(a) > 3 && as.integer(a[4]) == 1L
cat(sprintf("cap %d | K %d | n %g | re %s\n", ns$se_reparam_np, K, n, re))
for (seed in seeds) {
  set.seed(seed)
  d <- data.frame(x = stats::runif(n, 1, 2),
                  g = factor(sample(1:50, n, TRUE)))
  for (k in 1:7) d[[paste0("x", k)]] <- d$x^k
  Z <- matrix(stats::rnorm(n * K), n, K)
  colnames(Z) <- paste0("z", seq_len(K))
  d <- cbind(d, Z)
  d$y <- sin(3 * d$x) + 0.05 * rowSums(Z[, 1:10]) +
    (if (re) stats::rnorm(50, 0, 0.3)[d$g] else 0) +
    stats::rnorm(n, 0, 0.5)
  rhs <- paste(c(paste0("x", 1:7), colnames(Z)), collapse = " + ")
  fo <- stats::as.formula(paste("y ~", rhs, if (re) "+ (1 | g)"))
  pt <- proc.time()
  r <- cap(frm(fo, data = d))
  tt <- proc.time() - pt
  f <- r$v
  if (is_err(f)) { cat("seed", seed, "ERROR", f, "\n"); next }
  lost <- ns$sdr_of(f)$se_lost
  fe <- suppressWarnings(fixef(f))
  out <- if (!re) {
    l <- stats::lm(stats::as.formula(paste("y ~", rhs)), data = d)
    k <- length(stats::coef(l))
    se_ref <- sqrt(diag(stats::vcov(l)) * (n - k) / n)
    sprintf("logLik - lm %+.1e | SE / lm %s",
            as.numeric(logLik(f)) - as.numeric(logLik(l)),
            ratio_str(unname(fe[sub("^[(]Intercept[)]$", "Intercept",
                                    names(se_ref)), "Est.Error"] / se_ref)))
  } else {
    sprintf("SE finite %d of %d", sum(is.finite(fe[, "Est.Error"])),
            nrow(fe))
  }
  cat(sprintf(paste0("seed %d | outer %d | %.1f s (cpu %.1f) | code %d | ",
                     "lost %s | root %s | %s | warn %s\n"),
              seed, length(f$opt$par), tt[["elapsed"]],
              tt[["user.self"]] + tt[["sys.self"]], f$opt$convergence,
              if (length(lost)) paste(names(table(lost)), table(lost),
                                      collapse = ",") else "none",
              !is.null(ns$sdr_of(f)$cov_fixed_root), out, wtag(r$w)))
}
