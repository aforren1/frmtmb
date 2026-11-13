# Short check of lane arspat's fix round, part B: sar() with complex
# eigenvalues, the sign of det(I - rho W) over the whole interval, and
# log|det| against determinant(). Part C: se() beside ar(cov = FALSE) on
# gaussian, student and skew_normal (run on the lane and on rellib-r20).
source("dev/reviews/arspat-review-common.R")
part <- commandArgs(trailingOnly = TRUE)[1]
if (identical(part, "B")) {
  suppressPackageStartupMessages(library(frmtmb.spatial))
  dirW <- function(n, k, seed) {
    set.seed(seed)
    xy <- cbind(runif(n), runif(n))
    D <- as.matrix(dist(xy)); diag(D) <- Inf
    K <- t(apply(D, 1, function(r) as.numeric(rank(r) <= k)))
    K / rowSums(K)
  }
  # the lane's construction: a directed 3-cycle (complex real part -0.5)
  # beside an edge scaled by 0.2 (real eigenvalues +-0.2)
  cyc <- matrix(0, 3, 3); cyc[cbind(1:3, c(2, 3, 1))] <- 1
  W5 <- matrix(0, 5, 5); W5[1:3, 1:3] <- cyc; W5[4, 5] <- W5[5, 4] <- 0.2
  Ws <- list(lane5 = W5)
  for (k in 1:4) for (s in 1:3) Ws[[sprintf("knn%d_s%d", k, s)]] <- dirW(40, k, s)
  # unstandardized directed graph, and a random nonnegative matrix
  set.seed(9); A <- matrix(rbinom(900, 1, 0.1), 30); diag(A) <- 0
  Ws$binary30 <- A
  set.seed(10); B <- matrix(runif(400), 20) * (matrix(runif(400), 20) < 0.2)
  diag(B) <- 0; Ws$rand20 <- B
  for (nm in names(Ws)) {
    W <- Ws[[nm]]; n <- nrow(W)
    aux <- tryCatch(frmtmb.spatial:::sar_build(list(M = quote(W), type = "lag"),
                                               seq_len(n), n, function(e, w) W),
                    error = function(e) conditionMessage(e))
    if (is.character(aux)) { cat(nm, "REFUSED:", aux, "\n"); next }
    ev <- eigen(W, only.values = TRUE)$values
    rr <- seq(aux$lo, aux$hi, length.out = 2003)[-c(1, 2003)]
    sg <- vapply(rr, function(r) determinant(diag(n) - r * W)$sign, 0)
    df <- vapply(rr, function(r) abs(as.numeric(frmtmb.spatial:::sar_logdet(r, aux)) -
                   as.numeric(determinant(diag(n) - r * W)$modulus)), 0)
    ends <- c(aux$lo * (1 - 1e-7), aux$hi * (1 - 1e-7))
    dfe <- vapply(ends, function(r) abs(as.numeric(frmtmb.spatial:::sar_logdet(r, aux)) -
                    as.numeric(determinant(diag(n) - r * W)$modulus)), 0)
    sge <- vapply(ends, function(r) determinant(diag(n) - r * W)$sign, 0)
    # just outside: the first real eigenvalue's root is the boundary
    out <- c(aux$lo * (1 + 1e-6), aux$hi * (1 + 1e-6))
    sgo <- vapply(out, function(r) determinant(diag(n) - r * W)$sign, 0)
    cat(sprintf("%-9s complex %2d interval (%.5f, %.5f) min sign %d, max |diff| %.2e interior, ends %.2e sign %d %d, just outside sign %d %d\n",
                nm, sum(abs(Im(ev)) > 1e-8 * max(Mod(ev))), aux$lo, aux$hi,
                min(sg), max(df), max(dfe), sge[1], sge[2], sgo[1], sgo[2]))
  }
}
if (identical(part, "C")) {
  set.seed(7)
  d <- expand.grid(t = 1:6, g = factor(1:20))
  d$x <- rnorm(nrow(d)); d$s <- runif(nrow(d), 0.2, 0.6)
  d$y <- 1 + 0.6 * d$x + rnorm(nrow(d))
  for (fam in list(gaussian(), student(), skew_normal())) {
    for (tm in c("ar(t, g)", "ar(t, g, p = 2)", "ma(t, g)", "arma(t, g)")) {
      chk(paste(fam$family, tm, "se(sigma = TRUE)"), {
        f <- frm(bf(stats::as.formula(paste("y | se(s, sigma = TRUE) ~ x +", tm))),
                 family = fam, data = d)
        cat("logLik", format(as.numeric(logLik(f)), digits = 10), "\n")
        NULL
      })
    }
    chk(paste(fam$family, "ar(t, g) se(sigma = FALSE)"), {
      f <- frm(y | se(s) ~ x + ar(t, g), family = fam, data = d)
      cat("logLik", format(as.numeric(logLik(f)), digits = 10), "\n")
      NULL
    })
  }
  if (nzchar(Sys.getenv("ARSPAT_BRMS"))) {
    for (fam in c("gaussian", "student", "skew_normal")) {
      chk(paste("brms", fam, "se + ar(cov = FALSE)"), {
        brms::make_stancode(y | se(s, sigma = TRUE) ~ x + ar(t, g),
                            family = fam, data = d)
        "fits"
      })
    }
    chk("brms skew_normal ar(cov = FALSE) without se: structure", {
      sc <- brms::make_stancode(y ~ x + ar(t, g), family = brms::skew_normal(),
                                data = d)
      cat(grep("err|Err|ar\\[|sderr|omega\\[n\\] =|mu\\[n\\] \\+=",
               strsplit(sc, "\n")[[1]], value = TRUE), sep = "\n")
      NULL
    })
  }
}
