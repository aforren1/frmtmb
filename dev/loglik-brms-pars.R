# Lane loglik: the Stan parameters each brms 2.23.0 program declares, for
# the residual covariance cases the probe uses.
source("C:/Users/adf44/source/r/frmtmb-wt-loglik/dev/loglik-brms-lib.R")
set.seed(1)
d <- do.call(rbind, lapply(1:5, function(g) {
  data.frame(g = g, t = seq_len(c(6, 5, 4, 6, 3)[g]))
}))
d$x <- rnorm(nrow(d))
d$y <- rnorm(nrow(d))
n <- nrow(d)
V <- 0.5^abs(outer(1:n, 1:n, "-"))
P <- matrix(0, n, n)
P[cbind(1:(n - 1), 2:n)] <- 1
P <- P + t(P)
dg <- 1 / sqrt(rowSums(P))
W <- P * outer(dg, dg)
forms <- list(
  ar = bf(y ~ x + ar(t, g, cov = TRUE)),
  ma = bf(y ~ x + ma(t, g, cov = TRUE)),
  arma = bf(y ~ x + arma(t, g, cov = TRUE)),
  arhet = bf(y ~ x + ar(t, g, cov = TRUE), sigma ~ x),
  cosy = bf(y ~ x + cosy(t, g)),
  unstr = bf(y ~ x + unstr(t, g)),
  fcor = bf(y ~ x + fcor(V)),
  lag = bf(y ~ x + sar(W)),
  err = bf(y ~ x + sar(W, type = "error")))
for (nm in names(forms)) {
  for (fam in c("gaussian", "student")) {
    code <- tryCatch(make_stancode(forms[[nm]], data = d,
                                   family = brmsfamily(fam),
                                   data2 = list(V = V, W = W)),
                     error = function(e) paste("ERROR", conditionMessage(e)))
    lines <- strsplit(code, "\n")[[1]]
    s <- grep("^parameters \\{", lines)
    e <- grep("^transformed parameters \\{", lines)
    cat("==", nm, fam, "\n")
    if (length(s)) cat(lines[(s + 1):(e - 2)], sep = "\n") else cat(code, "\n")
    lk <- grep("_lpdf\\(Y|target \\+=.*Y", lines, value = TRUE)
    cat("LIK:", lk, sep = "\n  ")
  }
}
