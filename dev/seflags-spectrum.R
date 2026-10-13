# Lane seflags, item 3: the unit-diagonal spectrum the SE check reads on
# the curved ridges (c0k, a + exp(b), a * b, the same with (1 | g)) and
# on mo() with an unobserved category, and the noise of each direction:
# where do the ridge's eigenvalues sit against the identified ones?
#   Rscript dev/seflags-spectrum.R <lib> [seeds]
a <- commandArgs(trailingOnly = TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
seeds <- if (length(a) > 1) eval(parse(text = a[2])) else 1:6
cat("lib", find.package("frmtmb"), "\n")
quiet <- function(e) suppressMessages(suppressWarnings(e))
spec <- function(f) {
  p <- f$opt$par
  if (!length(f$obj$env$random)) {
    H <- f$obj$he(p); E <- abs(H - t(H)) / 2; H <- (H + t(H)) / 2
  } else {
    h <- f$cache$hessian_fixed %||% ns$fit_outer_hessian(f)
    H <- h$H; E <- h$E
  }
  D <- sqrt(abs(diag(H)))
  e <- eigen(H / outer(D, D), symmetric = TRUE)
  Es <- E / outer(D, D)
  nk <- sqrt(colSums((Es %*% e$vectors)^2))
  k <- length(e$values)
  g <- f$obj$gr(p)
  list(rel = e$values / max(abs(e$values)), ev = e$values, nk = nk,
       gmax = max(abs(g)), vmin = e$vectors[, k])
}
row <- function(des, s, f) {
  sp <- tryCatch(spec(f), error = function(e) NULL)
  if (is.null(sp)) { cat(des, s, "no spectrum\n"); return() }
  k <- length(sp$ev)
  lost <- quiet(ns$sdr_of(f))$se_lost
  cat(sprintf("%-8s seed %2d code %d | min ev %.3g rel %.3g noise %.2g | next %.3g | grad %.2g | lost %s\n",
              des, s, f$opt$convergence, sp$ev[k], sp$rel[k], sp$nk[k],
              sp$ev[k - 1], sp$gmax,
              if (length(lost)) paste(names(lost), collapse = ",") else "none"))
}
for (s in seeds) {
  set.seed(s)
  n <- 60
  dd <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:6, 10)))
  dd$y <- 3 + 0.5 * dd$x + 0.3 * dd$z + rnorm(6, 0, 0.5)[dd$g] +
    rnorm(n, 0, 0.4)
  row("aexpb", s, quiet(frm(bf(y ~ a + exp(b), a ~ 1 + x, b ~ 1, nl = TRUE),
                            data = dd, start = list(beta = c(1, 0, -1)))))
  row("ab", s, quiet(frm(bf(y ~ a * b, a ~ 1, b ~ 1, nl = TRUE), data = dd,
                         start = list(beta = c(1, 1)))))
  row("aexpbre", s, quiet(frm(bf(y ~ a + exp(b), a ~ 1 + x + (1 | g), b ~ 1,
                                 nl = TRUE), data = dd,
                              start = list(beta = c(1, 0, -1)))))
  set.seed(s)
  G <- 8
  dn <- data.frame(g = factor(rep(seq_len(G), each = 10)))
  dn$x <- rnorm(nrow(dn))
  u <- rnorm(G, 0, 0.6)[dn$g]
  dn$yn <- 1 + 0.3 * dn$x + exp(0.3 + u)^0.8 + rnorm(nrow(dn), 0, 0.3)
  f <- tryCatch(quiet(frm(bf(yn ~ c0 + exp(a)^k, c0 ~ 1 + x,
                             a ~ 1 + (1 | g), k ~ 1, nl = TRUE), data = dn)),
                error = function(e) NULL)
  if (!is.null(f)) row("c0k", s, f)
  set.seed(s)
  inc <- sample(c(0L, 1L, 3L), 100, TRUE)
  ls <- c(30, 60, 70, 75)[inc + 1L] + rnorm(100, sd = 7)
  dm <- data.frame(inc, ls, age = rnorm(100, mean = 40, sd = 10))
  row("segap", s, quiet(frm(ls ~ mo(inc), data = dm)))
  row("segap+a", s, quiet(frm(ls ~ mo(inc) + age, data = dm)))
}
