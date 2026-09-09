# An INDEPENDENT high-precision reference for frm_lincmt(), written so
# that it does not share a line of code with the lane's own
# dev/lincmt/lincmt-mpfr.R. Two arms:
#
#   expm_ss()  scaling and squaring, 300 bits, 60 Taylor terms, s + 8
#              squarings. Never forms an eigenvalue.
#   eig_pf()   the TEXTBOOK partial-fraction form, with the eigenvalues
#              found by bisection and Newton at 300 bits. It DOES form
#              eigenvalues and it DOES divide by their differences,
#              which at 300 bits still leaves about 80 digits where the
#              double-precision spelling has none.
#
# Agreement between the two is what makes either usable as a reference.
# Measured: 2.46e-86 over 36 random points (dev/rev-lincmt-cancel.R,
# seed 4041).
suppressPackageStartupMessages(library(Rmpfr))
PB <- 300L
mp <- function(x) mpfr(x, PB)

expm_ss <- function(A) {
  n <- nrow(A)
  nrm <- max(as.numeric(apply(abs(A), 1, sum)))
  s <- max(0L, as.integer(ceiling(log2(max(nrm, 1e-300)))) + 8L)
  B <- A / mp(2)^s
  I <- mpfrArray(0, PB, c(n, n))
  for (i in seq_len(n)) I[i, i] <- mp(1)
  T <- I
  P <- I
  for (k in 1:60) {
    P <- P %*% B / mp(k)
    T <- T + P
  }
  for (k in seq_len(s)) T <- T %*% T
  T
}

# rate matrix, depot first when present; dA/dt = M A
mk_M <- function(ncmt, depot, p) {
  n <- ncmt + as.integer(depot)
  M <- mpfrArray(0, PB, c(n, n))
  ic <- if (depot) 2L else 1L
  d <- -mp(p[["ke"]])
  if (depot) {
    M[1L, 1L] <- -mp(p[["ka"]])
    M[ic, 1L] <- mp(p[["ka"]])
  }
  if (ncmt >= 2L) {
    d <- d - mp(p[["k12"]])
    M[ic + 1L, ic] <- mp(p[["k12"]])
    M[ic, ic + 1L] <- mp(p[["k21"]])
    M[ic + 1L, ic + 1L] <- -mp(p[["k21"]])
  }
  if (ncmt == 3L) {
    d <- d - mp(p[["k13"]])
    M[ic + 2L, ic] <- mp(p[["k13"]])
    M[ic, ic + 2L] <- mp(p[["k31"]])
    M[ic + 2L, ic + 2L] <- -mp(p[["k31"]])
  }
  M[ic, ic] <- d
  M
}

# central amount at lag u after a unit bolus into the depot (or into
# central when depot = FALSE), by scaling and squaring
bolus_ss <- function(ncmt, depot, p, u) {
  M <- mk_M(ncmt, depot, p)
  ic <- if (depot) 2L else 1L
  if (u <= 0) return(if (depot) mp(0) else mp(1))
  expm_ss(M * mp(u))[ic, 1L]
}

# The disposition eigenvalues, by bisection on the interlacing
# brackets and Newton polishing, all at 300 bits.
disp_roots <- function(ncmt, p) {
  ke <- mp(p[["ke"]])
  if (ncmt == 1L) return(list(ke))
  k12 <- mp(p[["k12"]]); k21 <- mp(p[["k21"]])
  if (ncmt == 2L) {
    b <- ke + k12 + k21
    dsc <- sqrt(b * b - 4 * ke * k21)
    return(list((b + dsc) / 2, (b - dsc) / 2))
  }
  k13 <- mp(p[["k13"]]); k31 <- mp(p[["k31"]])
  a2 <- ke + k12 + k13 + k21 + k31
  a1 <- ke * k21 + ke * k31 + k21 * k31 + k12 * k31 + k13 * k21
  a0 <- ke * k21 * k31
  f <- function(x) x^3 - a2 * x^2 + a1 * x - a0
  fp <- function(x) 3 * x^2 - 2 * a2 * x + a1
  ends <- sort(c(mp(0), min(k21, k31), max(k21, k31), a2 * 2))
  out <- list()
  for (j in 1:3) {
    lo <- ends[[j]]; hi <- ends[[j + 1L]]
    flo <- f(lo)
    for (it in 1:400) {
      mid <- (lo + hi) / 2
      if (as.numeric(sign(f(mid)) * sign(flo)) >= 0) lo <- mid else
        hi <- mid
    }
    x <- (lo + hi) / 2
    for (it in 1:80) {
      d <- fp(x)
      if (as.numeric(abs(d)) == 0) break
      x <- x - f(x) / d
    }
    out[[j]] <- x
  }
  out
}

bolus_pf <- function(ncmt, depot, p, u) {
  lam <- disp_roots(ncmt, p)
  zs <- if (ncmt == 1L) list() else if (ncmt == 2L)
    list(mp(p[["k21"]])) else
      list(mp(p[["k21"]]), mp(p[["k31"]]))
  nn <- length(lam)
  cf <- vector("list", nn)
  for (i in seq_len(nn)) {
    num <- mp(1)
    for (z in zs) num <- num * (lam[[i]] - z)
    den <- mp(1)
    for (j in seq_len(nn)) if (j != i)
      den <- den * (lam[[i]] - lam[[j]])
    cf[[i]] <- num / den
  }
  uu <- mp(u)
  if (!depot) {
    s <- mp(0)
    for (i in seq_len(nn)) s <- s + cf[[i]] * exp(-lam[[i]] * uu)
    return(s)
  }
  ka <- mp(p[["ka"]])
  s <- mp(0)
  for (i in seq_len(nn)) {
    d <- lam[[i]] - ka
    term <- if (as.numeric(abs(d)) == 0) uu * exp(-ka * uu) else
      (exp(-lam[[i]] * uu) - exp(-ka * uu)) / (ka - lam[[i]])
    s <- s + cf[[i]] * ka * term
  }
  s
}

relerr <- function(a, b) as.numeric(abs(mp(a) - b) / abs(b))
