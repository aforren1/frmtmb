# Short check, item 3: the margin slot. (a) a slot that forgets sigma in G
# is refused by the term's name; (b) a term without the slot (fallback, one
# whiten() per missing value) against the observed-rows fit; (c) is the
# build-time check a proof? A sar() slot with W_mm + W_mm' written 2 W_mm
# (right only for a symmetric W), on a row-standardized W, missing rows
# whose first and last are not neighbors. Lane only; registrations live in
# this process only.
source("dev/reviews/arspat2-short-common.R")
myfix_build <- function(spec, rows, n_all, lookup) {
  V <- lookup(spec$V, "myfix()")[rows, rows, drop = FALSE]
  L <- t(chol(V))
  list(L = L, Linv = forwardsolve(L, diag(nrow(L))),
       logdet = sum(log(diag(L))))
}
myfix_whiten <- function(y, mu, sigma, theta, aux) {
  u <- as.vector(aux$Linv %*% ((y - mu) / sigma))
  n <- length(y)
  lsig <- if (length(sigma) == 1L) n * log(sigma) else sum(log(sigma))
  list(u = u, logjac = -aux$logdet - lsig)
}
myfix_unwhiten <- function(u, mu, sigma, theta, aux) mu + sigma * (aux$L %*% u)
reg_myfix <- function(...) {
  frmtmb_register_autocor("myfix", args = "V", build = myfix_build,
                          whiten = myfix_whiten, unwhiten = myfix_unwhiten,
                          replace = TRUE, ...)
}
n <- 60
V <- 0.5^abs(outer(1:n, 1:n, "-"))
set.seed(4)
d <- data.frame(x = rnorm(n), h = factor(rep(1:6, each = 10)))
d$y <- 1 + 0.5 * d$x + rnorm(6)[d$h] +
  as.vector(t(chol(V)) %*% rnorm(n)) * sqrt(4 / rchisq(1, 4))
miss <- c(2, 8, 15, 22, 29, 36, 43, 50, 57)
d$y[miss] <- NA
ok <- !is.na(d$y)

# (a)
wrong <- function(rows, aux) {
  Lm <- aux$Linv[, rows, drop = FALSE]
  G <- crossprod(Lm)
  function(u, sigma, theta) list(G = G, v = as.vector(t(Lm) %*% u) / sigma)
}
reg_myfix(margin = wrong, sigma_varies = FALSE)
e <- tryCatch(frm(y | mi() ~ x + (1 | h) + myfix(V), data = d,
                  data2 = list(V = V), family = student()),
              error = function(e) conditionMessage(e))
cat("(a) sigma forgotten:", if (is.character(e)) e else "FITTED (not refused)", "\n")
# the same slot with sigma = 1 cannot be told apart, but the check point
# is sigma 1.37, so it must still refuse; and with no mi() the slot is
# never read
f_nomiss <- tryCatch(quiet(frm(y ~ x + myfix(V), data = d[ok, ],
                               data2 = list(V = V[ok, ok]))),
                     error = function(e) conditionMessage(e))
cat("(a) same wrong slot, no missing value:",
    if (is.character(f_nomiss)) f_nomiss else "fits (slot unused)", "\n")

# (b) fallback: no slot, sigma per row and one sigma, gaussian and student
reg_myfix(sigma_varies = TRUE)
for (fam in list(gaussian(), student())) {
  for (sg in c("one", "per row")) {
    sf <- if (sg == "one") sigma ~ 1 else sigma ~ x
    a <- quiet(frm(bf(y | mi() ~ x + (1 | h) + myfix(V), sf), family = fam,
                   data = d, data2 = list(V = V)))
    b <- quiet(frm(bf(y ~ x + (1 | h) + fcor(Vo), sf), family = fam,
                   data = d[ok, ], data2 = list(Vo = V[ok, ok])))
    stopifnot(identical(names(a$opt$par), names(b$opt$par)))
    r <- vapply(0:2, function(k) {
      p <- pts(a$opt$par, k)
      reld(a$obj$fn(p), b$obj$fn(p))
    }, 0)
    cat(sprintf("(b) fallback %s sigma %s: unread %d, slot NULL %s, rel vs observed-rows fcor %s, logLik %.8f vs %.8f\n",
                fam$family, sg, length(a$frame$mi_map$y$unread),
                is.null(ns$ac_block_entry(a$frame$autocor$y)[["margin"]]),
                paste(sprintf("%.1e", r), collapse = " "),
                as.numeric(logLik(a)), as.numeric(logLik(b))))
  }
}

# (c) a sar slot right only for a symmetric W
sp <- asNamespace("frmtmb.spatial")
k <- 7
Wb <- rook(k)
W <- Wb / rowSums(Wb)
nn <- nrow(W)
set.seed(21)
ds <- data.frame(x = rnorm(nn))
ds$y <- as.vector(solve(diag(nn) - 0.6 * W, 1 + 0.5 * ds$x + rnorm(nn)))
# missing rows: interior neighbors of different degree; first and last
# are not neighbors of each other
ms <- c(1, 2, 8, 9, 17, 24, 25, 49)
cat("(c) W[first, last]", W[ms[1], ms[length(ms)]], W[ms[length(ms)], ms[1]],
    "; interior asymmetry max |W_mm - W_mm'|",
    max(abs(W[ms, ms] - t(W[ms, ms]))), "\n")
ds$y[ms] <- NA
fit_right <- quiet(frm(y | mi() ~ x + sar(W, type = "error"), data = ds,
                       data2 = list(W = W)))
sar_margin_sym <- function(rows, aux) {
  W <- aux[["W"]]
  m <- length(rows)
  Wm <- W[, rows, drop = FALSE]
  Wmm <- as.matrix(W[rows, rows, drop = FALSE])
  I <- as.vector(diag(m))
  S1 <- as.vector(2 * Wmm)
  S2 <- as.vector(as.matrix(Matrix::crossprod(Wm)))
  WmT <- Matrix::t(Wm)
  function(u, sigma, theta) {
    rho <- sp$sar_rho(theta, aux)
    G <- RTMB::matrix((I - rho * S1 + rho^2 * S2) / sigma^2, m, m)
    v <- (u[rows] - rho * as.vector(WmT %*% u)) / sigma
    list(G = G, v = v)
  }
}
frmtmb_register_autocor(
  "sar", args = c("M", "type"), parse = sp$sar_parse, build = sp$sar_build,
  whiten = sp$sar_whiten, unwhiten = sp$sar_unwhiten,
  jacobian = sp$sar_jacobian, margin = sar_margin_sym, npar = 1L,
  start = sp$sar_start, natural = sp$sar_natural, sigma_varies = FALSE,
  bounds = sp$sar_bounds, usage = "sar(M, type = \"lag\")",
  compat_note = sp$sar_compat_note, replace = TRUE,
  package = "frmtmb.spatial")
fit_wrong <- tryCatch(quiet(frm(y | mi() ~ x + sar(W, type = "error"),
                                data = ds, data2 = list(W = W))),
                      error = function(e) conditionMessage(e))
if (is.character(fit_wrong)) {
  cat("(c) wrong-for-asymmetric slot: REFUSED:", fit_wrong, "\n")
} else {
  p <- fit_right$opt$par
  cat(sprintf("(c) wrong-for-asymmetric slot: ACCEPTED; objective at the right fit's optimum %.6f vs %.6f (diff %.3e); logLik %.6f vs right %.6f\n",
              fit_wrong$obj$fn(p), fit_right$obj$fn(p),
              fit_wrong$obj$fn(p) - fit_right$obj$fn(p),
              as.numeric(logLik(fit_wrong)), as.numeric(logLik(fit_right))))
}
