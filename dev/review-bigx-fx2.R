# Confirmation of the M2 fix: products and powers of shifted columns at
# 1e4 to 1e6, the capped warning on raw polynomials, and every
# "false convergence (8)" stop judged against the centered fit.
source("C:/Users/adf44/source/r/frmtmb-wt-bigx/dev/review-bigx-common.R")
wtxt <- function(r) if (!length(r$w)) "none" else
  paste(unique(paste0("[", r$wc, "] ", r$w)), collapse = " || ")
gain_of <- function(r) {
  x <- regmatches(r$w, regexpr("gains [0-9.e+-]+", r$w))
  if (length(x)) as.numeric(sub("gains ", "", x[1])) else NA
}
rows <- NULL
rec <- function(lbl, r, fz) {
  conv <- if (r$err) NA else r$f$opt$convergence
  c8 <- any(grepl("false convergence", r$w))
  rows <<- rbind(rows, data.frame(case = lbl, conv = conv, c8 = c8,
    atmax = any(grepl("at the maximum", r$w)), gain = gain_of(r),
    dll = if (r$err) NA else ll_of(r$f) - ll_of(fz),
    other = sum(!grepl("false convergence", r$w))))
}
cat("== products and powers\n")
for (s in 1:5) for (m in c(1e4, 1e5, 1e6)) {
  set.seed(40 + s); n <- 300
  z1 <- rnorm(n); z2 <- rnorm(n)
  d <- data.frame(z1 = z1, z2 = z2, w1 = z1 + m, w2 = z2 + m)
  d$y <- 0.5 + 0.4 * z1 - 0.3 * z2 + 0.3 * z1 * z2 + rnorm(n, 0, 0.5)
  d$yp <- rpois(n, exp(0.2 + 0.3 * z1 + 0.05 * z1^2))
  fz <- quiet(frm(bf(y ~ z1 * z2), data = d))
  r <- rfit(bf(y ~ w1 * w2), data = d)
  g <- glm(y ~ w1 * w2, data = d)
  # the same model on columns centered after forming the product
  dc <- d; dc$p <- d$w1 * d$w2; dc$p <- dc$p - mean(dc$p)
  dc$a <- d$w1 - mean(d$w1); dc$b <- d$w2 - mean(d$w2)
  lc <- lm(y ~ a + b + p, data = dc)
  Xs <- model.matrix(~ w1 * w2, d)
  sc <- function(X) X / rep(sqrt(colSums(X^2) / nrow(X)), each = nrow(X))
  Xc <- Xs; Xc[, -1] <- Xc[, -1] - rep(colMeans(Xs[, -1]), each = n)
  cat(sprintf(paste0("s%d gaussian w1*w2 m=%g: kept %s dll %+.3e | glm rank %d dll %+.3e | ",
                     "lm on centered product dll %+.3e, product t %.1f | scaled rank raw %d centered %d\n   %s\n   msg: %s\n"),
              s, m, paste(names(beta_of(r$f)), collapse = ","), ll_of(r$f) - ll_of(fz),
              g$rank, as.numeric(logLik(g)) - ll_of(fz),
              as.numeric(logLik(lc)) - ll_of(fz), summary(lc)$coefficients["p", 3],
              qr(sc(Xs), tol = 1e-11)$rank, qr(sc(Xc), tol = 1e-11)$rank,
              wtxt(r), paste(unique(r$msg), collapse = " | ")))
  rec(sprintf("gauss w1*w2 m=%g s%d", m, s), r, fz)
  fz <- quiet(frm(bf(yp ~ z1 + I(z1^2)), family = poisson(), data = d))
  r <- rfit(bf(yp ~ w1 + I(w1^2)), family = poisson(), data = d)
  cat(sprintf("s%d poisson w1+I(w1^2) m=%g: kept %s dll %+.3e\n   %s\n   msg: %s\n", s, m,
              paste(names(beta_of(r$f)), collapse = ","), ll_of(r$f) - ll_of(fz),
              wtxt(r), paste(unique(r$msg), collapse = " | ")))
  rec(sprintf("pois w1+I(w1^2) m=%g s%d", m, s), r, fz)
}
cat("== raw polynomials on [1000, 1010]\n")
for (s in 501:503) {
  set.seed(s); n <- 200
  xx <- runif(n, 1000, 1010)
  d <- data.frame(xx = xx, yy = sin((xx - 1000) / 2) + rnorm(n, 0, 0.3))
  for (deg in c(5, 6, 8)) {
    fo <- as.formula(sprintf("yy ~ poly(xx, %d, raw = TRUE)", deg))
    r <- rfit(bf(fo), data = d)
    kept <- names(beta_of(r$f))
    X <- model.matrix(fo, d)
    # is the warning's claim true: full rank at 1e-11 scaled, centered
    # kappa of all kept-at-rank columns above 1e10, and of the kept ones below
    sc <- function(X) X / rep(sqrt(colSums(X^2) / nrow(X)), each = nrow(X))
    Xc <- X; Xc[, -1] <- Xc[, -1] - rep(colMeans(X[, -1]), each = n)
    qa <- qr(sc(Xc), tol = 1e-11)
    kk <- function(M) { q <- qr(sc(M), tol = 1e-11); 1 / rcond(qr.R(q)[seq_len(q$rank), seq_len(q$rank)], triangular = TRUE) }
    colkeep <- match(sub("poly\\(xx, [0-9], raw = TRUE\\)", "p", kept), sub("poly\\(xx, [0-9], raw = TRUE\\)", "p", colnames(X)))
    g <- lm(yy ~ ., data = data.frame(yy = d$yy, X[, colkeep[-1], drop = FALSE]))
    cat(sprintf("s%d deg %d: kept %d cols, raw scaled rank %d, centered scaled rank %d, centered kappa of kept %.2g; ll %.6f lm on kept %.6f; any NaN SE %s\n   %s\n   msg: %s\n",
                s, deg, length(kept), qr(sc(X), tol = 1e-11)$rank, qa$rank,
                kk(Xc[, colkeep, drop = FALSE]), ll_of(r$f), as.numeric(logLik(g)),
                anyNA(se_beta(r$f)), wtxt(r), paste(unique(r$msg), collapse = " | ")))
  }
}
cat("== y ~ w far out (code 8 battery)\n")
for (s in 1:10) for (m in c(1e7, 1e8, 1e9)) for (fam in c("gaussian", "poisson")) {
  set.seed(900 + s); n <- 150; z <- rnorm(n)
  y <- if (fam == "gaussian") 0.5 + 0.4 * z + rnorm(n, 0, 0.5) else rpois(n, exp(0.5 + 0.4 * z))
  d <- data.frame(y = y, z = z, w = z + m)
  F <- if (fam == "gaussian") gaussian() else poisson()
  fz <- quiet(frm(bf(y ~ z), family = F, data = d))
  r <- rfit(bf(y ~ w), family = F, data = d)
  rec(sprintf("%s y~w m=%g s%d", fam, m, s), r, fz)
}
cat("== summary of stops\n")
print(table(code8 = rows$c8, at_max_text = rows$atmax))
c8 <- rows[rows$c8, ]
cat("code-8 fits:", nrow(c8), " gain range", format(range(c8$gain, na.rm = TRUE), digits = 3),
    " dll vs centered range", format(range(c8$dll, na.rm = TRUE), digits = 3),
    " short by more than 1e-6:", sum(c8$dll < -1e-6, na.rm = TRUE),
    " with another warning:", sum(c8$other > 0), "\n")
print(c8[c8$dll < -1e-7 | c8$other > 0 | is.na(c8$gain), ], row.names = FALSE)
