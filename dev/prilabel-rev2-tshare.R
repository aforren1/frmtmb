# Re-check of prilabel fix round 1, finding 1: Student-t gr(by = )
# terms and the per-level shared scale. My own reference density (a
# closed-form multivariate t per level, stacked by level NAME across
# blocks), on the tape at several points; uneven level sets; an empty
# factor level; simulate(), predict() and quadrature. Data seed 21 (as
# dev/prilabel-rev-tby.R); perturbation seed 33; draw seed 34.
.libPaths(c("C:/Users/adf44/source/r/wt-prilabel-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-prilabel"
suppressPackageStartupMessages(library(frmtmb))
sink(file.path(wt, "dev/prilabel-rev-log/rev2-tshare.txt"), split = TRUE)
cat("frmtmb from", dirname(find.package("frmtmb")), "\n")
set.seed(21)
n <- 360
g <- factor(sample(sprintf("g%02d", 1:24), n, TRUE))
d <- data.frame(g, x = rnorm(n), z = rnorm(n),
                hg = factor(ifelse(as.integer(g) <= 12L, "u", "v")))
w <- sqrt(5 / rchisq(24, 5))
d$y <- 1 + 0.5 * d$x + (w * rnorm(24, 0, 0.8))[d$g] +
  (w * rnorm(24, 0, 0.4))[d$g] * d$z + rnorm(n, 0, 0.6)
d$yz <- d$y + rnorm(n)
try_ <- function(expr) tryCatch(expr, error = function(e) {
  paste("ERROR:", conditionMessage(e))
})
ll <- function(f) format(as.numeric(logLik(f)), digits = 10)
sharing <- function(f) deparse(lapply(f$frame$re_blocks, `[[`, "t_share"))
st <- function(...) NULL

cat("\n######## by-term alone (rellib-r9: -361.6239067, -364.5278394)\n")
for (f in list(y ~ x + (1 | gr(g, by = hg, dist = "student")),
               y ~ x + (1 + x | gr(g, by = hg, dist = "student")))) {
  d0 <- d
  d0$y <- 1 + 0.5 * d$x + rnorm(24, 0, 0.8)[d$g]
  # the rellib-r9 numbers came from dev/prilabel-rev-tby.R's data
  set.seed(21)
  g0 <- factor(sample(sprintf("g%02d", 1:24), n, TRUE))
  d0 <- data.frame(g = g0, x = rnorm(n), z = rnorm(n),
                   hg = factor(ifelse(as.integer(g0) <= 12L, "u", "v")))
  d0$y <- 1 + 0.5 * d0$x + rnorm(24, 0, 0.8)[d0$g] + rnorm(n, 0, 0.6)
  r <- try_(frm(bf(f), data = d0))
  q <- try_(frm(bf(f), data = d0, quadrature = TRUE))
  cat(deparse1(f), "\n  logLik", if (is.character(r)) r else ll(r),
      " t_share", if (!is.character(r)) sharing(r), "\n  quadrature:",
      if (is.character(q)) substr(q, 1, 100) else paste("fits", ll(q)), "\n")
}

cat("\n######## my reference: per level, stack every block that holds it\n")
vc <- function(bk, th) {
  frmtmb:::covstruct_registry[[bk$covstruct]]$vcov(th[bk$theta_idx], bk)
}
mvt1 <- function(b, S, nu) {
  D <- length(b)
  lgamma((nu + D) / 2) - lgamma(nu / 2) - D / 2 * log(nu * pi) -
    0.5 * as.numeric(determinant(S)$modulus) -
    (nu + D) / 2 * log1p(drop(crossprod(b, solve(S, b))) / nu)
}
split_par <- function(p) {
  out <- list()
  for (nm in unique(names(p))) out[[nm]] <- unname(p[names(p) == nm])
  out
}
ref_t <- function(fit, pars, data) {
  bks <- Filter(function(bk) !is.null(bk$dist_nu), fit$frame$re_blocks)
  nu <- bks[[1]]$dist_nu
  per <- list()
  for (bk in bks) {
    B <- matrix(pars$b[bk$b_idx], bk$dim)
    S <- as.matrix(vc(bk, pars$theta))
    for (j in seq_along(bk$levels)) {
      lv <- as.character(bk$levels[j])
      per[[lv]] <- c(per[[lv]], list(list(b = B[, j], S = S)))
    }
  }
  sum(vapply(per, function(pl) {
    b <- unlist(lapply(pl, `[[`, "b"))
    S <- as.matrix(Matrix::bdiag(lapply(pl, `[[`, "S")))
    mvt1(b, S, nu)
  }, 0))
}
gauss_nll <- function(fit, pars, data, yname = "y") {
  lp <- Filter(function(l) identical(l$dpar, "mu"), fit$frame$linpreds)[[1]]
  eta <- as.vector(lp$X %*% pars[[lp$par]][lp$idx]) +
    as.vector(lp$Z %*% pars$b)
  -sum(dnorm(data[[yname]], eta, exp(pars$betad[1]), log = TRUE))
}
tape_check <- function(fit, data, label) {
  cat("\n--", label, " t_share", sharing(fit), "\n")
  p0 <- fit$obj$env$last.par.best
  set.seed(33)
  for (r in 1:4) {
    p <- p0 + if (r == 1) 0 else rnorm(length(p0), 0, 0.4)
    pars <- split_par(p)
    tape <- as.numeric(fit$obj$env$f(p))
    ref <- gauss_nll(fit, pars, data) - ref_t(fit, pars, data)
    cat(sprintf("  point %d: tape %.10f ref %.10f rel %.3e\n", r, tape, ref,
                (tape - ref) / abs(ref)))
  }
}
F1 <- frm(bf(y ~ x + (1 | gr(g, by = hg, dist = "student")) +
               (0 + z | gr(g, dist = "student"))), data = d)
tape_check(F1, d, "scalar by-term + plain term on g")
F2 <- frm(bf(y ~ x + (1 + x | gr(g, by = hg, dist = "student")) +
               (0 + z | gr(g, dist = "student"))), data = d)
tape_check(F2, d, "correlated by-term + plain term on g")

cat("\n######## an empty factor level, and uneven level sets\n")
de <- d
de$g <- factor(de$g, levels = c(levels(d$g), "gZZ"))
Fe <- try_(frm(bf(y ~ x + (1 | gr(g, by = hg, dist = "student")) +
                    (0 + z | gr(g, dist = "student"))), data = de))
cat("unused level gZZ: logLik", if (is.character(Fe)) Fe else ll(Fe),
    " (without it", ll(F1), ")\n")
if (!is.character(Fe)) {
  cat("  levels per block:",
      vapply(Fe$frame$re_blocks, function(bk) length(bk$levels), 0L), "\n")
}
# a level of g held by one block only: the plain term sees only rows
# with z != 0, i.e. levels g01..g06 lose every row of that term
du <- d
du$z[as.integer(du$g) <= 6L] <- 0
Fu <- frm(bf(y ~ x + (1 | gr(g, by = hg, dist = "student")) +
               (0 + z | gr(g, dist = "student"))), data = du)
cat("z = 0 on g01..g06: blocks' levels",
    vapply(Fu$frame$re_blocks, function(bk) length(bk$levels), 0L), "\n")
tape_check(Fu, du, "plain term without data on g01..g06")
# two responses, the second observed on a subset of the levels
dm <- d
dm$s2 <- as.integer(dm$g) > 6L
Fm <- try_(frm(bf(y ~ x + (1 | gr(g, dist = "student"))) +
                 bf(yz | subset(s2) ~ (1 | gr(g, dist = "student"))) +
                 set_rescor(FALSE), data = dm))
if (is.character(Fm)) {
  cat("mv subset:", Fm, "\n")
} else {
  cat("mv subset: logLik", ll(Fm), " levels per block",
      vapply(Fm$frame$re_blocks, function(bk) length(bk$levels), 0L),
      " t_share", sharing(Fm), " t_n",
      vapply(Fm$frame$re_blocks, function(bk) bk$t_n %||% NA_integer_, 0L),
      "\n")
}

cat("\n######## draws: one mixing variable per level of g\n")
bks <- F1$frame$re_blocks
set.seed(34)
D <- replicate(4000, frmtmb:::draw_b(F1))
for (lv in c("g02", "g20")) {
  i_by <- unlist(lapply(bks[1:2], function(bk) {
    j <- match(lv, bk$levels)
    if (is.na(j)) NULL else bk$b_idx[j]
  }))
  i_z <- bks[[3]]$b_idx[match(lv, bks[[3]]$levels)]
  cat(sprintf("  level %s: cor(|b_by|, |b_z|) %.3f\n", lv,
              cor(abs(D[i_by, ]), abs(D[i_z, ]))))
}
i_a <- bks[[1]]$b_idx[1]
i_b <- bks[[3]]$b_idx[match(bks[[1]]$levels[2], bks[[3]]$levels)]
cat(sprintf("  different levels: cor %.3f\n", cor(abs(D[i_a, ]),
                                                  abs(D[i_b, ]))))

cat("\n######## simulate / predict on the shared model\n")
s1 <- try_(simulate(F1, nsim = 2, seed = 1))
cat("simulate():", if (is.character(s1)) s1 else "runs", "\n")
s2 <- try_(simulate(F1, nsim = 2, seed = 1,
                    re_formula = ~ (1 | gr(g, by = hg, dist = "student"))))
cat("simulate(keep by, drop plain):", substr(s2, 1, 100), "\n")
s3 <- try_(simulate(F1, nsim = 2, seed = 1, re_formula = NA))
cat("simulate(re_formula = NA):", if (is.character(s3)) s3 else "runs", "\n")
nd <- d[1:3, ]
nd$g <- factor(c("gN1", "gN2", "gN3"))
p1 <- try_(predict(F1, newdata = nd, allow_new_levels = TRUE)[, 1])
cat("predict new levels:", if (is.character(p1)) p1 else round(p1, 3), "\n")
s4 <- try_(simulate(F1, nsim = 2, seed = 1, newdata = nd,
                    allow_new_levels = TRUE))
cat("simulate new levels:", if (is.character(s4)) substr(s4, 1, 100) else
  "runs", "\n")

cat("\n######## quadrature = TRUE\n")
for (f in list(
  y ~ x + (1 | gr(g, by = hg, dist = "student")) + (0 + z | gr(g, dist = "student")),
  y ~ x + (1 | gr(g, by = hg, dist = "student")) + (1 | h0),
  y ~ x + (1 | gr(g, dist = "student")) + (1 | gr(h0, dist = "student")),
  y ~ x + (1 | gr(g, by = hg, dist = "student")) + (0 + z | g))) {
  d$h0 <- factor(sample(1:5, n, TRUE))
  q <- try_(frm(bf(f), data = d, quadrature = TRUE))
  cat(format(deparse1(f), width = 92), ":",
      if (is.character(q)) substr(q, 8, 90) else "fits", "\n")
}
sink()
