# Re-check of prilabel fix round 1: the shared Student-t density where
# D_j varies by level (two responses, the second on a subset of the
# levels), on the tape against my per-level reference. Data seed 21 as
# dev/prilabel-rev2-tshare.R; perturbation seed 35.
.libPaths(c("C:/Users/adf44/source/r/wt-prilabel-lib",
            "C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-prilabel"
suppressPackageStartupMessages(library(frmtmb))
sink(file.path(wt, "dev/prilabel-rev-log/rev2-tuneven.txt"), split = TRUE)
set.seed(21)
n <- 360
g <- factor(sample(sprintf("g%02d", 1:24), n, TRUE))
d <- data.frame(g, x = rnorm(n), z = rnorm(n))
w <- sqrt(5 / rchisq(24, 5))
d$y <- 1 + 0.5 * d$x + (w * rnorm(24, 0, 0.8))[d$g] + rnorm(n, 0, 0.6)
d$yz <- d$y + rnorm(n)
d$s2 <- as.integer(d$g) > 6L
fit <- frm(bf(y ~ x + (1 | gr(g, dist = "student"))) +
             bf(yz | subset(s2) ~ (1 + x | gr(g, dist = "student"))) +
             set_rescor(FALSE), data = d)
bks <- fit$frame$re_blocks
cat("levels per block", vapply(bks, function(bk) length(bk$levels), 0L),
    " t_n", vapply(bks, function(bk) bk$t_n %||% NA_integer_, 0L), "\n")
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
ref <- function(pars) {
  per <- list()
  for (bk in bks) {
    B <- matrix(pars$b[bk$b_idx], bk$dim)
    S <- as.matrix(vc(bk, pars$theta))
    for (j in seq_along(bk$levels)) {
      lv <- as.character(bk$levels[j])
      per[[lv]] <- c(per[[lv]], list(list(b = B[, j], S = S)))
    }
  }
  Dj <- vapply(per, function(pl) sum(lengths(lapply(pl, `[[`, "b"))), 0)
  cat("  D_j table:", paste(names(table(Dj)), table(Dj), sep = "x"), "\n")
  t_ll <- sum(vapply(per, function(pl) {
    mvt1(unlist(lapply(pl, `[[`, "b")),
         as.matrix(Matrix::bdiag(lapply(pl, `[[`, "S"))), bks[[1]]$dist_nu)
  }, 0))
  nll <- 0
  for (lp in fit$frame$linpreds) {
    if (!identical(lp$dpar, "mu")) next
    eta <- as.vector(lp$X %*% pars[[lp$par]][lp$idx]) +
      as.vector(lp$Z %*% pars$b)
    r <- lp$resp
    y <- fit$frame$y[[r]]
    k <- match(r, names(fit$spec$responses))
    nll <- nll - sum(dnorm(y, eta, exp(pars$betad[k]), log = TRUE))
  }
  nll - t_ll
}
p0 <- fit$obj$env$last.par.best
set.seed(35)
for (r in 1:4) {
  p <- p0 + if (r == 1) 0 else rnorm(length(p0), 0, 0.4)
  tape <- as.numeric(fit$obj$env$f(p))
  rf <- ref(split_par(p))
  cat(sprintf("  point %d: tape %.10f ref %.10f rel %.3e\n", r, tape, rf,
              (tape - rf) / abs(rf)))
}
sink()
