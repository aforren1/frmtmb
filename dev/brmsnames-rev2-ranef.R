## Reviewer recheck, MAJOR 4 breadth: ranef(ds) and coef(ds) on every
## random-effect structure frm_sample() takes, against the ML accessors of
## a fit set to the same draw (draws_fit_at()), plus a count of NA cells.
## Also: does each r_ column hold the coefficient expand_b() builds?
##   Rscript dev/brmsnames-rev2-ranef.R
## Data seed 71, sampler seed 4, chains 1, iter 80; draws 1, 20, 40.
.libPaths(c("C:/Users/adf44/source/r/brmsnames-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
q <- function(e) suppressWarnings(suppressMessages(e))
try1 <- function(e) {
  tryCatch(e, error = function(err) paste("ERROR:", conditionMessage(err)))
}
q(library(frmtmb)); q(library(frmtmb.sample))
sns <- asNamespace("frmtmb.sample")
set.seed(71)
ng <- 12; nt <- 5
times <- c(0, 0.4, 1, 1.7, 3)
dd <- data.frame(g = factor(rep(sprintf("g%02d", seq_len(ng)), each = nt)),
                 tim = num_factor(rep(times, ng)),
                 x1 = rnorm(ng * nt), x2 = rnorm(ng * nt))
dd$f <- factor(rep(c("a", "b", "a", "b", "a"), ng))
dd$h <- factor(sample(1:4, ng * nt, TRUE))
u <- rnorm(ng, 0, 0.8)
dd$y <- 1 + 0.5 * dd$x1 + u[as.integer(dd$g)] +
  0.3 * dd$x2 * u[as.integer(dd$g)] + rnorm(ng * nt, 0, 0.5)
Q <- Matrix::bandSparse(ng, k = c(-1, 0, 1),
                        diagonals = list(rep(-0.4, ng - 1), rep(1.2, ng),
                                         rep(-0.4, ng - 1)))
dimnames(Q) <- list(levels(dd$g), levels(dd$g))
A <- solve(as.matrix(Q))
dimnames(A) <- dimnames(Q)
V <- matrix(c(1.2, 0.5, 0.5, 0.8), 2, 2)
xc <- round(runif(ng), 2) * 10
yc <- round(runif(ng), 2) * 10
dd$pos <- num_factor(xc, yc)[as.integer(dd$g)]
dd$one <- factor(1)
ms <- list(
  us = list(y ~ x1 + (1 + x2 | g)),
  us_ncp = list(y ~ x1 + (1 + x2 | g), ncp = TRUE),
  diag = list(y ~ x1 + (1 + x2 || g)),
  rr = list(y ~ x1 + rr(x1 + x2 | g, d = 1)),
  ar1 = list(y ~ x1 + ar1(tim + 0 | g)),
  ou = list(y ~ x1 + ou(tim + 0 | g)),
  cs = list(y ~ x1 + cs(tim + 0 | g)),
  homcs = list(y ~ x1 + homcs(tim + 0 | g)),
  hetar1 = list(y ~ x1 + hetar1(tim + 0 | g)),
  homdiag = list(y ~ x1 + homdiag(tim + 0 | g)),
  toep = list(y ~ x1 + toep(tim + 0 | g)),
  exp = list(y ~ x1 + exp(pos + 0 | one)),
  gr_cov = list(y ~ x1 + (1 | gr(g, cov = A))),
  gr_prec = list(y ~ x1 + (1 | gr(g, prec = Q))),
  equalto = list(y ~ x1 + equalto(f + 0 | g, V)),
  two_blocks = list(y ~ x1 + (1 | g) + (0 + x2 | h)),
  smooth_re = list(y ~ s(x1) + (1 | g))
)
d2 <- list(A = A, Q = Q, V = V)
flat <- function(r) {
  # the ML ranef: a list of level x coefficient tables
  unlist(lapply(r, function(t) as.numeric(as.matrix(t))))
}
for (nm in names(ms)) {
  M <- ms[[nm]]
  fit <- try1(q(frm(bf(M[[1]]), family = gaussian(), data = dd, data2 = d2)))
  if (is.character(fit)) {
    cat(sprintf("%-11s frm: %s\n", nm, substr(fit, 1, 120)))
    next
  }
  set.seed(4)
  ds <- try1(q(frm_sample(fit, chains = 1, iter = 80, refresh = 0, seed = 4,
                          reparameterize = isTRUE(M$ncp))))
  if (is.character(ds)) {
    cat(sprintf("%-11s frm_sample: %s\n", nm, substr(ds, 1, 120)))
    next
  }
  rd <- try1(ranef(ds, summary = FALSE))
  cd <- try1(coef(ds, summary = FALSE))
  if (is.character(rd)) {
    cat(sprintf("%-11s ranef(ds): %s\n", nm, substr(rd, 1, 140)))
    next
  }
  na_r <- sum(vapply(rd, function(a) sum(is.na(a)), 1))
  na_c <- if (is.list(cd)) sum(vapply(cd, function(a) sum(is.na(a)), 1)) else
    NA
  idx <- sns$draws_par_index(fit)
  worst <- 0; worst_c <- 0; nmis <- 0
  for (i in c(1, 20, 40)) {
    sh <- sns$draws_fit_at(ds, i, idx)
    mr <- try1(flat(ranef(sh)))
    dr <- unlist(lapply(rd, function(a) as.numeric(a[i, , ])))
    if (is.character(mr) || length(mr) != length(dr)) {
      nmis <- nmis + 1
      next
    }
    worst <- max(worst, max(abs(sort(mr) - sort(dr)) /
                              pmax(abs(sort(mr)), 1e-12)))
    if (is.list(cd)) {
      mc <- try1(flat(coef(sh)))
      dc <- unlist(lapply(cd, function(a) as.numeric(a[i, , ])))
      if (!is.character(mc) && length(mc) == length(dc)) {
        worst_c <- max(worst_c, max(abs(sort(mc) - sort(dc)) /
                                      pmax(abs(sort(mc)), 1e-12)))
      } else nmis <- nmis + 1
    }
  }
  cat(sprintf(paste0("%-11s groups %-14s ranef NA %3d  coef NA %3s  ",
                     "worst rel diff ranef %.2g coef %.2g  shape ",
                     "mismatches %d\n"),
              nm, paste(names(rd), collapse = ","), na_r,
              as.character(na_c), worst, worst_c, nmis))
}
