# Merge check 9 (ci082 x arspat2): confint_varcorr() on the cases ci082
# added to tests/testthat/test-covstruct-seam.R ("a fixed correlation has
# no width on every rounding": xcs/xch k 0..19, the estimated (0 + t | g),
# the ou() grids) and on arspat2's sar cases (the short check's coverage
# set: 5x5 row-standardized rook W, rho 0.9, seeds 1..100), on one
# library. Writes an rds of every row for dev/rel083-kt/arspat2/mc9-cmp.R.
#   Rscript dev/rel083-kt/arspat2/mc9-confint-dump.R <lib> <tag>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
has_sp <- requireNamespace("frmtmb.spatial", quietly = TRUE) &&
  startsWith(normalizePath(find.package("frmtmb.spatial"), "/"),
             normalizePath(a[1], "/"))
cat("lib:", dirname(getNamespaceInfo("frmtmb", "path")),
    as.character(packageVersion("frmtmb")), "| spatial from this lib:",
    has_sp, "\n")
src <- parse("tests/testthat/test-covstruct-seam.R")
keep <- vapply(src, function(e) {
  is.call(e) && identical(e[[1L]], as.name("<-")) &&
    as.character(e[[2L]])[1] %in% c("xcs_C", "register_xcs", "xcs_data")
}, NA)
for (e in src[keep]) eval(e)
quiet <- function(expr) suppressWarnings(suppressMessages(expr))
row <- function(case, cv) {
  data.frame(case = case, term = cv$term, type = cv$type,
             estimate = cv$estimate, lwr = cv$lwr, upr = cv$upr)
}
out <- list()
register_xcs()
register_xcs("xch", 0.999)
for (st in c("xcs", "xch")) {
  d0 <- xcs_data(rho = if (st == "xch") 0.999 else 0.5)
  f <- stats::as.formula(paste0("y ~ x + ", st, "(0 + t | g)"))
  for (k in 0:19) {
    d <- d0
    d$y <- d0$y * (1 + 2e-15 * k)
    fit <- quiet(frm(f, family = gaussian(), data = d))
    out[[length(out) + 1L]] <- row(paste(st, "k", k),
                                   quiet(confint_varcorr(fit)))
  }
}
d0 <- xcs_data()
fu <- quiet(frm(bf(y ~ x + (0 + t | g)), family = gaussian(), data = d0))
out[[length(out) + 1L]] <- row("estimated (0 + t | g)",
                               quiet(confint_varcorr(fu)))
nf <- if (exists("num_factor")) num_factor else function(x) {
  factor(x, levels = sort(unique(x)))
}
grids <- list(list(seed = 1L, times = c(0, 0.5, 1, 2, 8, 20, 40)))
for (s in 1:3) {
  grids[[length(grids) + 1L]] <-
    list(seed = s, times = c(0, 0.5, 1, 2, 8, 20, 40, 1000))
}
for (gr in grids) {
  set.seed(gr$seed)
  times <- gr$times
  b <- stats::rnorm(50L * length(times))
  dou <- data.frame(y = rep(b, each = 2L) +
                      stats::rnorm(2L * length(b), 0, 0.5),
                    g = factor(rep(seq_len(50L), each = 2L * length(times))),
                    tim = nf(rep(rep(times, each = 2L), 50L)))
  fo <- quiet(frm(y ~ 1 + ou(tim + 0 | g), family = gaussian(), data = dou))
  out[[length(out) + 1L]] <- row(paste("ou times", length(times), "seed",
                                       gr$seed), quiet(confint_varcorr(fo)))
}
if (has_sp) {
  suppressPackageStartupMessages(library(frmtmb.spatial))
  rook <- function(k) {
    g <- expand.grid(r = seq_len(k), c = seq_len(k))
    n <- nrow(g)
    outer(seq_len(n), seq_len(n), function(i, j) {
      as.numeric(abs(g$r[i] - g$r[j]) + abs(g$c[i] - g$c[j]) == 1)
    })
  }
  rho <- 0.9
  W <- rook(5)
  W <- W / rowSums(W)
  n <- nrow(W)
  for (seed in 1:100) {
    set.seed(seed)
    d <- data.frame(x = rnorm(n))
    d$y <- as.vector(solve(diag(n) - rho * W, 1 + 0.5 * d$x + rnorm(n)))
    fit <- quiet(frm(y ~ x + sar(W), data = d, data2 = list(W = W)))
    out[[length(out) + 1L]] <- row(paste("sar seed", seed),
                                   quiet(confint_varcorr(fit)))
  }
}
res <- do.call(rbind, out)
saveRDS(res, sprintf("dev/rel083-log/kt/arspat2/mc9-%s.rds", a[2]))
cat("rows", nrow(res), "cases", length(unique(res$case)), "\n")
