# Fix round 2, N2: the build-time cost of checking a margin slot in full
# (every column of G and v against whiten(), off the tape) against the
# two-column check it replaces. Seconds per call of
# autocor_seam_margin_setup() (full check) and of the two-column work
# alone, minimum of 5 rounds, on the cost designs: sar error n = 400,
# m = 300 and fcor n = 400, m = 200 (data seed 11, missing seed m).
lib <- "C:/Users/adf44/source/r/wt-arspat2-lib"
.libPaths(c(lib, "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
  library(frmtmb.spatial)
})
ns <- asNamespace("frmtmb")
g <- expand.grid(r = 1:20, c = 1:20)
n <- nrow(g)
Aj <- outer(seq_len(n), seq_len(n), function(i, j) {
  as.numeric(abs(g$r[i] - g$r[j]) + abs(g$c[i] - g$c[j]) == 1)
})
W <- Aj / rowSums(Aj)
A0 <- diag(n) - 0.4 * W
M <- solve(crossprod(A0))
set.seed(11)
d <- data.frame(x = rnorm(n), h = factor(rep(1:8, length.out = n)))
d$y <- as.vector(solve(A0, 1 + 0.5 * d$x + rnorm(8)[d$h] + rnorm(n)))
best <- function(f) min(replicate(5, system.time(f())[["elapsed"]]))
for (cs in list(list(ty = "sar", m = 300), list(ty = "fcor", m = 200))) {
  dm <- d
  set.seed(cs$m)
  dm$y[sort(sample(n, cs$m))] <- NA
  fo <- if (cs$ty == "fcor") y | mi() ~ x + (1 | h) + fcor(M) else
    y | mi() ~ x + (1 | h) + sar(W, type = "error")
  ob <- suppressWarnings(suppressMessages(
    frm(fo, data = dm, data2 = list(M = M, W = W), family = student(),
        dry_run = "objective")))
  ac <- ob$frame$autocor$y
  e <- ns$ac_block_entry(ac)
  rows <- which(is.na(dm$y))
  i <- seq_len(ac$d)
  sg <- if (isTRUE(e$sigma_varies)) exp(0.3 * sin(i)) else 1.37
  th <- e$start(ac$aux) + 0.3
  th <- th[seq_along(e$start(ac$aux))]
  u <- sin(1.3 * i)
  t_full <- best(function() ns$autocor_seam_margin_setup(e, ac, rows))
  t_two <- best(function() {
    e$margin(rows, ac$aux)(u, sg, th)
    ns$autocor_seam_margin_cols(e, ac$aux, sg, th,
                                rows[unique(c(1L, length(rows)))], u)
  })
  t_fit <- system.time(suppressWarnings(suppressMessages(
    frm(fo, data = dm, data2 = list(M = M, W = W), family = student(),
        prior = set_prior("gamma(2, 0.1)", class = "nu")))))[["elapsed"]]
  cat(sprintf("%s n %d m %d: setup with full check %.3f s, two-column work %.3f s, one fit %.2f s\n",
              cs$ty, n, length(rows), t_full, t_two, t_fit))
}
