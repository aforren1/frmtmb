# Lane wcheck: what an mvbf() frame holds for the weights check.
#   Rscript dev/wcheck-probe-frame.R <lib|base>
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r20"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (a[1] == "base") c(base, ul) else c(a[1], base, ul))
suppressMessages(library(frmtmb))
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")

set.seed(20261021)
n <- 80
d <- data.frame(obs = factor(seq_len(n)), x = rnorm(n), xs = rnorm(n),
                w = runif(n, 0.5, 2), g = factor(rep(1:20, each = 4)))
eta <- 0.5 + 0.4 * d$x + rnorm(n, 0, 0.5)
d$y <- eta + rnorm(n, 0, exp(0.3 * d$xs))
d$y2 <- d$y + rnorm(n)
d$sub <- seq_len(n) <= 60
G <- gaussian()
fr <- suppressWarnings(frm(mvbf(bf(y | weights(w) ~ x + (1 | obs)),
                                bf(y2 | weights(w) ~ x + (1 | p | obs))),
                           family = G, data = d, dry_run = "frame"))
cat("responses:", names(fr$spec$responses), " n_obs:", fr$n_obs, "\n")
cat("linpreds:", names(fr$linpreds), "\n")
for (k in names(fr$linpreds)) {
  lp <- fr$linpreds[[k]]
  cat(k, " X:", paste(dim(lp$X), collapse = "x"), " Z:",
      if (is.null(lp$Z)) "NULL" else paste(dim(lp$Z), collapse = "x"), "\n")
}
for (bk in fr$re_blocks) {
  cat("block", block_written_label <- frmtmb:::block_written_label(bk),
      " cs:", bk$covstruct, " dim:", bk$dim, " c_idx:",
      paste(range(bk$c_idx), collapse = "-"), " resp:", bk$resp %||% "?",
      " dpar:", bk$dpar %||% "?", "\n")
}
cat("names(block 1):", names(fr$re_blocks[[1]]), "\n")
cat("aterm_values names:", names(fr$aterm_values), "\n")
str(lapply(fr$aterm_values, function(v) lapply(v, head, 3)))
cat("y lengths:", vapply(fr$y, length, 1L), "\n")

fs <- suppressWarnings(frm(mvbf(bf(y | weights(w) + subset(sub) ~ x +
                                     (1 | obs)),
                                bf(y2 | weights(w) ~ x + (1 | obs))),
                           family = G, data = d, dry_run = "frame"))
cat("subset frame n_obs:", fs$n_obs, " y lengths:", vapply(fs$y, length, 1L),
    " w lengths:", vapply(fs$aterm_values, function(v) length(v$weights), 1L),
    "\n")
for (k in names(fs$linpreds)) {
  lp <- fs$linpreds[[k]]
  cat(k, " Z:", if (is.null(lp$Z)) "NULL" else paste(dim(lp$Z), collapse = "x"),
      "\n")
}
cat("DONE\n")
