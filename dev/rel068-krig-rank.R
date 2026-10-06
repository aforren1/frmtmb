# gp_krig_factor() after the 0.68.0 dense fallback (gpby review r2):
# its time against the lane's pivoted-only factor in both regimes, the
# low rank it was built for and the high rank where it lost, with a
# control; and its law against gp_krig_cov() in both.
#
#   Rscript dev/rel068-krig-rank.R > dev/rel068-log/krig-rank.txt
#
# Instrument (dev/lane-rules.md): the arms run interleaved in one
# process, each timed over a block of repeats grown past 1.2 s, and the
# minimum of 5 rounds is reported. The control is the lane's function
# against a second copy of itself, which must read 1.00.
.libPaths(c("C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("frmtmb", format(packageVersion("frmtmb")), "from",
    dirname(find.package("frmtmb")), "\n")
# the lane's pivoted-only factor, read from lane gpby's worktree
lane_src <- "C:/Users/adf44/source/r/frmtmb-wt-gpby/R/predict.R"
ex <- parse(lane_src, keep.source = FALSE)
def <- Filter(function(e) {
  is.call(e) && identical(as.character(e[[2]]), "gp_krig_factor")
}, as.list(ex))[[1]]
old <- eval(def[[3]], ns)
environment(old) <- ns
old2 <- eval(def[[3]], ns)
environment(old2) <- ns
new <- ns$gp_krig_factor
set.seed(5)
d <- data.frame(x = round(runif(60, 0, 6), 1))
d$y <- sin(d$x) + rnorm(60, 0, 0.3)
f1 <- frm(bf(y ~ gp(x)), data = d)
bk <- Filter(function(b) b$covstruct == "gp", f1$frame$re_blocks)[[1]]
krig_at <- function(ell, n) {
  fs <- f1
  if (!is.na(ell)) fs$estimates$theta[bk$theta_idx[2]] <- log(ell)
  nd <- data.frame(x = seq(-1, 9, length.out = n) + 1e-3)
  ed <- ns$lp_eta_design(fs, fs$frame$linpreds[["y.mu"]], nd, FALSE, FALSE)
  Filter(Negate(is.null), lapply(ed$sm_parts, `[[`, "krig"))[[1]]
}
block <- function(f, kg) {
  reps <- 1L
  repeat {
    t <- system.time(for (i in seq_len(reps)) f(kg))[[3]]
    if (t > 1.2) return(t / reps)
    reps <- reps * 2L
  }
}
law_err <- function(fc, kg) {
  S <- ns$gp_krig_cov(kg)
  M <- tcrossprod(fc$L) + diag(fc$white, nrow(fc$L))
  M <- kg$sd2 * M[fc$idx, fc$idx] * outer(kg$w, kg$w)
  max(abs(M - S)) / kg$sd2
}
cat(sprintf("%-6s %-5s %-5s | %-9s %-9s %-9s | %-7s %-7s | %s\n", "ell",
            "rows", "rank", "lane s", "new s", "ctrl s", "new/ln",
            "ctrl/ln", "law err new / lane (of sd^2)"))
for (cfg in list(c(NA, 1000), c(0.05, 100), c(0.05, 400), c(0.05, 1000))) {
  kg <- krig_at(cfg[1], cfg[2])
  tl <- tn <- tc <- Inf
  for (r in 1:5) {
    tl <- min(tl, block(old, kg))
    tn <- min(tn, block(new, kg))
    tc <- min(tc, block(old2, kg))
  }
  fo <- old(kg)
  fn <- new(kg)
  cat(sprintf("%-6s %-5d %-5s | %-9.4f %-9.4f %-9.4f | %-7.2f %-7.2f | %.2g / %.2g\n",
              if (is.na(cfg[1])) "fit" else format(cfg[1]),
              length(kg$rows), paste0(ncol(fo$L), "/", ncol(fn$L)), tl, tn,
              tc, tn / tl, tc / tl, law_err(fn, kg), law_err(fo, kg)))
}
cat("KRIG DONE\n")
