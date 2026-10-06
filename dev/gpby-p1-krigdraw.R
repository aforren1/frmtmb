# Punch round 1, m3 and m4: frmtmb.sample's kriging draw on the final
# build. (1) Rows at one unseen position carry one value, exactly.
# (2) posterior_epred() at 300 unseen positions over 300 draws, against
# a control of 300 rows at observed positions (no kriging), with the
# number of factorizations and where the time goes.
# Usage: Rscript dev/gpby-p1-krigdraw.R <lane|base>
arm <- commandArgs(TRUE)[1]
.libPaths(c(if (arm == "lane") "C:/Users/adf44/source/r/wt-gpby-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("lib:", find.package("frmtmb"), find.package("frmtmb.sample"), "\n")
options(mc.cores = 1)
# second argument "dense": the draw before the pivoted factor, a dense
# Cholesky of gp_krig_cov() at the distinct positions, as a control
if (identical(commandArgs(TRUE)[2], "dense")) {
  dense <- function(krig) {
    key <- pos_rowkey(krig$X)
    u <- !duplicated(key)
    ku <- krig
    ku$X <- krig$X[u, , drop = FALSE]
    ku$Xw <- krig$Xw[u, , drop = FALSE]
    ku$Ks <- krig$Ks[u, , drop = FALSE]
    if (!is.null(krig$P)) ku$P <- krig$P[u, , drop = FALSE]
    ku$rs <- krig$rs[u]
    ku$w <- rep(1, sum(u))
    S <- gp_krig_cov(ku)
    R <- chol(S)
    krig$w * as.vector(crossprod(R, stats::rnorm(nrow(S))))[
      match(key, key[u])]
  }
  environment(dense) <- asNamespace("frmtmb")
  utils::assignInNamespace("gp_krig_draw", dense, "frmtmb")
  cat("dense Cholesky draw\n")
}
set.seed(5)
n <- 60
d <- data.frame(x = round(stats::runif(n, 0, 6), 1))
d$y <- 0.5 + sin(d$x) + stats::rnorm(n, 0, 0.3)
ds <- suppressWarnings(suppressMessages(
  frm_sample(bf(y ~ gp(x)), data = d, family = gaussian(), chains = 1,
             iter = 600, refresh = 0, seed = 4)))
nd <- data.frame(x = c(d$x[1], 6.7, 6.7, 7.4, 7.4 + 1e-9, d$x[2]))
set.seed(1)
e <- posterior_epred(ds, newdata = nd)
cat(sprintf(paste0("rows 2,3 same unseen x: max |diff| %.3e, identical ",
                   "%s (sd of row 2 %.3e) | rows 4,5 1e-9 apart: max ",
                   "|diff| %.3e\n"), max(abs(e[, 2] - e[, 3])),
            identical(e[, 2], e[, 3]), sd(e[, 2]), max(abs(e[, 4] - e[, 5]))))
ns <- asNamespace("frmtmb")
cnt <- new.env()
cnt$draw <- 0L
cnt$cov <- 0L
if (arm == "lane") {
  suppressMessages({
    trace("gp_krig_draw", where = ns, print = FALSE,
          tracer = quote(assign("draw", get("draw", cnt) + 1L, cnt)))
    trace("gp_krig_cov", where = ns, print = FALSE,
          tracer = quote(assign("cov", get("cov", cnt) + 1L, cnt)))
  })
}
big <- data.frame(x = seq(6.05, 9, length.out = 300))
ctl <- data.frame(x = rep(sort(unique(d$x)), length.out = 300))
# CPU time as well as wall time: other lanes share the machine
time3 <- function(nd) {
  vapply(1:3, function(i) {
    st <- system.time(posterior_epred(ds, newdata = nd, ndraws = 300))
    c(wall = st[["elapsed"]], cpu = st[["user.self"]] + st[["sys.self"]])
  }, c(wall = 0, cpu = 0))
}
show <- function(lab, t) {
  cat(sprintf("TIME %s %s: wall %s s, cpu %s s (median cpu %.2f)\n", arm,
              lab, paste(sprintf("%.2f", t["wall", ]), collapse = " "),
              paste(sprintf("%.2f", t["cpu", ]), collapse = " "),
              median(t["cpu", ])))
}
tb <- time3(big)
nb <- c(cnt$draw, cnt$cov)
cnt$draw <- 0L
cnt$cov <- 0L
tc <- time3(ctl)
show("300 unseen rows x 300 draws", tb)
show("300 observed rows x 300 draws (control)", tc)
if (arm == "lane") {
  cat(sprintf("calls over 3 unseen runs: gp_krig_draw %d, gp_krig_cov %d\n",
              nb[1], nb[2]))
  cat(sprintf("calls over 3 control runs: gp_krig_draw %d, gp_krig_cov %d\n",
              cnt$draw, cnt$cov))
  suppressMessages({
    untrace("gp_krig_draw", where = ns)
    untrace("gp_krig_cov", where = ns)
  })
}
if (arm == "lane" && !identical(commandArgs(TRUE)[2], "dense")) {
  rk <- new.env()
  rk$r <- integer(0)
  suppressMessages(trace("gp_krig_factor", where = ns, print = FALSE,
                         exit = quote(assign("r", c(get("r", rk),
                                                    ncol(returnValue()$L)),
                                             rk))))
  invisible(posterior_epred(ds, newdata = big, ndraws = 300))
  suppressMessages(untrace("gp_krig_factor", where = ns))
  cat(sprintf("pivoted factor rank over 300 draws at 300 positions: %s\n",
              paste(names(summary(rk$r)), summary(rk$r), collapse = ", ")))
}
# again with nothing traced: the counts above wrap the functions
show("300 unseen rows x 300 draws, untraced", time3(big))
show("300 observed rows x 300 draws (control), untraced", time3(ctl))
pf <- tempfile()
Rprof(pf, interval = 0.005)
invisible(posterior_epred(ds, newdata = big, ndraws = 300))
Rprof(NULL)
sp <- summaryRprof(pf)$by.self
print(utils::head(sp, 12))
cat("DONE\n")
