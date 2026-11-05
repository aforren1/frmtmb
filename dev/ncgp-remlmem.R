# What an exact gp() under REML costs on the plain Hessian tape, by hand from the
# frame, and with b ahead of beta in the template under the compressed one, where
# TMB's reverse sweep read right on the hand tape. Seed 20261103.
#   Rscript dev/ncgp-remlmem.R <n> <plain|reorder>
a <- commandArgs(TRUE)
n <- as.integer(a[1]); arm <- a[2]
.libPaths(c("C:/Users/adf44/source/r/wt-ncgp-lib", "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
set.seed(20261103)
x <- sort(runif(n, 0, 10))
d <- data.frame(x = x, y = sin(x) + 0.05 * rnorm(n))
fr <- frm(bf(y ~ gp(x)), data = d, family = gaussian(), dry_run = "frame")
tpl <- fr$par_template
nll <- frmtmb::build_objective(fr)
t0 <- proc.time()[[3]]
tc <- RTMB::TapeConfig(); RTMB::TapeConfig(matmul = "atomic")
if (arm == "reorder") {
  tpl <- tpl[c(setdiff(names(tpl), "beta"), "beta")]
  TMB::config(tmbad.sparse_hessian_compress = 1L, DLL = "RTMB")
}
print(names(tpl))
o <- RTMB::MakeADFun(nll, tpl, random = c("b", "beta"), silent = TRUE)
force(o$env$spHess)
RTMB::TapeConfig(tc)
TMB::config(tmbad.sparse_hessian_compress = 0L, DLL = "RTMB")
tb <- proc.time()[[3]] - t0
t0 <- proc.time()[[3]]
op <- nlminb(o$par, o$fn, o$gr, control = list(iter.max = 200, eval.max = 300))
topt <- proc.time()[[3]] - t0
p1 <- op$par + 0.02
h <- 1e-4
fd <- vapply(seq_along(p1), function(i) {
  e <- replace(numeric(length(p1)), i, h); (o$fn(p1 + e) - o$fn(p1 - e)) / (2 * h)
}, 0)
g <- o$gr(p1)
cat(sprintf("NCGPRM n=%d %s build %.1fs opt %.1fs %s obj %.5f; max |gr - fd| %.2e of %.3g\n", n, arm, tb, topt, op$message, op$objective,
            max(abs(g - fd)), max(abs(fd))))
