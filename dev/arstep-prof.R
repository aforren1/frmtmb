# Lane arstep, fix round: where the AR marginal's time goes. Rprof on the
# reviewer's p1000 and a300 designs (dev/review-arstep-cost.R), lane build,
# and counts of cell and projection calls. Args: design
args <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-arstep/dev/arstep-lib.R")
rev_libs("lane")
ns <- asNamespace("frmtmb")
id <- args[1]
if (id == "p600") {
  set.seed(9600); S <- 12; T <- 50
  d <- expand.grid(t = seq_len(T), s = factor(seq_len(S)))
  d$x <- d$t / 5
  d$y <- sin(d$x) + rnorm(S, 0, 0.3)[d$s] +
    unlist(lapply(seq_len(S), function(i) as.numeric(arima.sim(list(ar = 0.5), T, sd = 0.2))))
  f <- y ~ gp(x) + ar(time = t, gr = s)
} else if (id == "p1000") {
  set.seed(9601); S <- 25; T <- 40
  d <- expand.grid(t = seq_len(T), s = factor(seq_len(S)))
  d$x <- d$t / 4
  d$y <- sin(d$x) + unlist(lapply(seq_len(S), function(i)
    as.numeric(arima.sim(list(ar = 0.5), T, sd = 0.2))))
  f <- y ~ gp(x) + ar(time = t, gr = s)
} else {
  set.seed(9300); n <- 300
  x <- sort(runif(n, 0, 10))
  d <- data.frame(x = x, t = seq_len(n),
                  y = sin(x) + as.numeric(arima.sim(list(ar = 0.6), n, sd = 0.1)))
  f <- y ~ gp(x) + ar(time = t)
}
cnt <- new.env(); cnt$cell <- 0L; cnt$proj <- 0L; cnt$desc <- 0L
trace("gp_ar_cell", tracer = quote(assign("cell", get("cell", cnt) + 1L, envir = cnt)), print = FALSE, where = ns)
trace("gp_ar_proj", tracer = quote(assign("proj", get("proj", cnt) + 1L, envir = cnt)), print = FALSE, where = ns)
trace("gp_ar_descend", tracer = quote(assign("desc", get("desc", cnt) + 1L, envir = cnt)), print = FALSE, where = ns)
t0 <- proc.time()[[3]]
fr0 <- frm(f, data = d, family = gaussian())$frame; cnt$cell <- 0L; cnt$proj <- 0L; cnt$desc <- 0L; t0 <- proc.time()[[3]]
m <- ns$gp_marginal_max(fr0)
cat(sprintf("PRF %s marginal alone %.2f s, cells %d, projections %d, descents %d\n", id,
            proc.time()[[3]] - t0, cnt$cell, cnt$proj, cnt$desc))
Rprof(file.path(REV_OUT, paste0("prof-", id, ".out")), interval = 0.01)
t0 <- proc.time()[[3]]
fit <- frm(f, data = d, family = gaussian())
Rprof(NULL)
cat(sprintf("PRF %s fit %.2f s\n", id, proc.time()[[3]] - t0))
sp <- summaryRprof(file.path(REV_OUT, paste0("prof-", id, ".out")))
print(head(sp$by.total[, c("total.time", "total.pct")], 25))
cat("DONE\n")
