# Lane arstep, fix round: fewer descents on the AR marginal? "full" starts
# from the best cells of 6 length scales, the best 3 cells and the best cell
# at each AR coefficient; "trim" from the best cells of 2 length scales, the
# best cell and the best cell at each AR coefficient. One process, arms
# interleaved by row. Rows: dev/arstep-phigrid.R's 32 ar() rows (seed
# 74000 + 31 seed + n) and the reviewer's p600 and p1000 panels
# (dev/review-arstep-cost.R). Args: none
source("C:/Users/adf44/source/r/frmtmb-wt-arstep/dev/arstep-lib.R")
rev_libs("lane")
ns <- asNamespace("frmtmb")
orig <- ns$gp_marginal_opt
txt <- paste(deparse(orig), collapse = "\n")
old <- "ord <- unique(c(ord, o1[!duplicated(pk[o1])]))"
new <- "ord <- unique(c(o0[seq_len(min(2L, length(o0)))], order(vals)[1L], o1[!duplicated(pk[o1])]))"
t2 <- sub(old, new, txt, fixed = TRUE)
stopifnot(t2 != txt)
trim <- eval(parse(text = t2)); environment(trim) <- ns
put <- function(f) { unlockBinding("gp_marginal_opt", ns); assign("gp_marginal_opt", f, envir = ns); lockBinding("gp_marginal_opt", ns) }
curves <- list(mix = function(x) 2 * sin(3 * x) + 0.5 * x, step = function(x) 2 * (x > 5) - 1)
des <- list()
for (cn in names(curves)) for (n in c(40, 80, 160)) for (seed in 1:5) des[[length(des) + 1]] <- c(cn, n, seed)
for (extra in list(c("step", 160, 2), c("step", 80, 22), c("step", 80, 23))) des[[length(des) + 1]] <- extra
mkd <- function(rw) {
  if (rw[1] == "p600") {
    set.seed(9600); S <- 12; T <- 50
    d <- expand.grid(t = seq_len(T), s = factor(seq_len(S))); d$x <- d$t / 5
    d$y <- sin(d$x) + rnorm(S, 0, 0.3)[d$s] +
      unlist(lapply(seq_len(S), function(i) as.numeric(arima.sim(list(ar = 0.5), T, sd = 0.2))))
    return(list(d = d, f = y ~ gp(x) + ar(time = t, gr = s)))
  }
  if (rw[1] == "p1000") {
    set.seed(9601); S <- 25; T <- 40
    d <- expand.grid(t = seq_len(T), s = factor(seq_len(S))); d$x <- d$t / 4
    d$y <- sin(d$x) + unlist(lapply(seq_len(S), function(i) as.numeric(arima.sim(list(ar = 0.5), T, sd = 0.2))))
    return(list(d = d, f = y ~ gp(x) + ar(time = t, gr = s)))
  }
  cn <- rw[1]; n <- as.integer(rw[2]); seed <- as.integer(rw[3])
  set.seed(74000 + seed * 31 + n)
  x <- sort(round(runif(n, 0, 10), 4))
  list(d = data.frame(x = x, t = seq_len(n), y = curves[[cn]](x) + as.numeric(arima.sim(list(ar = 0.4), n, sd = 0.03))),
       f = y ~ gp(x) + ar(time = t))
}
des <- c(list(c("p600"), c("p1000"), c("p600"), c("p1000")), des)
tot <- c(full = 0, trim = 0)
for (rw in des) {
  m <- mkd(rw)
  out <- sprintf("TR %-12s", paste(rw, collapse = " "))
  for (arm in if (runif(1) < 2) c("full", "trim")) {
    put(if (arm == "full") orig else trim)
    t0 <- proc.time()[[3]]
    r <- rev_fit(m$f, m$d, family = gaussian())
    el <- proc.time()[[3]] - t0
    tot[arm] <- tot[arm] + el
    out <- paste(out, sprintf("| %s ll %.6f w %d wall %.2f", arm, as.numeric(logLik(r$fit)), length(r$warn), el))
  }
  cat(out, "\n"); flush.console()
}
put(orig)
cat(sprintf("TRSUM wall full %.1f trim %.1f\n", tot["full"], tot["trim"]))
cat("DONE\n")
