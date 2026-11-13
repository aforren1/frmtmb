# Lane gpstop, short check S1: the margin of a span rule on the median
# nearest-neighbour gap, read before the rule changes. For every fit of the
# review's 122-row sweep (fx-sweep-*.txt, ROW lines, both arms) and of
# fitedge's 192 rows (gpf-*-lane.tsv, gpf-*-base.tsv), the fitted length
# scale over the median nearest-neighbour gap of the distinct x, and for
# every kept relinearized descent (REG ... verdict TRUE) the same ratio.
# Ratios use the data's x (one-dimensional positions); for gp(by = ) the
# gap of all x. No frmtmb code. Args: none.
o <- "C:/Users/adf44/source/r/frmtmb-wt-gpstop/dev/gpstop-out"
mednn <- function(x) {
  u <- sort(unique(x))
  g <- diff(u)
  stats::median(pmin(c(Inf, g), c(g, Inf)))
}
# the sweep's data (gpstop-fx-sweep.R mk()), x only
bodies <- c("exp_a", "exp_b", "logis_a", "mm_V", "mm_K", "pow_a", "pow_b")
sweep_x <- function(tag) {
  p <- regmatches(tag, regexec("^(.+)_(gaussian|poisson|bernoulli|student)_n([0-9]+)_s([0-9]+)", tag))[[1]]
  set.seed(81100 + 1000 * match(p[3], c("gaussian", "poisson", "bernoulli",
                                         "student")) +
             100 * as.integer(p[5]) + as.integer(p[4]) +
             10 * match(p[2], bodies))
  sort(round(runif(as.integer(p[4]), 0, 10), 3))
}
rows <- list()
for (f in file.path(o, c("fx-sweep-gauss.txt", "fx-sweep-pb.txt",
                         "fx-sweep-student.txt"))) {
  l <- readLines(f)
  for (r in grep("^ROW ", l, value = TRUE)) {
    v <- strsplit(r, " +")[[1]]
    if (grepl("_k[0-9]+$", v[2]) || "ERROR" %in% v) next
    th <- as.numeric(strsplit(v[match("theta", v) + 1L], ",")[[1]])
    x <- sweep_x(v[2])
    conv <- as.integer(v[match("conv", v) + 1L])
    nw <- as.integer(v[match("nwarn", v) + 1L])
    rows[[length(rows) + 1L]] <- data.frame(
      set = "sweep", tag = v[2], arm = v[3], conv = conv, nwarn = nw,
      ratio = exp(th[2]) / mednn(x), range = exp(th[2]) / diff(range(x)))
  }
  for (r in grep("^REG .*verdict TRUE", l, value = TRUE)) {
    v <- strsplit(r, " +")[[1]]
    x <- sweep_x(v[2])
    lr <- as.numeric(v[match("ls/range", v) + 1L])
    rows[[length(rows) + 1L]] <- data.frame(
      set = "sweep-kept", tag = v[2], arm = "relin", conv = 0L, nwarn = 0L,
      ratio = lr * diff(range(x)) / mednn(x), range = lr)
  }
}
# fitedge's 192: data as dev/gpstop-gpfresh.R builds them
source("C:/Users/adf44/source/r/frmtmb-wt-gpstop/dev/review-fitedge-designs.R")
fe_x <- function(seed, n, ci) {
  set.seed(18000 + 100 * seed + n + ci)
  sort(runif(n, 0, 10))
}
for (s in c("fresh", "lane96")) for (arm in c("lane", "base")) {
  t <- utils::read.delim(file.path(o, sprintf("gpf-%s-%s.tsv", s, arm)),
                         colClasses = c(theta = "character",
                                        wclass = "character"))
  names(t)[1] <- "design"
  for (i in seq_len(nrow(t))) {
    th <- as.numeric(strsplit(t$theta[i], ",")[[1]])
    if (s == "fresh") {
      x <- rv_mk(t$design[i], t$seed[i], t$n[i],
                 match(t$curve[i], c("sin", "wig", "exp")))$x
    } else {
      x <- fe_x(t$seed[i], t$n[i],
                match(t$curve[i], c("sin", "wig", "exp", "step")))
    }
    # gp blocks: (1 | g) comes first in H; gp(by = ) holds two blocks in I
    ls <- switch(t$design[i], H = th[3], I = th[c(2, 4)], th[2])
    for (l in ls) {
      rows[[length(rows) + 1L]] <- data.frame(
        set = paste0("fe-", s), tag = paste(t$design[i], t$curve[i], t$n[i],
                                            t$seed[i]),
        arm = arm, conv = t$conv[i], nwarn = t$nwarn[i],
        ratio = exp(l) / mednn(x), range = exp(l) / diff(range(x)))
    }
  }
}
d <- do.call(rbind, rows)
clean <- d$conv == 0 & d$nwarn == 0
for (s in unique(d$set)) {
  i <- d$set == s
  q <- stats::quantile(d$ratio[i & clean], c(0, 0.01, 0.5), names = FALSE)
  cat(sprintf(paste0("MARGIN %-11s fits %3d clean %3d | clean ls/mednn min %.3g",
                     " 1%% %.3g median %.3g | below 0.1: clean %d, all %d | ",
                     "below 0.5: clean %d\n"),
              s, sum(i), sum(i & clean), q[1], q[2], q[3],
              sum(i & clean & d$ratio < 0.1), sum(i & d$ratio < 0.1),
              sum(i & clean & d$ratio < 0.5)))
}
lo <- d[clean & d$ratio < 1, ]
lo <- lo[order(lo$ratio), ]
for (k in seq_len(min(15L, nrow(lo)))) {
  cat(sprintf("MARGINLOW %s %s %s ls/mednn %.4g ls/range %.4g\n", lo$set[k],
              lo$tag[k], lo$arm[k], lo$ratio[k], lo$range[k]))
}
nc <- d[!clean & d$ratio < 0.1, ]
for (k in seq_len(nrow(nc))) {
  cat(sprintf("MARGINDEG %s %s %s conv %d nwarn %d ls/mednn %.4g\n",
              nc$set[k], nc$tag[k], nc$arm[k], nc$conv[k], nc$nwarn[k],
              nc$ratio[k]))
}
# the 192 rows' lane fits that the second reading moved (bcost-fx-*.txt)
for (k in c("fe-fresh E wig 120 2", "fe-fresh L wig 120 2",
            "fe-lane96 D step 60 1", "fe-lane96 D step 60 3")) {
  i <- which(paste(d$set, d$tag) == k & d$arm == "lane")
  cat(sprintf("MARGINMOVED %s ls/mednn %.3g ls/range %.4g\n", k, d$ratio[i],
              d$range[i]))
}
# the white-noise points of the short check, on MAJOR 2's positions
cat(sprintf("MARGINWHITE plain ls/mednn %.3g | tie ls/mednn %.3g\n",
            0.0094 / 0.0665, exp(-5.25) / 0.0665))
