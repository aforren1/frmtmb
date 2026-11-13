# Fix round, decision 1b: where to put draw_edge_c. For every parameter
# side of a set of fits, the unscreened scale r = d / (k * csd), d the
# distance at which the log-likelihood falls by k^2 / 2 (r = 1: the normal
# law's own; capped at 1), from draw_edge_side() with draw_edge_c = 1.
# (a steep side's scale at z = 2). Then whether the shipped law (draw_edge_c as built) changes the fit's
# draws at all. Prints one line per fit: its smallest r, where, and the
# verdict. Usage: Rscript dev/ppdraw-fix-c.R <lib> <out.tsv>
args <- commandArgs(TRUE)
.libPaths(c(args[1], "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib:", find.package("frmtmb"), " draw_edge_c",
    frmtmb:::draw_edge_c, "\n")
ns <- asNamespace("frmtmb")
e1 <- new.env(parent = ns)
e1$draw_edge_c <- 1
side1 <- frmtmb:::draw_edge_side
environment(side1) <- e1
q <- function(expr) suppressMessages(suppressWarnings(expr))
rows <- list()
probe <- function(label, group, fit) {
  if (is.null(fit)) return(invisible())
  ds <- q(frmtmb:::fit_draw_space(fit))
  V <- ds$V
  if (is.null(V) || !all(is.finite(V))) {
    cat(sprintf("%-34s %-10s no usable covariance\n", label, group))
    return(invisible())
  }
  p <- fit$opt$par
  off <- length(ds$map$names) - length(p)
  sdv <- sqrt(pmax(diag(V), 0))
  csd <- frmtmb:::draw_edge_csd(fit, V, sdv, off, !is.null(ds$root))
  if (is.null(csd)) {
    cat(sprintf("%-34s %-10s no trusted csd (factor, no Hessian)\n", label,
                group))
    return(invisible())
  }
  obj <- fit$obj
  f0 <- obj$fn(p)
  r <- matrix(1, length(p), 2)
  for (k in seq_along(p)) {
    j <- off + k
    if (!is.finite(sdv[j]) || sdv[j] <= 0) next
    for (s in 1:2) {
      v <- side1(obj, p, k, c(-1, 1)[s], sdv[j], min(csd[j], sdv[j]),
                 f0)
      # a steep side's scale at z = 2: where the fall is 2, over 2 csd
      r[k, s] <- if (length(v) > 1L) v[2] / 2 else v
    }
  }
  obj$fn(p)
  nm <- rownames(V)[off + seq_along(p)]
  law <- frmtmb:::draw_edge_build(fit, ds$map, V, diag(V), !is.null(ds$root))
  rr <- r
  rr[rr == 0] <- NA
  m <- suppressWarnings(min(rr, na.rm = TRUE))
  w <- which(rr == m, arr.ind = TRUE)[1, ]
  changes <- !isFALSE(law)
  rows[[length(rows) + 1L]] <<- data.frame(
    label = label, group = group, min_r = m,
    where = paste0(nm[w[1]], if (w[2] == 1) " -" else " +"),
    sd_over_csd = sdv[off + w[1]] / csd[off + w[1]],
    flat_sides = sum(r == 0), changes = changes,
    narrowed = if (changes) paste(nm[law$narrowed - off], collapse = " ")
      else "",
    edge = if (changes) paste(nm[law$edge - off], collapse = " ") else "")
  cat(sprintf("%-34s %-10s min_r %6.3f at %-30s flat %d changes %s\n",
              label, group, m, rows[[length(rows)]]$where, sum(r == 0),
              changes))
}
n <- 300L
source("dev/ppdraw-sweep-gen.R")
for (fam in fams) for (case in c("interior", "boundary")) {
  s <- spec(fam)
  probe(paste(fam, case), "sweep",
        q(tryCatch(frm(s[[1]], family = s[[2]], data = gen(fam, case)),
                   error = function(e) NULL)))
}
for (seed in 1:30) {
  set.seed(seed)
  d <- data.frame(x = stats::rnorm(200))
  d$y <- stats::rpois(200, exp(1 + 0.3 * d$x))
  probe(paste("zip null seed", seed), "zipnull",
        q(frm(bf(y ~ x), family = zero_inflated_poisson(), data = d)))
}
for (seed in 1:40) {
  set.seed(seed)
  d <- data.frame(g = factor(rep(1:6, each = 20)), x = stats::rnorm(120))
  d$y <- 1 + 0.5 * d$x + stats::rnorm(6, 0, 0.1)[d$g] + stats::rnorm(120)
  probe(paste("gauss (1|g) 6 grp seed", seed), "g6",
        q(frm(bf(y ~ x + (1 | g)), data = d)))
}
for (sd_s in c(0.05, 0.1, 0.2)) for (seed in 1:5) {
  set.seed(seed)
  dd <- data.frame(g = factor(rep(1:20, each = 10)), x = stats::rnorm(200))
  dd$y <- 1 + 0.5 * dd$x + stats::rnorm(20, 0, 0.7)[dd$g] +
    stats::rnorm(20, 0, sd_s)[dd$g] * dd$x + stats::rnorm(200)
  probe(sprintf("slope sd %.2f seed %d", sd_s, seed), "slope",
        q(frm(bf(y ~ x + (1 + x | g)), data = dd)))
}
# correlated interior coefficients: an uncentered covariate
for (seed in 1:5) {
  set.seed(seed)
  d <- data.frame(x = stats::rnorm(200, 100, 1))
  d$y <- stats::rpois(200, exp(-29 + 0.3 * d$x))
  probe(paste("poisson uncentered x seed", seed), "corr",
        q(frm(bf(y ~ x), family = poisson(), data = d)))
  d$yb <- stats::rbinom(200, 1, stats::plogis(-30 + 0.3 * d$x))
  probe(paste("logit uncentered x seed", seed), "corr",
        q(frm(bf(yb ~ x), family = bernoulli(), data = d)))
}
data(sleepstudy, package = "lme4", envir = environment())
probe("sleepstudy (Days | Subject)", "corr",
      q(frm(Reaction ~ Days + (Days | Subject), data = sleepstudy)))
out <- do.call(rbind, rows)
utils::write.table(out, args[2], sep = "\t", row.names = FALSE, quote = FALSE)
cat("done\n")
