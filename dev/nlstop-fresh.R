# Review of lane nlstart: fresh nonlinear designs, not the lane's 450 rows.
# Copy for lane nlstop (rows unchanged). Usage: NLSTOP_ARM=after|before Rscript dev/review-nlstart-fresh.R <fam,...> [ref]
# Prints ROW lines: the default fit (and, on the after arm, the default with
# the probe off), and with `ref` the best of 17 started fits (16 random,
# seed seed*13+7, plus the truth), started fits with se = FALSE.
args <- commandArgs(trailingOnly = TRUE)
fams <- strsplit(args[1], ",")[[1]]
do_ref <- length(args) > 1 && args[2] == "ref"
source("dev/nlstop-common.R")
F <- list(gaussian = gaussian(), poisson = poisson(), bernoulli = bernoulli(),
          binomial = binomial(), gamma = Gamma("log"),
          cumulative = cumulative())
pos <- function(v) pmax(v, 0.05)
# each body: formula pieces, truth, the scale of random starts per nlpar
B <- list(
  mm = list(f = "Vm * xp / (K + xp)", p = c(Vm = "1", K = "1"),
            lat = function(d) 3 * d$xp / (1 + d$xp),
            tr = c(Vm = 3, K = 1), sc = c(Vm = 2, K = 2)),
  pow = list(f = "a * xp^b", p = c(a = "1", b = "1"),
             lat = function(d) 1 * d$xp^0.5, tr = c(a = 1, b = 0.5),
             sc = c(a = 1, b = 1)),
  logis3 = list(f = "Asym / (1 + exp((xmid - x) / scal))",
                p = c(Asym = "1", xmid = "1", scal = "1"),
                lat = function(d) 3 / (1 + exp((0.3 - d$x) / 0.7)),
                tr = c(Asym = 3, xmid = 0.3, scal = 0.7),
                sc = c(Asym = 2, xmid = 1, scal = 1)),
  gomp = list(f = "Asym * exp(-b2 * b3^xq)",
              p = c(Asym = "1", b2 = "1", b3 = "1"),
              lat = function(d) 3 * exp(-2 * 0.5^d$xq),
              tr = c(Asym = 3, b2 = 2, b3 = 0.5),
              sc = c(Asym = 2, b2 = 2, b3 = 0.5)),
  sig4 = list(f = "lo + (hi - lo) / (1 + exp((xmid - x) / scal))",
              p = c(lo = "1", hi = "1", xmid = "1", scal = "1"),
              lat = function(d) 0.5 + 2.5 / (1 + exp((0 - d$x) / 0.6)),
              tr = c(lo = 0.5, hi = 3, xmid = 0, scal = 0.6),
              sc = c(lo = 1, hi = 2, xmid = 1, scal = 1)),
  sat_z = list(f = "b1 * (1 - exp(-b2 * x))", p = c(b1 = "1 + z", b2 = "1"),
               lat = function(d) 1.5 * exp(0.6 * d$x) + 0.3 * d$z,
               tr = NULL, sc = c(b1 = 1, b2 = 1)),
  two_sat = list(f = "b1 * (1 - exp(-b2 * x)) + b3 * (1 - exp(-b4 * z))",
                 p = c(b1 = "1", b2 = "1", b3 = "1", b4 = "1"),
                 lat = function(d) 1.5 * exp(0.6 * d$x) +
                   1.5 * (1 - exp(-0.6 * d$z)),
                 tr = NULL, sc = c(b1 = 1, b2 = 1, b3 = 1, b4 = 1)),
  two_sat_rev = list(f = "b1 * (1 - exp(-b2 * z)) + b3 * (1 - exp(-b4 * x))",
                     p = c(b1 = "1", b2 = "1", b3 = "1", b4 = "1"),
                     lat = function(d) 1.5 * exp(0.6 * d$x) +
                       1.5 * (1 - exp(-0.6 * d$z)),
                     tr = NULL, sc = c(b1 = 1, b2 = 1, b3 = 1, b4 = 1)),
  bigx = list(f = "b1 * exp(b2 * xb)", p = c(b1 = "1", b2 = "1"),
              lat = function(d) 0.5 * exp(6e-4 * d$xb),
              tr = c(b1 = 0.5, b2 = 6e-4), sc = c(b1 = 1, b2 = 1e-3)),
  five = list(f = "b1 * exp(b2 * x) + b3 * (1 - exp(-b4 * z)) + b5 * w",
              p = c(b1 = "1", b2 = "1", b3 = "1", b4 = "1", b5 = "1"),
              lat = function(d) 1 * exp(0.5 * d$x) -
                1 * (1 - exp(0.5 * d$z)) + 0.3 * d$w,
              tr = c(b1 = 1, b2 = 0.5, b3 = -1, b4 = -0.5, b5 = 0.3),
              sc = c(b1 = 1, b2 = 1, b3 = 1, b4 = 1, b5 = 1)),
  sat_re = list(f = "b1 * (1 - exp(-b2 * x))",
                p = c(b1 = "1 + (1 | g)", b2 = "1"),
                lat = function(d) (1.5 + d$ug) * exp(0.6 * d$x),
                tr = NULL, sc = c(b1 = 1, b2 = 1)),
  sat_sig = list(f = "b1 * (1 - exp(-b2 * x))", p = c(b1 = "1", b2 = "1"),
                 sig = TRUE,
                 lat = function(d) 1.5 * exp(0.6 * d$x) + 0.3 * d$z,
                 tr = NULL, sc = c(b1 = 1, b2 = 1)))
gen <- function(fam, body, seed) {
  set.seed(seed)
  n <- 300
  d <- data.frame(x = round(rnorm(n), 3), z = round(rnorm(n), 3),
                  w = round(rnorm(n), 3), xp = round(runif(n, 0.1, 5), 3),
                  xq = round(runif(n, 0, 5), 3),
                  xb = round(runif(n, 0, 3000), 1),
                  g = factor(sample(1:15, n, TRUE)))
  d$ug <- rnorm(15, 0, 0.4)[d$g]
  lat <- B[[body]]$lat(d)
  sdv <- if (isTRUE(B[[body]]$sig)) exp(-1 + 0.5 * d$z) else 0.5
  d$nt <- 10L
  d$y <- switch(fam,
    gaussian = lat + rnorm(n, 0, sdv),
    poisson = rpois(n, pos(lat)),
    bernoulli = rbinom(n, 1, plogis(lat - 2)),
    binomial = rbinom(n, 10, plogis(lat - 2)),
    gamma = rgamma(n, 3, 3 / pos(lat)),
    cumulative = ordered(cut(lat + rlogis(n), c(-Inf, 0.5, 1.5, 3, Inf),
                             labels = 1:4)))
  d
}
form <- function(fam, body) {
  b <- B[[body]]
  lhs <- if (fam == "binomial") "y | trials(nt)" else "y"
  args <- c(list(as.formula(paste(lhs, "~", b$f))),
            lapply(names(b$p), function(p) as.formula(paste(p, "~", b$p[[p]]))))
  if (isTRUE(b$sig)) args <- c(args, list(sigma ~ z))
  do.call(bf, c(args, list(nl = TRUE)))
}
cat("ROWHEAD fam body seed ll conv t nwarn wclass multi ll_off conv_off",
    "t_off ref nref b_nl\n")
only <- Sys.getenv("REV_BODIES", "")
bodies <- if (nzchar(only)) strsplit(only, ",")[[1]] else names(B)
seeds <- as.numeric(strsplit(Sys.getenv("REV_SEEDS", "11,22,33,44,55"),
                             ",")[[1]])
for (fam in fams) for (body in bodies) for (seed in seeds) {
  if (isTRUE(B[[body]]$sig) && fam != "gaussian") next
  d <- gen(fam, body, seed * 1000 + match(body, names(B)))
  fm <- form(fam, body)
  t0 <- proc.time()[["elapsed"]]
  r <- tt(frm(fm, family = F[[fam]], data = d))
  t1 <- proc.time()[["elapsed"]] - t0
  a <- llv(r)
  conv <- if (iserr(r)) NA else r$opt$convergence
  wcl <- attr(r, "wc")
  multi <- any(attr(r, "mc") == "frmtmb_nl_multistart") ||
    any(grepl("sign", attr(r, "m")))
  ll_off <- NA; conv_off <- NA; t_off <- NA
  if (arm != "before") {
    t0 <- proc.time()[["elapsed"]]
    ro <- withr::with_options(list(frmtmb.nl_sign_probe = FALSE),
                              tt(frm(fm, family = F[[fam]], data = d)))
    t_off <- proc.time()[["elapsed"]] - t0
    ll_off <- llv(ro)
    conv_off <- if (iserr(ro)) NA else ro$opt$convergence
  }
  best <- NA; nref <- 0
  bnl <- if (!iserr(r)) {
    bb <- r$estimates$beta
    paste(sprintf("%s=%.4g", names(bb), bb), collapse = ";")
  } else paste("ERR", gsub("\\s+", "_", substr(r$msg, 1, 80)))
  if (do_ref && !iserr(r)) {
    nm <- names(r$estimates$beta)
    np <- sub("_.*$", "", nm)
    sc <- B[[body]]$sc[np]
    sc[is.na(sc)] <- 0.5
    set.seed(seed * 13 + 7)
    starts <- lapply(1:16, function(i) {
      setNames(rnorm(length(nm)) * sc * if (i > 8) 3 else 1, nm)
    })
    if (!is.null(B[[body]]$tr)) {
      tr <- setNames(rep(0, length(nm)), nm)
      ic <- grepl("_\\(Intercept\\)$", nm)
      tr[ic] <- B[[body]]$tr[np[ic]]
      starts <- c(starts, list(tr))
    }
    best <- -Inf
    for (st in starts) {
      v <- llv(tt(frm(fm, family = F[[fam]], data = d, se = FALSE,
                      start = list(beta = st))))
      if (!is.na(v)) {
        nref <- nref + 1
        if (v > best) best <- v
      }
    }
  }
  cat(sprintf("ROW %s %s %d %.6f %s %.3f %d %s %d %.6f %s %.3f %.6f %d %s\n",
              fam, body, seed, a, conv, t1, length(wcl),
              if (length(wcl)) paste(unique(wcl), collapse = "+") else "-",
              as.integer(multi), ll_off, conv_off, t_off, best, nref, bnl))
  flush.console()
}
cat("DONE\n")
