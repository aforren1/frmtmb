# Round 2, item 1, part four. Two questions the clamp raises.
#
#   A  How reachable is the NaN, as a function of the one thing a user
#      can see: the SPREAD of the rate constants. The lane's Boundary
#      section describes one construction, k13 at 1e-300 with k31 on a
#      root. dev/rev-lincmt-f6c.R found the NaN at k13 = 1.8e-05, so
#      the described region is not the region.
#   B  Does the clamp fix the RANDOM cases too, and what does it cost
#      the value on ordinary parameters, where nothing is colliding?
#
# Script path: dev/rev-lincmt-f6d.R. Seed 4242.
source("C:/Users/adf44/source/r/frmtmb-wt-lincmt/dev/rev-lincmt-src.R")

tt <- c(0.5, 1, 2, 4, 8, 12, 24, 48)
ev <- data.frame(time = 0, state = "depot", value = 100, ii = 8,
                 addl = 3L)
gradof <- function(p) {
  f <- function(th) sum(frm_lincmt(
    parms = list(ke = exp(th[1]), k12 = exp(th[2]), k21 = exp(th[3]),
                 k13 = exp(th[4]), k31 = exp(th[5]), ka = exp(th[6]),
                 V = 10),
    times = tt, ncmt = 3, depot = TRUE, events = ev))
  x <- log(p)
  as.numeric(MakeTape(f, x)$jacobian(x))
}
valof <- function(p) frm_lincmt(
  parms = list(ke = p[1], k12 = p[2], k21 = p[3], k13 = p[4],
               k31 = p[5], ka = p[6], V = 10),
  times = tt, ncmt = 3, depot = TRUE, events = ev)

ns <- asNamespace("frmtmb.ode")
orig <- get("lincmt_disp", envir = ns)
mk_clamped <- function(d) {
  src <- deparse(orig)
  hit <- grep("arg/sqrt(1 + lincmt_relu(arg * arg - 1))", src,
              fixed = TRUE)
  stopifnot(length(hit) == 1L)
  src[[hit]] <- paste0(
    "    arg <- arg/sqrt(1 + lincmt_relu(arg * arg - 1)); ",
    "arg <- arg - lincmt_relu(arg - (1 - ", format(d, digits = 17),
    ")) + lincmt_relu(-(1 - ", format(d, digits = 17), ") - arg)")
  fn <- eval(parse(text = paste(src, collapse = "\n")))
  environment(fn) <- ns
  fn
}

cat("\n=== A. NaN rate against the spread of the rate constants ===\n")
cat("Rates drawn log-uniform over a window of the given width,",
    "\ncentred so the geometric mean is 0.2 per hour. 3000 draws each.",
    "\nA spread of four decades is 1e-2 to 1e2 relative to that mean,",
    "\nwhich a poorly started or weakly identified fit reaches.\n\n")
set.seed(4242)
cat(sprintf("%10s %10s %8s %12s\n", "decades", "draws", "NaN",
            "rate"))
for (dec in c(2, 3, 4, 5, 6, 8, 10, 12)) {
  N <- 3000L
  bad <- 0L
  for (i in seq_len(N)) {
    p <- 0.2 * 10^runif(6, -dec / 2, dec / 2)
    if (!all(is.finite(gradof(p)))) bad <- bad + 1L
  }
  cat(sprintf("%10d %10d %8d %11.2f%%\n", dec, N, bad,
              100 * bad / N))
}

cat("\n=== B. the clamp on the random cases, and its cost ===\n")
set.seed(4242)
# collect draws that make the shipped code return NaN
bad <- list()
while (length(bad) < 40L) {
  p <- 0.2 * 10^runif(6, -5, 5)
  if (!all(is.finite(gradof(p)))) bad[[length(bad) + 1L]] <- p
}
# and a matched set of ordinary draws, where nothing is colliding
set.seed(99)
ok <- lapply(1:400, function(i) 10^runif(6, -2, 1))

unlockBinding("lincmt_disp", ns)
cat(sprintf("%10s %18s %22s %14s\n", "clamp d", "NaN cases fixed",
            "value change on those", "on 400 ordinary"))
for (dd in c(1e-12, 1e-14, 1e-16)) {
  v_before_bad <- lapply(bad, valof)
  v_before_ok <- lapply(ok, valof)
  assign("lincmt_disp", mk_clamped(dd), envir = ns)
  fixed <- sum(vapply(bad, function(p) all(is.finite(gradof(p))),
                      TRUE))
  chg <- function(pl, vb) max(vapply(seq_along(pl), function(i) {
    va <- valof(pl[[i]])
    s <- max(abs(vb[[i]]))
    if (!(s > 0)) 0 else max(abs(va - vb[[i]])) / s
  }, 0))
  cb <- chg(bad, v_before_bad)
  co <- chg(ok, v_before_ok)
  assign("lincmt_disp", orig, envir = ns)
  cat(sprintf("%10.0e %11d of %3d %22.3e %14.3e\n", dd, fixed,
              length(bad), cb, co))
}
assign("lincmt_disp", orig, envir = ns)
stopifnot(identical(get("lincmt_disp", envir = ns), orig))

cat("\n=== C. the clamped gradient against a central difference ===\n")
cat("On the 40 cases the shipped code returns NaN for, with the",
    "\nclamp at 1e-14: the tape against a central difference of the",
    "\nsame function.\n\n")
fd <- function(f, x, h) vapply(seq_along(x), function(j) {
  xp <- x; xp[[j]] <- xp[[j]] + h
  xm <- x; xm[[j]] <- xm[[j]] - h
  (f(xp) - f(xm)) / (2 * h)
}, 0)
assign("lincmt_disp", mk_clamped(1e-14), envir = ns)
rel <- numeric(0)
for (p in bad) {
  f <- function(th) sum(frm_lincmt(
    parms = list(ke = exp(th[1]), k12 = exp(th[2]), k21 = exp(th[3]),
                 k13 = exp(th[4]), k31 = exp(th[5]), ka = exp(th[6]),
                 V = 10),
    times = tt, ncmt = 3, depot = TRUE, events = ev))
  x <- log(p)
  g <- as.numeric(MakeTape(f, x)$jacobian(x))
  g1 <- fd(f, x, 1e-5)
  g2 <- fd(f, x, 1e-4)
  s <- max(abs(g1), abs(g1 - g2))
  if (s > 0 && all(is.finite(g)) && all(is.finite(g1)))
    rel <- c(rel, max(abs(g - g1)) / s)
}
cat("  cases compared:", length(rel), "  worst relative difference",
    format(max(rel), digits = 4), "\n")
cat("  (the bar is the central difference's own step sensitivity,",
    "\n   so a value near 1 means the two agree as well as the",
    "\n   difference knows its own answer)\n")
assign("lincmt_disp", orig, envir = ns)
