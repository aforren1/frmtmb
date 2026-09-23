# A 9.5e-6 logLik gap seen in the lane's first exploratory run, where the
# two responses shared one group effect exactly. sim_v0() is that run's
# generator: y2 drawn right after yo, with 0.7 times y's group effect, so
# the generating group correlation is 1.
source("C:/Users/adf44/source/r/frmtmb-wt-drmtmb/dev/drmtmb-agree-lib.R")
sim_v0 <- function(seed, ng = 40, nper = 10) {
  set.seed(seed)
  n <- ng * nper
  g <- factor(rep(seq_len(ng), each = nper))
  u <- rnorm(ng, 0, 0.6)[g]
  v <- rnorm(ng, 0, 0.3)[g]
  x <- rnorm(n)
  z <- rnorm(n)
  d <- data.frame(g = g, x = x, z = z, w = rnorm(ng)[g])
  d$y <- 1 + 0.5 * x + u + rnorm(n, 0, exp(-0.2 + 0.3 * z + v))
  m <- plogis(0.2 + 0.4 * x + u)
  d$yb <- rbeta(n, m * 20, (1 - m) * 20)
  d$yc <- rnbinom(n, mu = exp(1 + 0.3 * x + u), size = 3)
  d$yt <- 1 + 0.5 * x + u + 0.8 * rt(n, df = 5)
  lat <- 0.8 * x + u + rlogis(n)
  d$yo <- factor(cut(lat, c(-Inf, -1, 0.5, 2, Inf), labels = FALSE),
                 ordered = TRUE)
  d$y2 <- 0.3 + 0.2 * x + 0.7 * u + 0.5 * (d$y - 1 - 0.5 * x - u) +
    rnorm(n, 0, 0.8)
  d
}
d <- sim_v0(101)
fd <- drmTMB::drmTMB(dbf(mu1 = y ~ x + (1 | p | g),
                         mu2 = y2 ~ x + (1 | p | g), rho12 = ~ 1),
                     family = drmTMB::biv_gaussian(), data = d)
ff <- frm(mvbf(bf(y ~ x + (1 | p | g)), bf(y2 ~ x + (1 | p | g))) +
            set_rescor(TRUE), family = gaussian(), data = d)
# rho_map(cc = 1) on purpose: this script is where the cap was found.
map <- c(lapply(1:6, function(i) lin(i, i)),
         list(rho_map(10, 7), lin(7, 8), lin(8, 9), rho_map(9, 10)))
r <- compare_fits("v0 boundary", fd, ff, map)
fmt_block(r)
eta <- fd$opt$par[[10]]
cat("drm group rho = tanh(eta):", format(tanh(eta), digits = 15),
    "  1 - rho:", format(1 - tanh(eta), digits = 4), "\n")
fz <- ff$opt$par[[9]]
rf <- fz / sqrt(1 + fz^2)
cat("frm group rho = f/sqrt(1+f^2):", format(rf, digits = 15),
    "  1 - rho:", format(1 - rf, digits = 4), "\n")
# Push drmTMB's own objective further along its correlation coordinate,
# re-optimizing everything else, to see whether its optimizer stopped
# early on a flat ridge or its link caps the correlation.
for (e in c(10, 11, 12, 14, 16, 18)) {
  o <- nlminb(fd$opt$par[-10], function(p) fd$obj$fn(append(p, e, 9)))
  cat(sprintf("drm objective, eta fixed at %4.1f: logLik %.10f\n", e,
              -o$objective))
}
cat("drm opt message:", fd$opt$message, " iterations:",
    fd$opt$iterations, "\n")

# Test the hypothesis that drmTMB's group correlation is c * tanh(eta)
# with c = 0.999999: evaluate frmtmb's objective at drm's optimum under
# each candidate map. The map that closes the gap is drmTMB's link.
pd <- fd$opt$par
pf <- ff$opt$par
at_c <- function(cc) {
  q <- pf
  q[1:6] <- pd[1:6]
  q[10] <- sinh(pd[[7]])
  q[7] <- pd[[8]]
  q[8] <- pd[[9]]
  r <- cc * tanh(pd[[10]])
  q[9] <- r / sqrt(1 - r^2)
  (-negll_at(ff$obj, q)) - (-negll_at(fd$obj, pd))
}
for (cc in c(1, 1 - 1e-7, 1 - 1e-6, 1 - 1e-5)) {
  cat(sprintf(paste("c = %.7f: frm objective minus drm objective",
                    "at drm optimum: %.3e\n"), cc, at_c(cc)))
}
ns <- asNamespace("drmTMB")
hit <- Filter(function(nm) {
  f <- get(nm, ns)
  is.function(f) && any(grepl("0.999999", deparse(f), fixed = TRUE))
}, ls(ns))
cat("drmTMB R functions containing the literal 0.999999:", length(hit),
    "\n")
print(head(hit, 20))
