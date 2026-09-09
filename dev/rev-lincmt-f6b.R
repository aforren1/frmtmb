# Round 2, item 1, part two: reachability, what a user sees, and
# whether the NaN has a fix that is not worse than the disease.
#
# Script path: dev/rev-lincmt-f6b.R. Seeds named at each section.
source("C:/Users/adf44/source/r/frmtmb-wt-lincmt/dev/rev-lincmt-src.R")

split_of <- function(ke, k12, k21, k13, k31) {
  d <- frmtmb.ode:::lincmt_disp(3L, ke, k12, k21, k13, k31)
  lam <- sort(vapply(d[["lam"]], as.numeric, 0))
  min(diff(lam)) / max(lam)
}

cat("\n=== A. is the collision reachable at ORDINARY rates? ===\n")
cat("Interlacing says the poles and the zeros of the central",
    "\ncompartment's transfer function strictly separate while k12 and",
    "\nk13 are both positive, so a coalescence forces one of them to",
    "\nzero. Tested rather than asserted: 60000 draws, worst (smallest)",
    "\nrelative eigenvalue split, and what the minimum of k12 and k13",
    "\nwas there.\n\n")
set.seed(6161)
cat(sprintf("%-34s %14s %14s\n", "box", "min rel split",
            "min(k12,k13)"))
boxes <- list(
  "rates 1e-2 to 10, a PK box" = c(-2, 1),
  "rates 1e-4 to 20" = c(-4, log10(20)),
  "rates 1e-8 to 50" = c(-8, log10(50)),
  "rates 1e-16 to 50" = c(-16, log10(50)))
for (nm in names(boxes)) {
  lo <- boxes[[nm]][[1L]]; hi <- boxes[[nm]][[2L]]
  worst <- Inf; wat <- NULL
  for (i in 1:15000) {
    p <- 10^runif(5, lo, hi)
    s <- split_of(p[1], p[2], p[3], p[4], p[5])
    if (is.finite(s) && s < worst) { worst <- s; wat <- p }
  }
  cat(sprintf("%-34s %14.3e %14.3e\n", nm, worst,
              min(wat[[2L]], wat[[4L]])))
}

cat("\n=== B. the k31 window, as a function of k13 ===\n")
cat("k31 swept around the slow root of the reduced quadratic; the",
    "\nhalf-width of the band where the tape's gradient is NaN.\n\n")
KE <- 0.2; K12 <- 0.4; K21 <- 0.1
b <- KE + K12 + K21
qlo <- (b - sqrt(b * b - 4 * KE * K21)) / 2
tt <- c(0.5, 1, 2, 4, 8, 12, 24, 48)
ev <- data.frame(time = 0, state = "depot", value = 100, ii = 8,
                 addl = 3L)
isnan_at <- function(k13, k31) {
  f <- function(th) sum(frm_lincmt(
    parms = list(ke = exp(th[1]), k12 = K12, k21 = K21, k13 = k13,
                 k31 = exp(th[2]), ka = exp(th[3]), V = 10),
    times = tt, ncmt = 3, depot = TRUE, events = ev))
  x <- c(log(KE), log(k31), log(1.1))
  !all(is.finite(as.numeric(MakeTape(f, x)$jacobian(x))))
}
cat(sprintf("%10s %16s %16s\n", "k13", "NaN at the root",
            "widest dk with NaN"))
for (k13 in c(1e-300, 1e-20, 1e-16, 1e-15, 1e-14, 1e-12)) {
  at <- isnan_at(k13, qlo)
  w <- 0
  if (at) for (dk in 10^-(12:5)) if (isnan_at(k13, qlo * (1 + dk)))
    w <- max(w, dk)
  cat(sprintf("%10.0e %16s %16.0e\n", k13, at, w))
}

cat("\n=== C. can a FIT walk into it? ===\n")
cat("Three compartments fitted to data a TWO-compartment model",
    "\ngenerated, so k13 and k31 are unidentified and the optimizer is",
    "\nfree to walk them anywhere. Six starts, seed 31415. Reported:",
    "\nwhere lk13 and lk31 ended, the smallest relative eigenvalue",
    "\nsplit seen at ANY point the objective was evaluated, and whether",
    "\nthe fit produced a NaN.\n\n")
set.seed(31415)
TT <- c(0.25, 0.5, 1, 2, 4, 6, 8, 12)
NS <- 12L
d <- data.frame(id = factor(rep(seq_len(NS), each = length(TT))),
                time = rep(TT, NS), dose = 100)
lke <- log(0.25) + rnorm(NS, 0, 0.2)
lka <- log(1.2) + rnorm(NS, 0, 0.3)
i <- as.integer(d$id)
mu <- numeric(nrow(d))
for (j in seq_len(NS)) {
  k <- which(i == j)
  mu[k] <- frm_lincmt(parms = list(ke = exp(lke[[j]]), k12 = 0.35,
                                   k21 = 0.12, ka = exp(lka[[j]]),
                                   V = 12),
                      times = d$time[k], ncmt = 2, depot = TRUE,
                      init = list(depot = 100))
}
d$conc <- mu + rnorm(nrow(d), 0, 0.08)

seen <- new.env(); seen$min <- Inf; seen$n <- 0L
ns <- asNamespace("frmtmb.ode")
orig_disp <- get("lincmt_disp", envir = ns)
spy <- function(ncmt, ke, k12, k21, k13, k31) {
  r <- orig_disp(ncmt, ke, k12, k21, k13, k31)
  if (ncmt == 3L && !inherits(r[["lam"]][[1L]], "advector")) {
    lam <- sort(vapply(r[["lam"]], as.numeric, 0))
    s <- min(diff(lam)) / max(lam)
    seen$n <- seen$n + 1L
    if (is.finite(s) && s < seen$min) seen$min <- s
  }
  r
}
environment(spy) <- ns
unlockBinding("lincmt_disp", ns)
assign("lincmt_disp", spy, envir = ns)

form <- bf(conc ~ frm_lincmt(parms = list(ke = exp(lke), k12 = exp(lk12),
                                          k21 = exp(lk21),
                                          k13 = exp(lk13),
                                          k31 = exp(lk31),
                                          ka = exp(lka), V = exp(lV)),
                             times = time, group = id, ncmt = 3,
                             depot = TRUE, init = list(depot = dose)),
           lke ~ 1 + (1 | id), lka ~ 1 + (1 | id), lk12 ~ 1,
           lk21 ~ 1, lk13 ~ 1, lk31 ~ 1, lV ~ 1, nl = TRUE)
starts <- list(
  "ordinary" = c(log(0.25), log(1.2), log(0.35), log(0.12),
                 log(0.3), log(0.1), log(12)),
  "k13 already small" = c(log(0.25), log(1.2), log(0.35), log(0.12),
                          -20, log(0.1), log(12)),
  "k13 tiny, k31 on the root" = c(log(0.25), log(1.2), log(0.35),
                                  log(0.12), -30,
                                  log(0.029843788128357585), log(12)),
  "k13 tiny, k31 near the root" = c(log(0.25), log(1.2), log(0.35),
                                    log(0.12), -30,
                                    log(0.0298437 * (1 + 1e-7)),
                                    log(12)),
  "both peripheral rates tiny" = c(log(0.25), log(1.2), log(0.35),
                                   log(0.12), -35, -35, log(12)),
  "k13 at the underflow edge" = c(log(0.25), log(1.2), log(0.35),
                                  log(0.12), -700, log(0.03),
                                  log(12)))
for (nm in names(starts)) {
  seen$min <- Inf; seen$n <- 0L
  r <- tryCatch(frm(form + gaussian(), data = d,
                    start = list(beta = starts[[nm]])),
                error = function(e) conditionMessage(e),
                warning = function(w) {
                  suppressWarnings(frm(form + gaussian(), data = d,
                                       start = list(beta =
                                                      starts[[nm]])))
                })
  if (is.character(r)) {
    cat(sprintf("%-28s ERROR after %d disp calls, min split %.2e\n",
                nm, seen$n, seen$min))
    cat("     ", substr(r, 1, 88), "\n")
  } else {
    fe <- fixef(r)
    e <- if (is.matrix(fe)) fe[, 1L] else unlist(fe)
    cat(sprintf(paste0("%-28s ok   calls %5d  min split %.2e  ",
                       "lk13 %8.2f  lk31 %8.4f\n"), nm, seen$n,
                seen$min, e[["lk13.(Intercept)"]],
                e[["lk31.(Intercept)"]]))
  }
}
assign("lincmt_disp", orig_disp, envir = ns)

cat("\n=== D. what a user sees if a fit DOES land on it ===\n")
f <- function(th) sum(frm_lincmt(
  parms = list(ke = exp(th[1]), k12 = K12, k21 = K21, k13 = 1e-300,
               k31 = exp(th[2]), ka = exp(th[3]), V = 10),
  times = tt, ncmt = 3, depot = TRUE, events = ev))
x0 <- c(log(KE), log(qlo), log(1.1))
tp <- MakeTape(f, x0)
obj <- function(p) as.numeric(tp(p))
gr <- function(p) as.numeric(tp$jacobian(p))
r <- tryCatch(stats::nlminb(x0, obj, gr),
              error = function(e) conditionMessage(e))
if (is.character(r)) cat("  nlminb ERROR:", substr(r, 1, 100), "\n")
else cat("  nlminb returned: code", r$convergence, " message '",
         r$message, "'  objective ", format(r$objective),
         "\n  par ", paste(format(r$par, digits = 8),
                           collapse = "  "), "\n  iterations ",
         r$iterations, "\n")
