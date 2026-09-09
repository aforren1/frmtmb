# The identity has to hold on the tape, not only on doubles.
#
# Three things are checked: that the value taped equals the value
# computed numerically, that the gradient equals a central difference of
# the SAME closed form, and that the gradient equals frm_ode()'s
# gradient through its adjoint solve. The tape node counts of the two
# paths are printed beside them.
source("lincmt-src.R")

pk_dyn <- function(t, y, p) {
  list(c(-p[1] * y[1], p[1] * y[1] / p[3] - p[2] * y[2]))
}

tt <- c(0.5, 1, 2, 4, 6, 8, 12, 18, 24, 30, 36, 48)
doses <- data.frame(time = 0, state = "depot", value = 100, ii = 12,
                    addl = 3L)
ss <- data.frame(time = c(0, 12), state = "depot", value = 100,
                 ii = 12, addl = c(0L, 3L), ss = c(TRUE, FALSE))

lin_f <- function(th, ev, n_ss = 20L) {
  sum(frm_lincmt(parms = list(ka = exp(th[1]), ke = exp(th[2]),
                              V = exp(th[3])),
                 times = tt, ncmt = 1, depot = TRUE, events = ev,
                 n_ss = n_ss))
}
ode_f <- function(th, ev, n_ss = 20L) {
  sum(frm_ode(pk_dyn, init = list(0, 0), times = tt,
              parms = list(exp(th[1]), exp(th[2]), exp(th[3])),
              states = c("depot", "central"), output = "central",
              events = ev, n_ss = n_ss, atol = 1e-12, rtol = 1e-12))
}

fd <- function(f, x, h = 1e-6) {
  vapply(seq_along(x), function(j) {
    xp <- x; xp[j] <- xp[j] + h
    xm <- x; xm[j] <- xm[j] - h
    (f(xp) - f(xm)) / (2 * h)
  }, 0)
}

report <- function(label, ev, th) {
  Fl <- MakeTape(function(p) lin_f(p, ev), th)
  Fo <- MakeTape(function(p) ode_f(p, ev), th)
  gl <- as.numeric(Fl$jacobian(th))
  go <- as.numeric(Fo$jacobian(th))
  gd <- fd(function(p) lin_f(p, ev), th)
  vl <- Fl(th); vo <- Fo(th)
  cat("\n--", label, "--\n")
  cat("value  lincmt", format(vl, digits = 12), " ode",
      format(vo, digits = 12), " rel",
      format(abs(vl - vo) / abs(vo)), "\n")
  cat("grad   lincmt tape ", paste(format(gl, digits = 10),
                                   collapse = "  "), "\n")
  cat("grad   ode    tape ", paste(format(go, digits = 10),
                                   collapse = "  "), "\n")
  cat("grad   lincmt fd   ", paste(format(gd, digits = 10),
                                   collapse = "  "), "\n")
  cat("rel tape-vs-ode   ", format(max(abs(gl - go) / max(abs(go)))),
      "\n")
  cat("rel tape-vs-fd    ", format(max(abs(gl - gd) / max(abs(gd)))),
      "\n")
  nl <- nrow(Fl$data.frame()); no <- nrow(Fo$data.frame())
  cat("tape nodes  lincmt", nl, " ode", no, " ratio",
      format(no / nl, digits = 4), "\n")
  invisible(NULL)
}

th0 <- c(log(1), log(0.2), log(10))
report("addl/ii, ordinary parameters", doses, th0)
report("ss + addl, ordinary parameters", ss, th0)
# the degenerate start: two log rates at the same value is what
# start = list(beta = c(0, 0, ...)) gives
report("addl/ii, ka == ke exactly", doses, c(log(0.2), log(0.2),
                                             log(10)))
report("ss + addl, ka == ke exactly", ss, c(log(0.2), log(0.2),
                                            log(10)))

cat("\n== gradient through ka - ke as it collapses ==\n")
cat(sprintf("%12s %18s %18s %10s\n", "ka - ke", "d/dlka (tape)",
            "d/dlka (fd)", "rel"))
for (d in c(0, 10^-(14:0))) {
  th <- c(log(0.2 + d), log(0.2), log(10))
  F <- MakeTape(function(p) lin_f(p, doses), th)
  g <- as.numeric(F$jacobian(th))[1]
  gd <- fd(function(p) lin_f(p, doses), th, h = 1e-5)[1]
  cat(sprintf("%12.1e %18.10e %18.10e %10.2e\n", d, g, gd,
              abs(g - gd) / abs(gd)))
}

cat("\n== second derivative: the Laplace inner problem needs it ==\n")
th <- c(log(0.2), log(0.2), log(10))
F <- MakeTape(function(p) lin_f(p, ss), th)
H <- F$jacfun()$jacobian(th)
print(round(H, 8))
cat("finite difference of the tape gradient:\n")
G <- function(p) as.numeric(F$jacobian(p))
Hfd <- t(vapply(seq_along(th), function(j) {
  xp <- th; xp[j] <- xp[j] + 1e-5
  xm <- th; xm[j] <- xm[j] - 1e-5
  (G(xp) - G(xm)) / 2e-5
}, numeric(3)))
print(round(Hfd, 8))
cat("max rel diff:", format(max(abs(H - Hfd)) / max(abs(Hfd))), "\n")
