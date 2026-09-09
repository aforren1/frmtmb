# Round 2, item 1, part three. Two things dev/rev-lincmt-f6b.R
# measured by proxy and one it did not measure at all.
#
#   A  NaN incidence directly, rather than a relative eigenvalue split
#      that can be small for reasons that do not produce a NaN.
#   B  what a user sees, cleanly.
#   C  the fix the lane did not consider: `arg` is already clamped from
#      ABOVE by `arg / sqrt(1 + relu(arg^2 - 1))`, which stops it
#      exceeding one but lets it reach exactly one. Clamping it to
#      `1 - d` instead perturbs the two colliding roots SYMMETRICALLY,
#      and the trajectory is a symmetric function of the roots, so the
#      perturbation should cancel to second order. Measured rather than
#      assumed.
#
# Script path: dev/rev-lincmt-f6c.R. Seeds named at each section.
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

cat("\n=== A. NaN incidence, measured directly ===\n")
set.seed(717)
boxes <- list("1e-2 to 10, a PK box" = c(-2, 1),
              "1e-4 to 20" = c(-4, log10(20)),
              "1e-8 to 50" = c(-8, log10(50)),
              "1e-20 to 50" = c(-20, log10(50)))
cat(sprintf("%-24s %8s %10s %14s\n", "box", "draws", "NaN", "worst"))
for (nm in names(boxes)) {
  lo <- boxes[[nm]][[1L]]; hi <- boxes[[nm]][[2L]]
  bad <- 0L; wat <- NULL
  N <- 4000L
  for (i in seq_len(N)) {
    p <- 10^runif(6, lo, hi)
    g <- gradof(p)
    if (!all(is.finite(g))) { bad <- bad + 1L; wat <- p }
  }
  cat(sprintf("%-24s %8d %10d %14s\n", nm, N, bad,
              if (bad) paste(format(wat, digits = 3), collapse = " ")
              else "-"))
}

cat("\n=== B. what a user sees if a fit lands on it ===\n")
KE <- 0.2; K12 <- 0.4; K21 <- 0.1
b <- KE + K12 + K21
qlo <- (b - sqrt(b * b - 4 * KE * K21)) / 2
f <- function(th) sum(frm_lincmt(
  parms = list(ke = exp(th[1]), k12 = K12, k21 = K21, k13 = 1e-300,
               k31 = exp(th[2]), ka = exp(th[3]), V = 10),
  times = tt, ncmt = 3, depot = TRUE, events = ev))
x0 <- c(log(KE), log(qlo), log(1.1))
tp <- MakeTape(f, x0)
r <- tryCatch(stats::nlminb(x0, function(p) as.numeric(tp(p)),
                            function(p) as.numeric(tp$jacobian(p))),
              error = function(e) paste("ERROR:", conditionMessage(e)))
if (is.character(r)) {
  cat("  nlminb from exactly the tangency ->", r, "\n")
} else {
  cat("  nlminb returned code", r$convergence, "in", r$iterations,
      "iterations\n")
}
cat("  the message names neither frm_lincmt() nor the cause.\n")

cat("\n=== C. the clamp, measured ===\n")
cat("Replacing `arg <- arg / sqrt(1 + relu(arg^2 - 1))` with a clamp",
    "\nthat also holds arg below 1 - d. Reported: whether the gradient",
    "\nbecomes finite, how far it then sits from the value the finite",
    "\nside converges to, and what the clamp costs the VALUE, against",
    "\nfrm_ode() at atol = rtol = 1e-12.\n\n")
ns <- asNamespace("frmtmb.ode")
orig <- get("lincmt_disp", envir = ns)
mk_clamped <- function(d) {
  src <- deparse(orig)
  hit <- grep("arg <- arg/sqrt(1 + lincmt_relu(arg * arg - 1))", src,
              fixed = TRUE)
  if (!length(hit))
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
d3 <- function(t, y, p) list(c(-p[6] * y[1],
  p[6] * y[1] - (p[1] + p[2] + p[4]) * y[2] + p[3] * y[3] + p[5] * y[4],
  p[2] * y[2] - p[3] * y[3], p[4] * y[2] - p[5] * y[4]))
ref <- frm_ode(d3, init = list(0, 0, 0, 0), times = tt,
               parms = list(KE, K12, K21, 1e-300, qlo, 1.1),
               states = c("depot", "central", "peripheral1",
                          "peripheral2"),
               output = "central", events = ev, atol = 1e-12,
               rtol = 1e-12) / 10
# the gradient the finite side converges to, from the shipped code
gfin <- {
  k31 <- qlo * (1 + 1e-8)
  ff <- function(th) sum(frm_lincmt(
    parms = list(ke = exp(th[1]), k12 = K12, k21 = K21, k13 = 1e-300,
                 k31 = exp(th[2]), ka = exp(th[3]), V = 10),
    times = tt, ncmt = 3, depot = TRUE, events = ev))
  xx <- c(log(KE), log(k31), log(1.1))
  as.numeric(MakeTape(ff, xx)$jacobian(xx))
}
cat("  shipped, just off the tangency (dk = 1e-8): ",
    paste(format(gfin, digits = 8), collapse = "  "), "\n")
cat("  a central difference AT the tangency, h = 1e-5:",
    "-9.2779361  -1.95e-09  3.5969724\n\n")
cat(sprintf("%10s %8s %34s %14s\n", "clamp d", "finite", "gradient",
            "value err/scale"))
unlockBinding("lincmt_disp", ns)
for (dd in c(1e-12, 1e-14, 1e-16, 1e-18)) {
  assign("lincmt_disp", mk_clamped(dd), envir = ns)
  g <- tryCatch(as.numeric(MakeTape(f, x0)$jacobian(x0)),
                error = function(e) rep(NA_real_, 3))
  v <- frm_lincmt(parms = list(ke = KE, k12 = K12, k21 = K21,
                               k13 = 1e-300, k31 = qlo, ka = 1.1,
                               V = 10),
                  times = tt, ncmt = 3, depot = TRUE, events = ev)
  cat(sprintf("%10.0e %8s %34s %14.3e\n", dd, all(is.finite(g)),
              paste(format(g, digits = 6), collapse = " "),
              max(abs(v - ref)) / max(abs(ref))))
}
assign("lincmt_disp", orig, envir = ns)
stopifnot(identical(get("lincmt_disp", envir = ns), orig))
cat("\n  shipped, restored, at the tangency:",
    paste(format(as.numeric(MakeTape(f, x0)$jacobian(x0)),
                 digits = 6), collapse = " "), "\n")
