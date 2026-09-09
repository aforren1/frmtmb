# Attack 3: the gradient. Three questions.
#
#   1. Does the pinning test actually FAIL against the unfixed
#      primitive, at the ratio the record claims?
#   2. Is the fix complete, or local? lincmt_phi() is one removable
#      singularity; the file has others.
#   3. What in this file is only VALUE tested? A value that is right
#      while its gradient is 1e286 is invisible to a value test.
source("C:/Users/adf44/source/r/frmtmb-wt-lincmt/dev/rev-lincmt-src.R")
ns <- asNamespace("frmtmb.ode")
tt <- c(0, 0.5, 1, 2, 4, 8, 12, 18, 24, 30, 36, 48)

cat("\n=== A. the pinning test, seen to fail ===\n")
ev <- data.frame(time = 0, state = "depot", value = 100, ii = 12,
                 addl = 3L)
f1 <- function(th) sum(frm_lincmt(
  parms = list(ka = exp(th[1]), ke = exp(th[2]), V = exp(th[3])),
  times = tt, ncmt = 1, depot = TRUE, events = ev))
th <- c(log(0.2), log(0.2), log(10))
near <- c(log(0.2 * (1 + 1e-4)), log(0.2), log(10))
report <- function(tag) {
  g <- as.numeric(MakeTape(f1, th)$jacobian(th))
  gn <- as.numeric(MakeTape(f1, near)$jacobian(near))
  cat(sprintf("%-32s value %.12g\n", tag, f1(th)))
  cat("   grad at ka == ke:",
      paste(format(g, digits = 6), collapse = "  "), "\n")
  cat("   the test's assertion max|g| < 2 max|g_near|:",
      if (max(abs(g)) < 2 * max(abs(gn))) "PASS" else "FAIL",
      " ratio", format(max(abs(g)) / max(abs(gn)), digits = 4), "\n")
}
fixed_phi <- get("lincmt_phi", envir = ns)
report("shipped lincmt_phi()")
old_phi <- function(x) { z <- x + 1e-300; -expm1(-z) / z }
environment(old_phi) <- ns
unlockBinding("lincmt_phi", ns)
assign("lincmt_phi", old_phi, envir = ns)
report("offset-only, the first version")
assign("lincmt_phi", fixed_phi, envir = ns)
stopifnot(identical(get("lincmt_phi", envir = ns), fixed_phi))

cat("\n=== B. the same singularity in every response branch ===\n")
# lincmt_phi() is reached through five distinct paths in
# lincmt_resp(). Each is exercised at exact coalescence and its
# gradient compared with a central difference of the same function.
mk <- function(build) function(th) sum(build(exp(th[1]), exp(th[2]),
                                             exp(th[3])))
cases <- list(
  "bolus, ka == ke" = function(ka, ke, V)
    frm_lincmt(parms = list(ka = ka, ke = ke, V = V), times = tt,
               ncmt = 1, depot = TRUE, init = list(depot = 100)),
  "addl block, ka == ke" = function(ka, ke, V)
    frm_lincmt(parms = list(ka = ka, ke = ke, V = V), times = tt,
               ncmt = 1, depot = TRUE, events = ev),
  "ss row, ka == ke" = function(ka, ke, V)
    frm_lincmt(parms = list(ka = ka, ke = ke, V = V), times = tt,
               ncmt = 1, depot = TRUE,
               events = data.frame(time = 0, state = "depot",
                                   value = 100, ii = 12, ss = TRUE)),
  "ss + addl, ka == ke" = function(ka, ke, V)
    frm_lincmt(parms = list(ka = ka, ke = ke, V = V), times = tt,
               ncmt = 1, depot = TRUE,
               events = data.frame(time = c(0, 12), state = "depot",
                                   value = 100, ii = c(12, 12),
                                   addl = c(0L, 3L),
                                   ss = c(TRUE, FALSE))),
  "n_ss = 20, ka == ke" = function(ka, ke, V)
    frm_lincmt(parms = list(ka = ka, ke = ke, V = V), times = tt,
               ncmt = 1, depot = TRUE, n_ss = 20,
               events = data.frame(time = 0, state = "depot",
                                   value = 100, ii = 12, ss = TRUE)),
  "infusion, running, lam u" = function(ka, ke, V)
    frm_lincmt(parms = list(ke = ke, V = V), times = c(0.5, 1, 1.5),
               ncmt = 1, depot = FALSE,
               events = data.frame(time = 0, state = "central",
                                   value = 100, duration = 2)),
  "infusion, finished" = function(ka, ke, V)
    frm_lincmt(parms = list(ke = ke, V = V), times = c(3, 6, 12),
               ncmt = 1, depot = FALSE,
               events = data.frame(time = 0, state = "central",
                                   value = 100, duration = 2)),
  "repeated infusion at ss" = function(ka, ke, V)
    frm_lincmt(parms = list(ke = ke, V = V), times = c(3, 6, 11),
               ncmt = 1, depot = FALSE,
               events = data.frame(time = 0, state = "central",
                                   value = 100, duration = 2,
                                   ii = 12, ss = TRUE)))
fd <- function(f, x, h) vapply(seq_along(x), function(j) {
  xp <- x; xp[[j]] <- xp[[j]] + h
  xm <- x; xm[[j]] <- xm[[j]] - h
  (f(xp) - f(xm)) / (2 * h)
}, 0)
cat(sprintf("%-26s %10s %12s %12s\n", "case", "finite", "max|g|",
            "rel vs FD"))
for (nm in names(cases)) {
  f <- mk(cases[[nm]])
  x <- c(log(0.2), log(0.2), log(10))
  g <- as.numeric(MakeTape(f, x)$jacobian(x))
  g1 <- fd(f, x, 1e-5)
  cat(sprintf("%-26s %10s %12.4g %12.3e\n", nm,
              all(is.finite(g)), max(abs(g)),
              max(abs(g - g1)) / max(abs(g1))))
}

cat("\n=== C. the THREE-compartment eigenvalue collision ===\n")
cat("lincmt_disp()'s cubic goes through acos(). At a double root its",
    "\nargument is +/- 1, where acos has an infinite derivative. The",
    "\nlane's gradient sweep sets every rate EQUAL, which for a",
    "\nmammillary model is NOT a double root of the disposition:",
    "\nke = k12 = k21 = k13 = k31 = r gives r, (2-sqrt3) r and",
    "\n(2+sqrt3) r. So this collision is not in it.\n\n")
KE <- 0.2; K12 <- 0.4; K21 <- 0.1
b <- KE + K12 + K21
qlo <- (b - sqrt(b * b - 4 * KE * K21)) / 2
f3 <- function(th) sum(frm_lincmt(
  parms = list(ke = exp(th[1]), k12 = exp(th[2]), k21 = exp(th[3]),
               k13 = exp(th[4]), k31 = exp(th[5]), ka = exp(th[6]),
               V = 10),
  times = tt, ncmt = 3, depot = TRUE, init = list(depot = 100)))
cat(sprintf("%-30s %8s %14s %14s %12s\n", "k13", "finite", "max|g|",
            "max|g| FD", "rel"))
for (k13 in c(1e-2, 1e-6, 1e-10, 1e-14, 1e-30, 1e-300, 0)) {
  x <- c(log(KE), log(K12), log(K21),
         if (k13 > 0) log(k13) else -745, log(qlo), log(1.1))
  if (k13 == 0) x[[4L]] <- -Inf
  ok <- tryCatch({
    g <- as.numeric(MakeTape(f3, x)$jacobian(x))
    g1 <- fd(f3, x, 1e-5)
    cat(sprintf("%-30.3g %8s %14.6g %14.6g %12.3e\n", k13,
                all(is.finite(g)), max(abs(g)), max(abs(g1)),
                max(abs(g - g1)) / max(abs(g1))))
    TRUE
  }, error = function(e) {
    cat(sprintf("%-30.3g ERROR %s\n", k13,
                substr(conditionMessage(e), 1, 50)))
    FALSE
  })
}

cat("\n=== D. the exact double root, k13 = 0 in the tape's own eyes ===\n")
# k13 enters as a bare rate constant, so a model can hold it at zero
f3b <- function(th) sum(frm_lincmt(
  parms = list(ke = exp(th[1]), k12 = exp(th[2]), k21 = exp(th[3]),
               k13 = th[4] * 0, k31 = exp(th[5]), ka = exp(th[6]),
               V = 10),
  times = tt, ncmt = 3, depot = TRUE, init = list(depot = 100)))
x <- c(log(KE), log(K12), log(K21), 0, log(qlo), log(1.1))
g <- as.numeric(MakeTape(f3b, x)$jacobian(x))
cat("k13 held at exactly 0, k31 on the slow root:\n  grad ",
    paste(format(g, digits = 5), collapse = "  "), "\n  finite:",
    all(is.finite(g)), "\n")
g1 <- fd(f3b, x, 1e-5)
cat("  central difference ",
    paste(format(g1, digits = 5), collapse = "  "), "\n")
cat("  rel", format(max(abs(g - g1)) / max(abs(g1)), digits = 4), "\n")

cat("\n=== E. three compartments, every rate underflowed to zero ===\n")
f3c <- function(th) sum(frm_lincmt(
  parms = list(ke = th[1] * 0, k12 = th[2] * 0, k21 = th[3] * 0,
               k13 = th[4] * 0, k31 = th[5] * 0, ka = exp(th[6]),
               V = 10),
  times = tt, ncmt = 3, depot = TRUE, init = list(depot = 100)))
x <- c(0, 0, 0, 0, 0, log(1.1))
v <- f3c(x)
g <- as.numeric(MakeTape(f3c, x)$jacobian(x))
cat("  value", format(v), " finite:", is.finite(v), "\n")
cat("  grad ", paste(format(g, digits = 5), collapse = "  "), "\n")

cat("\n=== F. lincmt_geo() where lam * ii goes to zero ===\n")
fg <- function(th) sum(frm_lincmt(
  parms = list(ka = exp(th[1]), ke = exp(th[2]), V = exp(th[3])),
  times = c(1, 6, 11), ncmt = 1, depot = TRUE,
  events = data.frame(time = 0, state = "depot", value = 100,
                      ii = 12, ss = TRUE)))
cat(sprintf("%10s %14s %14s %14s\n", "ke", "value", "max|g|",
            "rel vs FD"))
for (ke in c(1e-2, 1e-6, 1e-10, 1e-14, 1e-16, 0)) {
  x <- c(log(1.1), if (ke > 0) log(ke) else -745, log(10))
  v <- fg(x)
  g <- as.numeric(MakeTape(fg, x)$jacobian(x))
  g1 <- fd(fg, x, 1e-4)
  cat(sprintf("%10.1e %14.6g %14.6g %14.3e\n", ke, v, max(abs(g)),
              max(abs(g - g1)) / max(abs(g1))))
}

cat("\n=== G. Hessian at each coalescence ===\n")
hs <- function(f, x) {
  tp <- MakeTape(f, x)
  H <- tp$jacfun()$jacobian(x)
  Hf <- vapply(seq_along(x), function(j) {
    xp <- x; xp[[j]] <- xp[[j]] + 1e-4
    xm <- x; xm[[j]] <- xm[[j]] - 1e-4
    (as.numeric(tp$jacobian(xp)) - as.numeric(tp$jacobian(xm))) / 2e-4
  }, numeric(length(x)))
  c(finite = all(is.finite(H)),
    rel = max(abs(H - Hf)) / max(abs(Hf)))
}
for (nm in names(cases)) {
  f <- mk(cases[[nm]])
  r <- hs(f, c(log(0.2), log(0.2), log(10)))
  cat(sprintf("%-26s finite %s   rel vs FD %.3e\n", nm,
              as.logical(r[["finite"]]), r[["rel"]]))
}
