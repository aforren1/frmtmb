"""Patch one property out of R/ode.R, for the see-it-fail harness.

Every anchor is cut out of the CURRENT source and the result is
compared with the original, so a patch that matches nothing raises
instead of passing. The first version of this script did not do that:
after the shape changed under it, all eight variants silently ran the
shipped code and reported a clean pass, which is the harness failing
open in exactly the way it exists to catch.
"""
import io
import sys

path, variant = sys.argv[1], sys.argv[2]
before = io.open(path, encoding="utf-8").read()
s = before


def cut(text, first, last):
    a = text.index(first)
    b = text.index(last, a) + len(last)
    return text[a:b]


shape = cut(s, "  h <- 1 - ode_ss_rcap",
            "  list(y = yc + d2 * fac, r = r, f = fac, d = d2)")

COMMON = """  h <- 1 - ode_ss_rcap
  u <- (1 - r) / h
  uc <- lo(hi(u, 1), 0)
  um <- lo(u, 1)
  smooth <- function(x) x^4 * (35 - x * (84 - x * (70 - 20 * x)))
  step <- function(x) x^3 * (35 - x * (84 - x * (70 - 20 * x)))
"""
TAIL = """  gate <- smooth(lo(hi(r / ode_ss_rlow, 1), 0))
  fac <- (r / (um * h)) * step(uc) * gate
  list(y = yc + d2 * fac, r = r, f = fac, d = d2)"""

shapes = {
    "none": COMMON + """  fac <- 0 * r
  list(y = yc, r = r, f = fac, d = d2)""",
    "joint": ("""  r <- rep(sum(d1 * d2) / (sum(d1 * d1) + sum(del)), length(d1))
""" + COMMON + TAIL),
    "nodamp": ("""  r <- (d1 * d2) / (d1 * d1)
""" + COMMON + TAIL),
    "hardcap": COMMON + """  w <- lo(hi((1 - r) / (1 - ode_ss_rcap), 1), 0)
  rr <- lo(hi(r, ode_ss_rcap), 0)
  fac <- w * rr / (1 - rr)
  list(y = yc + d2 * fac, r = r, f = fac, d = d2)""",
    "oldgate": COMMON + """  gate <- step(lo(hi(r / ode_ss_rlow, 1), 0))
  fac <- (r / (um * h)) * step(uc) * gate
  list(y = yc + d2 * fac, r = r, f = fac, d = d2)""",
}

finish = cut(s, "ode_run_in_finish <- function(keep, y, ss_tol",
             "                    ss_extrapolate, atol, rtol)\n}")

finish_outside = '''ode_run_in_finish <- function(keep, y, ss_tol, ss_acc, label,
                              ss_extrapolate, atol, rtol) {
  have <- !vapply(keep, is.null, TRUE)
  last <- y
  three <- have[2L] && have[3L]
  ext <- if (three) ode_ss_extrapolate(keep[[2L]], keep[[3L]], last,
                                       atol, rtol) else NULL
  if (isTRUE(ss_extrapolate) && three) y <- ext[["y"]]
  if (is.null(ss_acc) || !have[3L] || inherits(y, "advector") ||
        inherits(keep[[3L]], "advector")) {
    return(y)
  }
  ode_run_in_report(keep, last, y, have, three, ss_tol, ss_acc, label,
                    ss_extrapolate, atol, rtol)
}'''

report = cut(s, "ode_run_in_report <- function(keep, last, y, have",
             "    stalled <- any(bad & r >= 1)\n  }")

report_old = '''ode_run_in_report <- function(keep, last, y, have, three, ss_tol,
                              ss_acc, label, ss_extrapolate, atol,
                              rtol) {
  scale <- max(1e-12, max(abs(as.numeric(y))))
  move <- abs(as.numeric(last) - as.numeric(keep[[3L]]))
  rel <- max(move) / scale
  stalled <- FALSE
  if (three) {
    ext <- ode_ss_extrapolate(keep[[2L]], keep[[3L]], last, atol, rtol)
    r <- as.numeric(ext[["r"]])
    tail <- abs(as.numeric(ext[["y"]]) - as.numeric(last))
    resid <- if (have[1L]) {
      prev <- ode_ss_extrapolate(keep[[1L]], keep[[2L]], keep[[3L]],
                                 atol, rtol)
      abs(as.numeric(ext[["y"]]) - as.numeric(prev[["y"]]))
    } else tail
    rel <- if (!isTRUE(ss_extrapolate)) max(tail) / scale
           else max(resid) / scale
    stalled <- any(r >= 1) && is.finite(rel) && rel > ss_tol
  }'''

undone_now = """    full <- ifelse(modelled, dd * abs(whole), move)
    undone <- ifelse(modelled, dd * abs(whole - fac), move)"""
undone_r1 = """    full <- dd * abs(fac)
    undone <- 0 * fac"""

if variant in shapes:
    s = s.replace(shape, shapes[variant], 1)
elif variant == "outside":
    s = s.replace(finish, finish_outside, 1)
elif variant == "oldreport":
    s = s.replace(report, report_old, 1)
elif variant == "noundone":
    if undone_now not in s:
        raise SystemExit("anchor drifted: undone block")
    s = s.replace(undone_now, undone_r1, 1)
else:
    raise SystemExit("unknown variant " + variant)

if s == before:
    raise SystemExit("PATCH WAS A NO-OP for " + variant)
io.open(path, "w", encoding="utf-8", newline="\n").write(s)
print("patched " + variant)
