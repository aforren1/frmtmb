## Reviewer re-check (a): every number and name in the new warning and
## the new diagnose() block, against an INDEPENDENT computation of the
## same quantities from the fit object. Seeds 9301 (bound held + short
## free set) and 9202 (autoscaled, two scales).
## usage: Rscript gradcheck-rev-18-recheck.R <core-lib>
LIB <- commandArgs(TRUE)[1]
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("CORE:", find.package("frmtmb"), "\n\n")
`%||%` <- function(x, y) if (is.null(x)) y else x
gw <- function(expr) {
  w <- character(0)
  val <- withCallingHandlers(expr, warning = function(cond) {
    w <<- c(w, conditionMessage(cond)); invokeRestart("muffleWarning")
  })
  list(fit = val, w = w)
}

## An independent recomputation of the five quantities, using only the
## fit's own tape, the box and base R. Nothing here calls grad_verdict().
indep <- function(f) {
  p <- f$opt$par
  nm <- frmtmb:::outer_par_names(f)
  g <- drop(f$obj$gr(p))
  pu <- f$par_units %||% rep(1, length(p))
  gu <- g * pu
  bx <- frmtmb:::fit_outer_box(f)
  slack <- sqrt(.Machine$double.eps) *
    pmax(1, abs(ifelse(is.finite(bx$lower), bx$lower, 0)),
         abs(ifelse(is.finite(bx$upper), bx$upper, 0)))
  at_lo <- is.finite(bx$lower) & (p - bx$lower) <= slack
  at_hi <- is.finite(bx$upper) & (bx$upper - p) <= slack
  act <- (at_lo & g > 0) | (at_hi & g < 0)
  gp <- gu; gp[act] <- 0
  free <- which(!act)
  H <- tryCatch(f$obj$he(p), error = function(e) NULL)
  if (is.null(H) || !identical(dim(H), c(length(p), length(p)))) {
    H <- tryCatch(stats::optimHess(p,
                                   function(q) as.numeric(f$obj$fn(q)),
                                   function(q) drop(f$obj$gr(q))),
                  error = function(e) NULL)
  }
  hr <- NA_real_
  if (!is.null(H)) {
    H <- (H + t(H)) / 2
    ch <- tryCatch(chol(H[free, free, drop = FALSE]),
                   error = function(e) NULL)
    gf <- g; gf[act] <- 0
    if (!is.null(ch)) {
      hr <- 0.5 * sum(backsolve(ch, gf[free], transpose = TRUE)^2)
    }
  }
  list(gmax = max(abs(gu)), gmax_par = nm[which.max(abs(gu))],
       held = nm[act],
       proj = max(abs(gp)),
       proj_par = if (all(act)) NA_character_ else nm[which.max(abs(gp))],
       headroom = hr, raw_max = max(abs(g)),
       raw_max_par = nm[which.max(abs(g))], par_units = pu)
}

check <- function(tag, r) {
  f <- r$fit
  i <- indep(f)
  d <- diagnose(f, quiet = TRUE)
  cat("=================== ", tag, "\n", sep = "")
  cat("INDEPENDENT:\n")
  cat("  natural-unit max|g| =", format(i$gmax, digits = 10),
      "at", i$gmax_par, "\n")
  cat("  raw        max|g|   =", format(i$raw_max, digits = 10),
      "at", i$raw_max_par, "\n")
  cat("  par_units           =",
      paste(format(i$par_units, digits = 6), collapse = " "), "\n")
  cat("  bound-held          = [", paste(i$held, collapse = ","), "]\n")
  cat("  projected max       =", format(i$proj, digits = 10),
      "at", i$proj_par, "\n")
  cat("  headroom            =", format(i$headroom, digits = 10), "\n")
  cat("PACKAGE (diagnose):\n")
  cat("  max_grad =", format(d$max_grad, digits = 10),
      "at", d$worst_grad, "\n")
  cat("  grad_proj =", format(d$grad_proj %||% NA, digits = 10),
      "at", d$grad_proj_par %||% NA, "\n")
  cat("  grad_bound_held = [",
      paste(d$grad_bound_held %||% character(0), collapse = ","), "]\n")
  cat("  grad_headroom =", format(d$grad_headroom %||% NA, digits = 10),
      "\n")
  ## the assertions, as identities where they should be identities
  cat("AGREEMENT:\n")
  cat("  diagnose max_grad == raw max      :",
      identical(d$max_grad, i$raw_max), "\n")
  cat("  diagnose worst_grad == raw argmax :",
      identical(d$worst_grad, i$raw_max_par), "\n")
  cat("  diagnose grad_proj == indep proj  :",
      identical(d$grad_proj, i$proj), "\n")
  cat("  diagnose grad_proj_par == indep   :",
      identical(d$grad_proj_par, i$proj_par), "\n")
  cat("  diagnose held == indep held       :",
      identical(d$grad_bound_held, i$held), "\n")
  hd <- d$grad_headroom %||% NA_real_
  cat("  headroom identical                :",
      identical(hd, i$headroom),
      " rel gap",
      if (is.finite(hd) && is.finite(i$headroom) && i$headroom != 0) {
        format(abs(hd / i$headroom - 1), digits = 4)
      } else "NA", "\n")
  cat("WARNING TEXT:\n")
  if (!length(r$w)) cat("  (no warning)\n")
  for (w in r$w) cat("  ", w, "\n", sep = "")
  cat("  prefix verbatim:",
      any(grepl("^Large maximum absolute gradient at the optimum",
                r$w)), "\n")
  cat("PRINTED diagnose():\n")
  writeLines(paste0("  ", capture.output(diagnose(f))))
  cat("\n")
  invisible(list(i = i, d = d))
}

## ---- seed 9301: a bound holds x, the free components are short
set.seed(9301)
n <- 3000
d1 <- data.frame(x = rnorm(n), z = rnorm(n))
d1$y <- rpois(n, exp(0.4 + 0.9 * d1$x - 0.6 * d1$z))
r1 <- gw(frm(bf(y ~ x + z), family = poisson(), data = d1,
             prior = set_prior("", class = "b", coef = "x", ub = 0.1),
             control = frmtmb_control(
               restarts = 0,
               optCtrl = list(rel.tol = 1e-2, x.tol = 1e-2,
                              iter.max = 2000, eval.max = 2000))))
check("seed 9301: ub on x binds, free set short", r1)

## ---- seed 9202: autoscaled, so raw and natural units differ
set.seed(9202)
n <- 3000
d2 <- data.frame(xs = rnorm(n) * 1e-5, z = rnorm(n))
d2$y <- rnorm(n, 1 + 2e5 * d2$xs + 0.5 * d2$z, 1)
r2 <- gw(frm(bf(y ~ xs + z), family = gaussian(), data = d2,
             prior = set_prior("", class = "b", coef = "z", ub = 0.1),
             control = frmtmb_control(
               restarts = 0,
               optCtrl = list(rel.tol = 1e-3, x.tol = 1e-3,
                              iter.max = 1000, eval.max = 1000))))
check("seed 9202: autoscaled AND ub on z", r2)

## ---- B1 again, since proj_par is new: every component bound-held, so
## proj_par must be NA rather than a wrong name
set.seed(9103)
n <- 600
d3 <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
d3$y <- rnorm(n, 1 + 2 * d3$x1 - 2 * d3$x2, 1)
r3 <- gw(frm(bf(y ~ x1 + x2), family = gaussian(), data = d3,
             prior = c(set_prior("", class = "b", coef = "x1", ub = 0.1),
                       set_prior("", class = "b", coef = "x2", lb = -0.1))))
check("seed 9103: two bounds bind, interior sigma", r3)

cat("DONE recheck-a\n")
