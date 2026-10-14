# Reviewer of lane optstart, round 2: restart_result() after fix round 1.
#   Rscript dev/optstart-rev2-restart.R <lib> <label>
# 1. Every restart_result() call on re1000 (20 seeds) and on the 20
#    loose rel.tol fits of dev/optstart-rev-restart.R: codes, gain
#    against the tolerance, largest move over its allowance, relabel.
# 2. fit$opt and diagnose() on a relabeled fit.
# 3. Direct calls on constructed results: where the move allowance
#    x.tol * max(|p|, u) is or is not a units-free scale.
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat(a[2], "lib", dirname(find.package("frmtmb")), "\n")
`%||%` <- function(a, b) if (is.null(a)) b else a
ns <- asNamespace("frmtmb")
orig <- get("restart_result", ns)
LOG <- list()
cur <- ""
wrap <- function(opt, opt2, control = list(), par_units = NULL) {
  out <- orig(opt, opt2, control, par_units)
  u <- if (length(par_units) == length(opt$par)) par_units else 1
  allow <- 1.5e-8 * pmax(abs(as.numeric(opt$par)), u)
  rt <- ns$restart_rel_tol(control$optimizer %||% "nlminb",
                           control$optCtrl %||% list())
  LOG[[length(LOG) + 1L]] <<- data.frame(case = cur,
    c1 = opt$convergence, c2 = opt2$convergence,
    relabeled = !is.null(out) && !identical(out$convergence,
                                            opt2$convergence),
    gain = opt$objective - opt2$objective,
    tol = rt * max(1, abs(opt$objective)),
    move_over_allow = max(abs(as.numeric(opt2$par) -
                                as.numeric(opt$par)) / allow),
    it2 = opt2$iterations %||% NA)
  out
}
unlockBinding("restart_result", ns)
assign("restart_result", wrap, ns)
fits <- list()
run <- function(label, expr) {
  cur <<- label
  w <- character()
  f <- tryCatch(withCallingHandlers(suppressMessages(expr),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    }), error = function(e) NULL)
  fits[[label]] <<- list(f = f, w = w)
}
for (s in 1:20) {
  set.seed(s)
  d <- data.frame(g = factor(rep(1:20, each = 5)), x = rnorm(100))
  d$y <- 1 + 0.5 * d$x + rnorm(20, 0, 1000)[d$g] + rnorm(100)
  run(paste0("re1000 s", s), frm(y ~ x + (1 | g), data = d))
}
for (s in 1:10) {
  set.seed(100 + s)
  n <- 3000
  d <- data.frame(g = factor(sample(1:60, n, TRUE)), x = rnorm(n),
                  z = rnorm(n))
  d$y <- rpois(n, exp(0.5 + 0.3 * d$x - 0.2 * d$z +
                        rnorm(60, 0, 0.8)[d$g]))
  run(paste0("pois reltol1e-4 s", s),
      frm(y ~ x + z + (1 + x | g), data = d, family = poisson(),
          control = frmtmb_control(optCtrl = list(rel.tol = 1e-4))))
  d$yg <- rnorm(n, 2 + d$x + rnorm(60, 0, 2)[d$g] +
                  rnorm(60, 0, 0.5)[d$g] * d$x, 0.5)
  run(paste0("gauss reltol1e-4 s", s),
      frm(yg ~ x + z + (1 + x | g), data = d,
          control = frmtmb_control(optCtrl = list(rel.tol = 1e-4))))
}
options(width = 200)
L <- do.call(rbind, LOG)
cat("\n## restart_result calls\n")
print(L, digits = 4, row.names = FALSE)
cat("\n## fits: code, warnings\n")
for (nm in names(fits)) {
  f <- fits[[nm]]$f
  cat(sprintf("%-22s code %s nwarn %d restart %s\n", nm,
              if (is.null(f)) "ERR" else f$opt$convergence,
              length(fits[[nm]]$w),
              if (!is.null(f$opt$restart)) {
                paste(f$opt$restart$convergence, f$opt$restart$message,
                      f$opt$restart$iterations)
              } else "-"))
}
rl <- names(fits)[vapply(fits, function(z) !is.null(z$f$opt$restart),
                         TRUE)]
if (length(rl)) {
  f <- fits[[rl[1]]]$f
  cat("\n## fit$opt of", rl[1], "\n")
  str(f$opt[setdiff(names(f$opt), "par")])
  cat("\n## diagnose()\n")
  print(diagnose(f))
}
cat("\n## constructed results (direct calls)\n")
rr <- orig
mk <- function(p, obj, code, msg, it) {
  list(par = p, objective = obj, convergence = code, message = msg,
       iterations = it)
}
show <- function(label, o1, o2, ctl = list(), u = NULL) {
  r <- rr(o1, o2, ctl, u)
  cat(sprintf("%-58s relabeled %s\n", label,
              !identical(r$convergence, o2$convergence)))
}
c0 <- function(p, obj = 1e3) mk(p, obj, 0L, "relative convergence (4)", 30L)
c8 <- function(p, obj = 1e3) mk(p, obj, 1L, "false convergence (8)", 1L)
show("the review's case: rel.tol 1e-4, gain 0.9, move 5",
     c0(c(a = 0), 1e4), c8(c(a = 5), 1e4 - 0.9),
     list(optCtrl = list(rel.tol = 1e-4)))
show("no move, no gain", c0(c(a = 1)), c8(c(a = 1)))
show("log-sd at -30, move 4e-7 (sd changes by 4e-7 relative)",
     c0(c(t = -30)), c8(c(t = -30 + 4e-7)))
show("log-sd at -30, move 6e-7", c0(c(t = -30)), c8(c(t = -30 + 6e-7)))
show("logit at 0, move 1e-8", c0(c(z = 0)), c8(c(z = 1e-8)))
show("coefficient 1e-6 (unit 1e-6), move 1e-8 = 1%, u given",
     c0(c(b = 1e-6)), c8(c(b = 1e-6 + 1e-8)), u = 1e-6)
show("coefficient 1e-6, move 1e-8 = 1%, no par_units (autoscale off)",
     c0(c(b = 1e-6)), c8(c(b = 1e-6 + 1e-8)))
show("caller's x.tol 1e-2, move 1% of a coefficient at 1",
     c0(c(b = 1)), c8(c(b = 1.01)), list(optCtrl = list(x.tol = 1e-2)))
show("optim L-BFGS-B, gain 1e-9 at objective 1",
     c0(c(a = 1), 1), c8(c(a = 1), 1 - 1e-9), list(optimizer = "optim"))
