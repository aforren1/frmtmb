# Reviewer of lane optstart: what restart_result() relabels, and where.
#   Rscript dev/optstart-rev-restart.R <lib> <label>
# On the lane build every call of restart_result() is logged with the
# first run's and the restart's codes, the gain, how far the restart
# moved (in natural units), its iterations, and the gradient and
# headroom at the point the fit reports. The same fits run on rellib-r8
# for the codes and warnings 0.70.0 gave.
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat(a[2], "lib", dirname(find.package("frmtmb")), "\n")
ns <- asNamespace("frmtmb")
LOG <- list()
cur <- ""
if (exists("restart_result", ns)) {
  orig <- get("restart_result", ns)
  wrap <- function(opt, opt2, optCtrl = NULL) {
    out <- orig(opt, opt2, optCtrl)
    pf <- parent.frame()
    u <- pf$par_units %||% 1
    LOG[[length(LOG) + 1L]] <<- data.frame(
      case = cur, c1 = opt$convergence, c2 = opt2$convergence,
      relabeled = !is.null(out) && !identical(out$convergence,
                                              opt2$convergence),
      gain = opt$objective - opt2$objective,
      tol = (optCtrl[["rel.tol"]] %||% 1e-10) * max(1, abs(opt$objective)),
      move = max(abs((opt2$par - opt$par) / u)),
      it2 = opt2$iterations %||% NA, msg2 = opt2$message %||% "")
    out
  }
  `%||%` <- function(a, b) if (is.null(a)) b else a
  environment(wrap) <- environment()
  unlockBinding("restart_result", ns)
  assign("restart_result", wrap, ns)
}
`%||%` <- function(a, b) if (is.null(a)) b else a
res <- list()
run <- function(label, expr) {
  cur <<- label
  w <- character()
  f <- tryCatch(withCallingHandlers(suppressMessages(expr),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    }), error = function(e) e)
  if (inherits(f, "error")) {
    res[[length(res) + 1L]] <<- data.frame(case = label, ll = NA, code = NA,
      msg = substr(conditionMessage(f), 1, 40), gmax = NA, headroom = NA,
      nwarn = NA, warn1 = "")
    return(invisible())
  }
  v <- tryCatch(frmtmb:::grad_verdict(f), error = function(e) NULL)
  hr <- if (!is.null(v) && is.finite(v$proj) && v$proj > 0) {
    g <- drop(f$obj$gr(f$opt$par))
    tryCatch(frmtmb:::grad_headroom(f, g, rep(FALSE, length(g))),
             error = function(e) NA)
  } else NA
  res[[length(res) + 1L]] <<- data.frame(case = label,
    ll = as.numeric(logLik(f)), code = f$opt$convergence,
    msg = substr(f$opt$message %||% "", 1, 28),
    gmax = if (is.null(v)) NA else v$gmax, headroom = hr,
    nwarn = length(w), warn1 = substr(paste(w, collapse = " | "), 1, 60))
}
# 1. a group sd 1000 times the residual sd, 20 seeds (the lane's case)
for (s in 1:20) {
  set.seed(s)
  d <- data.frame(g = factor(rep(1:20, each = 5)), x = rnorm(100))
  d$y <- 1 + 0.5 * d$x + rnorm(20, 0, 1000)[d$g] + rnorm(100)
  run(paste0("re1000 s", s), frm(y ~ x + (1 | g), data = d))
}
# 2. a caller's loose rel.tol: the restart rule's tolerance is rel.tol
# times |objective|, so it loosens with it
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
# 3. a flat direction: two nearly collinear covariates and a group term
for (s in 1:10) {
  set.seed(200 + s)
  n <- 500
  d <- data.frame(g = factor(sample(1:25, n, TRUE)), x1 = rnorm(n))
  d$x2 <- d$x1 + rnorm(n, 0, 1e-3)
  d$y <- 1 + d$x1 + rnorm(25, 0, 1)[d$g] + rnorm(n)
  run(paste0("collinear s", s), frm(y ~ x1 + x2 + (1 | g), data = d))
  # a negative binomial whose shape runs to infinity (Poisson data)
  d$c <- rpois(n, exp(0.3 + 0.5 * d$x1))
  run(paste0("nb-poisson s", s),
      frm(c ~ x1, data = d, family = negbinomial()))
}
L <- if (length(LOG)) do.call(rbind, LOG) else data.frame()
R <- do.call(rbind, res)
options(width = 200)
cat("\n## restart_result calls\n")
print(L, digits = 4, row.names = FALSE)
cat("\n## fits\n")
print(R, digits = 10, row.names = FALSE)
saveRDS(list(L = L, R = R), file.path("dev/optstart-rev-log",
                                      paste0("restart-", a[2], ".rds")))
