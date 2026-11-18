# rel084-kt copy of dev/nbrobust-repro.R: arm lib = rellib-r23 then user lib (lane lib and rellib-r22 removed); arm base unchanged. No other change.
# Lane nbrobust: the negbinomial fit at its poisson limit, over ulp
# perturbations of the data, on the base build or the lane build.
#   Rscript dev/nbrobust-repro.R <base|lib> <case> [k0] [k1] [cse]
# case nb200: seed 2, n = 200 (backlog item 4903, dev/ci082-probe.R nbk)
# case nb300: seed 4, n = 300 (backlog item 5224,
#   dev/reviews/famend-review-nbpert.R), negbinomial and negbinomial2
# x is scaled by 1 + 2e-15 * k for k = k0..k1 (default 0..20); [cse] is
# frmtmb_control(check_se = ), default "warning".
args <- commandArgs(TRUE)
arm <- args[1]
case <- args[2]
k0 <- if (length(args) >= 3) as.integer(args[3]) else 0L
k1 <- if (length(args) >= 4) as.integer(args[4]) else 20L
cse <- if (length(args) >= 5) args[5] else "warning"
base <- c("C:/Users/adf44/source/r/rellib-r22",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (arm == "lib") c("C:/Users/adf44/source/r/rellib-r23", "C:/Users/adf44/AppData/Local/R/win-library/4.6")
          else base)
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(system.file(package = "frmtmb")), " version",
    as.character(packageVersion("frmtmb")), "\n")
cat("R:", R.home(), " CORETYPE:", Sys.getenv("OPENBLAS_CORETYPE"),
    " threads:", Sys.getenv("OPENBLAS_NUM_THREADS"), "\n")
cat("La:", La_library(), "\n")
conds <- function(expr) {
  w <- list()
  m <- list()
  v <- withCallingHandlers(expr, warning = function(x) {
    w[[length(w) + 1L]] <<- x
    invokeRestart("muffleWarning")
  }, message = function(x) {
    m[[length(m) + 1L]] <<- x
    invokeRestart("muffleMessage")
  })
  list(value = v, warnings = w, messages = m)
}
has <- function(l, cls) any(vapply(l, function(x) inherits(x, cls), NA))
lim <- function(r) has(r$messages, "frmtmb_limit_fit") ||
  has(r$messages, "frmtmb_boundary_fit")
gradw <- function(r) any(grepl("Large maximum absolute gradient",
                               vapply(r$warnings, conditionMessage, "")))
if (case == "nb200") {
  set.seed(2)
  d0 <- data.frame(x = stats::rnorm(200))
  d0$y <- stats::rpois(200, exp(1 + 0.4 * d0$x))
  fams <- list(nb = negbinomial())
} else if (case == "nb300") {
  set.seed(4)
  d0 <- data.frame(x = stats::rnorm(300))
  d0$y <- stats::rpois(300, exp(1 + 0.3 * d0$x))
  fams <- list(nb = negbinomial(), nb2 = brmsfamily("negbinomial2"))
} else stop("case")
g0 <- glm(y ~ x, family = poisson(), data = d0)
cat(sprintf("poisson glm logLik %.12f\n", as.numeric(logLik(g0))))
t0 <- proc.time()
for (k in k0:k1) {
  d <- d0
  d$x <- d$x * (1 + 2e-15 * k)
  for (fn in names(fams)) {
    r <- conds(frm(y ~ x, family = fams[[fn]], data = d,
                   control = frmtmb_control(check_se = cse)))
    fit <- r$value
    p <- fit$opt$par
    cat(sprintf(paste("RESULT %s k %2d %-3s code %d (%s) shape %.4f",
                      "logLik %.12f lim %s gradwarn %s nwarn %d nmsg %d\n"),
                case, k, fn, fit$opt$convergence, fit$opt$message,
                p[length(p)], as.numeric(logLik(fit)), lim(r), gradw(r),
                length(r$warnings), length(r$messages)))
    for (x in r$warnings) cat("   W [", class(x)[1], "]",
                              substr(conditionMessage(x), 1, 220), "\n")
    for (x in r$messages) cat("   M [", class(x)[1], "]",
                              substr(conditionMessage(x), 1, 160), "\n")
  }
}
cat("elapsed", (proc.time() - t0)[["elapsed"]], "\n")
cat("END\n")
