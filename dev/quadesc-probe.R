# Lane quadesc: the rule's stationary escape on the three skew normals of
# item 5015 and the CI failure, one BLAS per process, every
# escape_stationary() call logged (where, starts, objective in and out,
# alpha in and out).
#   Rscript dev/quadesc-probe.R <lib|base> <case> [k...]
# case: alphax  test-quad-escape-record.R's fit (seed 1011, data drawn
#               after the two designs before it, as there)
#       symsmall SYMSMALL seed 1 of dev/reviews/quadseam2-review-esc2.R
#       weak     WEAK seed 4 of the same script
# k: the response is scaled by 1 + 2e-15 * k (default 0)
# Each FIT line carries the exact marginal log-likelihood (exact_ll()).
# QUADESC_REF=lap adds the Laplace fit from the same start, scored the
# same way (LAPX line). Unless QUADESC_NOREF=1, two more references: the
# quadrature fit started with alpha off its point (alpha intercept 3),
# and the Laplace fit.
a <- commandArgs(TRUE)
lib <- "C:/Users/adf44/source/r/wt-quadesc-lib"
.libPaths(c(if (a[1] == "lib") lib, "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
case <- a[2]
ks <- if (length(a) > 2) eval(parse(text = a[3])) else 0
refs <- !identical(Sys.getenv("QUADESC_NOREF"), "1")
cat("lib:", dirname(find.package("frmtmb")), "version",
    format(packageVersion("frmtmb")), "\nBLAS:",
    extSoftVersion()[["BLAS"]], "LAPACK:", La_library(), "CORETYPE:",
    Sys.getenv("OPENBLAS_CORETYPE"), "\n")
ns <- asNamespace("frmtmb")
orig <- get("escape_stationary", ns)
unlockBinding("escape_stationary", ns)
alpha_of <- function(par) {
  b <- par[names(par) == "betad"]
  paste(format(b[-1], digits = 4), collapse = ",")
}
assign("escape_stationary", function(obj, opt, frame, control, bounds,
                                     par_units = NULL, verbose = 0L,
                                     tally = NULL) {
  where <- if (is.null(tally)) "laplace" else "rule"
  st <- ns$stationary_escapes(obj, opt, frame)
  o <- orig(obj, opt, frame, control, bounds, par_units, verbose, tally)
  if (length(st)) {
    cat(sprintf("  ESC %s starts %d in %.8f out %.8f alpha_in %s alpha_out %s\n",
                where, length(st), opt$objective, o$objective,
                alpha_of(opt$par), alpha_of(o$par)))
  }
  o
}, envir = ns)
sn <- function(n, alpha) {
  dl <- alpha / sqrt(1 + alpha^2)
  dl * (abs(stats::rnorm(n)) - sqrt(2 / pi)) +
    sqrt(1 - dl^2) * stats::rnorm(n)
}
make <- function(case) {
  if (case == "alphax") {
    set.seed(1011)
    stats::rnorm(60 + 15 + 60)
    stats::rnorm(50 + 10)
    sn(50, 1.5)
    d <- data.frame(x = stats::rnorm(80), g = factor(rep(1:16, each = 5)))
    d$y <- 0.5 * d$x + stats::rnorm(16, 0, 1)[d$g] + sn(80, 2)
    list(d = d, fo = bf(y ~ x + (1 | g), sigma ~ 1, alpha ~ x),
         start = c(0, 0, 0), off = c(0, 3, 0))
  } else if (case == "symsmall") {
    set.seed(1001)
    d <- data.frame(x = rnorm(60), g = factor(rep(1:15, each = 4)))
    d$y <- 0.5 * d$x + rnorm(15, 0, 0.8)[d$g] + rnorm(60)
    list(d = d, fo = bf(y ~ x + (1 | g), sigma ~ 1, alpha ~ 1),
         start = c(0, 0), off = c(0, 3))
  } else if (case == "weak") {
    set.seed(1004)
    stats::rnorm(60 + 15 + 60)
    d <- data.frame(x = rnorm(50), g = factor(rep(1:10, each = 5)))
    d$y <- 0.5 * d$x + rnorm(10, 0, 1.5)[d$g] + sn(50, 1.5)
    list(d = d, fo = bf(y ~ x + (1 | g), sigma ~ 1, alpha ~ 1),
         start = c(0, 0), off = c(0, 3))
  } else stop("unknown case")
}
fit1 <- function(m, d, start, quad) {
  w <- character()
  t0 <- proc.time()[["elapsed"]]
  f <- withCallingHandlers(
    frm(m$fo, family = skew_normal(), data = d, quadrature = quad,
        start = list(betad = start)),
    warning = function(x) {
      w <<- c(w, paste0(class(x)[1], ": ", substr(conditionMessage(x), 1,
                                                  120)))
      invokeRestart("muffleWarning")
    },
    message = function(x) invokeRestart("muffleMessage"))
  list(f = f, w = w, t = proc.time()[["elapsed"]] - t0)
}
# The exact marginal log-likelihood at an outer vector (beta, betad, theta),
# by the trapezoid rule over 20001 points spanning +-12 group sds: no
# quadrature rule, no Laplace approximation.
exact_ll <- function(par, d, ax, n = 20001L) {
  b <- par[names(par) == "beta"]
  bd <- par[names(par) == "betad"]
  sdg <- exp(par[names(par) == "theta"][1])
  mu <- b[1] + b[2] * d$x
  sig <- exp(bd[1])
  al <- if (ax) bd[2] + bd[3] * d$x else rep(bd[2], nrow(d))
  u <- seq(-12 * sdg, 12 * sdg, length.out = n)
  h <- u[2] - u[1]
  tot <- 0
  for (g in levels(d$g)) {
    l <- stats::dnorm(u, 0, sdg, log = TRUE)
    for (i in which(d$g == g)) {
      l <- l + RTMBdist::dskewnorm2(d$y[i], mu[i] + u, sig, al[i],
                                    log = TRUE)
    }
    mx <- max(l)
    w <- exp(l - mx)
    tot <- tot + mx + log(h * (sum(w) - 0.5 * (w[1] + w[n])))
  }
  tot
}
m <- make(case)
ax <- case == "alphax"
for (k in ks) {
  d <- m$d
  d$y <- d$y * (1 + 2e-15 * k)
  r <- fit1(m, d, m$start, TRUE)
  f <- r$f
  esc <- f$opt[["stationary_escape"]]
  lap <- f$opt[["stationary_escape_laplace"]]
  al <- f$estimates$betad[grep("^alpha", names(f$estimates$betad))]
  ll <- as.numeric(logLik(f))
  ex <- exact_ll(f$opt$par, d, ax)
  cat(sprintf(paste0("FIT %s k %d ll %.10f alpha %s code %d esc %s lap %s",
                     " settle %s time %.1f exact %.6f\n"),
              case, k, ll, paste(format(al, digits = 5), collapse = ","),
              f$opt$convergence,
              paste(format(esc, digits = 4), collapse = "/"),
              paste(format(lap, digits = 4), collapse = "/"),
              paste(format(f$opt$quad_settle, digits = 3), collapse = "/"),
              r$t, ex))
  for (x in r$w) cat("  WARN", x, "\n")
  if (Sys.getenv("QUADESC_REF") %in% c("lap", "all")) {
    rl <- fit1(m, d, m$start, FALSE)
    exl <- exact_ll(rl$f$opt$par, d, ax)
    cat(sprintf("LAPX %s k %d lap_ll %.10f lap_exact %.6f margin %.6f\n",
                case, k, as.numeric(logLik(rl$f)), exl, ex - exl))
  }
  if (refs) {
    ro <- fit1(m, d, m$off, TRUE)
    llo <- as.numeric(logLik(ro$f))
    alo <- ro$f$estimates$betad[grep("^alpha", names(ro$f$estimates$betad))]
    cat(sprintf("REFOFF %s k %d ll %.10f alpha %s short %.6f nwarn %d\n",
                case, k, llo, paste(format(alo, digits = 5), collapse = ","),
                llo - ll, length(ro$w)))
    rl <- fit1(m, d, m$start, FALSE)
    all_ <- rl$f$estimates$betad[grep("^alpha",
                                      names(rl$f$estimates$betad))]
    cat(sprintf("REFLAP %s k %d ll %.10f alpha %s nwarn %d\n", case, k,
                as.numeric(logLik(rl$f)),
                paste(format(all_, digits = 5), collapse = ","),
                length(rl$w)))
  }
}
cat("END\n")
