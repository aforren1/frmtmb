# Review of lane quadesc: the NaN-raise-drop chain of the rule's
# stationary escape, traced restart by restart, and every endpoint scored by
# an exact marginal log-likelihood computed here (own skew-normal density in
# base R, own piecewise stats::integrate() split at every observation's
# cliff), independent of the lane's trapezoid and of RTMBdist.
#   Rscript dev/reviews/quadesc-review-chain.R <lib|base> <k...> [iso]
# Data: the ALPHAX fit of tests/testthat/test-quad-escape-record.R (seed
# 1011), response scaled by 1 + 2e-15 * k.
a <- commandArgs(TRUE)
lib <- "C:/Users/adf44/source/r/wt-quadesc-lib"
.libPaths(c(if (a[1] == "lib") lib, "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
ks <- eval(parse(text = a[2]))
iso <- length(a) > 2 && a[3] == "iso"
cat("lib:", dirname(find.package("frmtmb")), "version",
    format(packageVersion("frmtmb")), "BLAS:", extSoftVersion()[["BLAS"]],
    "LAPACK:", La_library(), "\n")
ns <- asNamespace("frmtmb")
st <- new.env()
st$where <- "none"
put <- function(name, f) {
  unlockBinding(name, ns)
  assign(name, f, envir = ns)
}
alpha_of <- function(par) {
  b <- par[names(par) == "betad"]
  paste(format(b[-1], digits = 5), collapse = ",")
}
o_opt <- ns$optimize_obj
put("optimize_obj", function(obj, control, bounds = list(lower = -Inf,
                                                         upper = Inf),
                             par_units = NULL, verbose = 0L,
                             start_par = obj$par, tally = NULL) {
  if (st$where == "rule") {
    st$rid <- st$rid + 1L
    st$hist[[st$rid]] <- list(f = list(), g = list())
  }
  tryCatch(o_opt(obj, control, bounds, par_units, verbose, start_par,
                 tally),
           error = function(e) {
             if (st$where != "none") {
               cat(sprintf("  OPTERR %s restart %d: %s\n", st$where,
                           st$rid, conditionMessage(e)))
             }
             stop(e)
           })
})
o_esc <- ns$escape_stationary
put("escape_stationary", function(obj, opt, frame, control, bounds,
                                  par_units = NULL, verbose = 0L,
                                  tally = NULL) {
  where <- if (is.null(tally)) "laplace" else "rule"
  starts <- ns$stationary_escapes(obj, opt, frame)
  ob <- obj
  if (where == "rule" && length(starts)) {
    st$rule_obj <- obj
    st$rid <- 0L
    st$hist <- list()
    # recorders under whatever the escape wraps around them
    ob$fn <- function(x, ...) {
      v <- obj$fn(x, ...)
      if (st$rid < 1L) return(v)
      h <- st$hist[[st$rid]]
      h$f[[length(h$f) + 1L]] <- list(x = x, v = v)
      st$hist[[st$rid]] <- h
      v
    }
    ob$gr <- function(x, ...) {
      g <- obj$gr(x, ...)
      if (st$rid < 1L) return(g)
      h <- st$hist[[st$rid]]
      h$g[[length(h$g) + 1L]] <- list(x = x, g = g)
      st$hist[[st$rid]] <- h
      g
    }
  }
  if (length(starts)) st$where <- where
  o <- o_esc(ob, opt, frame, control, bounds, par_units, verbose, tally)
  st$where <- "none"
  if (length(starts)) {
    cat(sprintf("  ESC %s starts %d in %.8f out %.8f alpha_in %s alpha_out %s\n",
                where, length(starts), opt$objective, o$objective,
                alpha_of(opt$par), alpha_of(o$par)))
  }
  if (where == "laplace") st$lap_end <- o$par
  if (where == "rule") st$rule_in <- opt$par
  o
})

sn <- function(n, alpha) {
  dl <- alpha / sqrt(1 + alpha^2)
  dl * (abs(stats::rnorm(n)) - sqrt(2 / pi)) +
    sqrt(1 - dl^2) * stats::rnorm(n)
}
set.seed(1011)
stats::rnorm(60 + 15 + 60)
stats::rnorm(50 + 10)
sn(50, 1.5)
d0 <- data.frame(x = stats::rnorm(80), g = factor(rep(1:16, each = 5)))
d0$y <- 0.5 * d0$x + stats::rnorm(16, 0, 1)[d0$g] + sn(80, 2)
fo <- bf(y ~ x + (1 | g), sigma ~ 1, alpha ~ x)

# skew normal with mean m, sd s, shape al (brms's parametrization), log
ldsn <- function(y, m, s, al) {
  dl <- al / sqrt(1 + al^2)
  om <- s / sqrt(1 - 2 * dl^2 / pi)
  xi <- m - om * dl * sqrt(2 / pi)
  z <- (y - xi) / om
  log(2) - log(om) + stats::dnorm(z, log = TRUE) +
    stats::pnorm(al * z, log.p = TRUE)
}
exact <- function(par, d) {
  b <- par[names(par) == "beta"]
  bd <- par[names(par) == "betad"]
  sdg <- exp(par[names(par) == "theta"][1])
  mu <- b[1] + b[2] * d$x
  sig <- exp(bd[1])
  al <- bd[2] + bd[3] * d$x
  L <- 12 * sdg
  tot <- 0
  for (g in levels(d$g)) {
    i <- which(d$g == g)
    lf <- function(u) {
      out <- stats::dnorm(u, 0, sdg, log = TRUE)
      for (j in i) out <- out + ldsn(d$y[j], mu[j] + u, sig, al[j])
      out
    }
    dl <- al[i] / sqrt(1 + al[i]^2)
    om <- sig / sqrt(1 - 2 * dl^2 / pi)
    edge <- d$y[i] - mu[i] + om * dl * sqrt(2 / pi)
    edge <- edge[edge > -L & edge < L]
    br <- sort(unique(c(-L, edge, L)))
    grid <- sort(c(seq(-L, L, length.out = 200001), edge,
                   edge + 1e-9, edge - 1e-9))
    lv <- lf(grid)
    mx <- max(lv[is.finite(lv)])
    s <- 0
    for (j in seq_len(length(br) - 1L)) {
      s <- s + stats::integrate(function(u) exp(lf(u) - mx), br[j],
                                br[j + 1], rel.tol = 1e-11, abs.tol = 0,
                                subdivisions = 2000L)$value
    }
    tot <- tot + mx + log(s)
  }
  tot
}

fitq <- function(d, quad) {
  w <- character()
  f <- withCallingHandlers(
    frm(fo, family = skew_normal(), data = d, quadrature = quad,
        start = list(betad = c(0, 0, 0))),
    warning = function(x) {
      w <<- c(w, paste0(class(x)[1], ": ",
                        substr(conditionMessage(x), 1, 110)))
      invokeRestart("muffleWarning")
    },
    message = function(x) invokeRestart("muffleMessage"))
  list(f = f, w = w)
}
summ_restart <- function(h) {
  fx <- vapply(h$f, function(e) as.numeric(e$v)[1], 1)
  gfin <- vapply(h$g, function(e) all(is.finite(e$g)), NA)
  # the lane's rule: lowest objective among points whose gradient was
  # read right after the objective there, both finite
  best <- NULL
  for (e in h$g) {
    if (!all(is.finite(e$g))) next
    m <- Filter(function(q) identical(unname(as.numeric(q$x)),
                                      unname(as.numeric(e$x))), h$f)
    if (!length(m)) next
    v <- as.numeric(m[[length(m)]]$v)
    if (is.finite(v) && (is.null(best) || v < best$v)) {
      best <- list(x = e$x, v = v)
    }
  }
  first_bad <- which(!gfin)[1]
  list(nf = length(fx), ng = length(gfin), nbad = sum(!gfin),
       best = best, bad = if (!is.na(first_bad)) h$g[[first_bad]])
}

for (k in ks) {
  d <- d0
  d$y <- d$y * (1 + 2e-15 * k)
  cat(sprintf("== k %d\n", k))
  r <- fitq(d, TRUE)
  f <- r$f
  par <- f$opt$par
  cat(sprintf(paste0("FIT k %d ll %.10f alpha %s code %d esc %s lap %s",
                     " exact %.6f\n"),
              k, as.numeric(logLik(f)), alpha_of(par), f$opt$convergence,
              paste(format(f$opt$stationary_escape, digits = 5),
                    collapse = "/"),
              paste(format(f$opt$stationary_escape_laplace, digits = 5),
                    collapse = "/"),
              exact(par, d)))
  for (x in r$w) cat("  WARN", x, "\n")
  if (!is.null(st$hist)) {
    for (j in seq_along(st$hist)) {
      s <- summ_restart(st$hist[[j]])
      cat(sprintf("  RESTART %d nfn %d ngr %d gr_nonfinite %d\n", j, s$nf,
                  s$ng, s$nbad))
      if (!is.null(s$best)) {
        cat(sprintf(paste0("    lowest finite-gradient point: obj %.6f",
                           " alpha %s exact %.6f\n"),
                    s$best$v, alpha_of(s$best$x), exact(s$best$x, d)))
      }
      if (!is.null(s$bad)) {
        cat(sprintf(paste0("    first non-finite gradient at alpha %s,",
                           " components %s\n"),
                    alpha_of(s$bad$x),
                    paste(names(s$bad$x)[!is.finite(s$bad$g)],
                          collapse = ",")))
        st$bad_x <- s$bad$x
      }
    }
  }
  st$hist <- NULL
  if (!is.null(st$lap_end)) {
    cat(sprintf("  LAPESC endpoint alpha %s exact %.6f\n",
                alpha_of(st$lap_end), exact(st$lap_end, d)))
  }
  if (!is.null(st$rule_in)) {
    cat(sprintf("  RULEIN (point the rule's escape started from) alpha %s exact %.6f\n",
                alpha_of(st$rule_in), exact(st$rule_in, d)))
  }
  rl <- fitq(d, FALSE)
  cat(sprintf("LAPFIT k %d ll %.10f alpha %s exact %.6f nwarn %d\n", k,
              as.numeric(logLik(rl$f)), alpha_of(rl$f$opt$par),
              exact(rl$f$opt$par, d), length(rl$w)))
  if (k == ks[1]) {
    cat("par names:", paste(names(par), collapse = " "), "\n")
    cat("theta sd", exp(par[names(par) == "theta"][1]), "VarCorr sd:\n")
    print(VarCorr(f))
  }
  if (iso && !is.null(st$bad_x)) {
    o <- st$rule_obj
    x <- st$bad_x
    v <- o$fn(x)
    g <- o$gr(x)
    cat(sprintf("ISO fn %.8f gr %s\n", v,
                paste(format(g, digits = 4), collapse = " ")))
    fd <- vapply(seq_along(x), function(j) {
      h <- 1e-7 * max(1, abs(x[j]))
      xp <- x; xm <- x
      xp[j] <- xp[j] + h; xm[j] <- xm[j] - h
      (o$fn(xp) - o$fn(xm)) / (2 * h)
    }, 1)
    cat("ISO central differences:", paste(format(fd, digits = 4),
                                          collapse = " "), "\n")
    tpl <- o$env$parameters
    cat("ISO anchor template:\n")
    str(tpl)
    saveRDS(list(x = x, tpl = tpl, d = d, lap_end = st$lap_end,
                 rule_in = st$rule_in),
            sprintf("dev/reviews/quadesc-review-iso-k%d.rds", k))
    iso <- FALSE
  }
  st$lap_end <- NULL
  st$rule_in <- NULL
  st$bad_x <- NULL
}
cat("END\n")
