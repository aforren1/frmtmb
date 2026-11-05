# Closed-form truth: the exact gp() on a gaussian response integrates to
#   y ~ N(mu, sd^2 (K + 1e-12 I)[J, J] + diag(sigma^2 or se^2)).
# My own multistart nlminb on that marginal, written here with no frmtmb code,
# is the reference; frm() must reach it. Designs are gpgap's and its review's.
# Seeds in each case. The same script on rellib-r16 shows the failures.
#   Rscript dev/ncgp-truth.R <label> <lib or "base">
args <- commandArgs(TRUE)
label <- args[1]
libs <- c("C:/Users/adf44/source/r/rellib-r16",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (args[2] != "base") libs <- c(args[2], libs)
.libPaths(libs)
suppressPackageStartupMessages(library(frmtmb))
cat("NCGPT lib", dirname(system.file(package = "frmtmb")), "version",
    as.character(packageVersion("frmtmb")), "\n")
REL <- 1e-12
if (nzchar(Sys.getenv("NCGP_PATCH"))) source(Sys.getenv("NCGP_PATCH"))

marginal <- function(p, y, pos2, J, se = NULL, nug = REL) {
  # p = mu, log sd, log rho, (log sigma)
  s2 <- exp(2 * p[2])
  K <- s2 * exp(-pos2 / (2 * exp(2 * p[3])))
  diag(K) <- diag(K) * (1 + nug)
  n <- length(y)
  m <- nrow(K)
  if (is.null(se) && m < n) {
    # replicated positions: the means and the within sum of squares are
    # sufficient, so the solve is over the positions
    sg2 <- exp(2 * p[4])
    nj <- tabulate(J, m)
    ybar <- as.vector(rowsum(y, J)) / nj
    ssw <- sum((y - ybar[J])^2)
    diag(K) <- diag(K) + sg2 / nj
    R <- tryCatch(chol(K), error = function(e) NULL)
    if (is.null(R)) return(1e12)
    r <- backsolve(R, ybar - p[1], transpose = TRUE)
    return(0.5 * sum(r^2) + sum(log(diag(R))) + 0.5 * ssw / sg2 +
             0.5 * (n - m) * log(sg2) + 0.5 * sum(log(nj)) +
             0.5 * n * log(2 * pi))
  }
  V <- K[J, J, drop = FALSE]
  nv <- if (is.null(se)) exp(2 * p[4]) else se^2
  diag(V) <- diag(V) + nv
  R <- tryCatch(chol(V), error = function(e) NULL)
  if (is.null(R)) return(1e12)
  r <- backsolve(R, y - p[1], transpose = TRUE)
  0.5 * sum(r^2) + sum(log(diag(R))) + 0.5 * length(y) * log(2 * pi)
}

closed_max <- function(y, x, se = NULL, extra_start = NULL) {
  u <- sort(unique(x))
  J <- match(x, u)
  D2 <- outer(u, u, "-")^2
  rx <- diff(range(u))
  ry <- max(sd(y), 1e-8 * max(abs(y)), 1e-300)
  noise0 <- log(ry) - c(1, 5)
  rho0 <- log(rx) + c(-3, -1.5, 0)
  grid <- expand.grid(rho = rho0, noise = noise0, off = c(0, 1.5))
  starts <- lapply(seq_len(nrow(grid)), function(i) {
    c(mean(y), log(ry) + grid$off[i], grid$rho[i],
      if (is.null(se)) grid$noise[i])
  })
  if (!is.null(extra_start)) starts <- c(list(extra_start), starts)
  best <- NULL
  scale <- rep(1, length(starts[[1]]))
  scale[1] <- ry
  for (st in starts) {
    fnl <- function(q) marginal(q * scale, y, D2, J, se)
    o <- tryCatch(stats::nlminb(st / scale, fnl,
                                control = list(iter.max = 800, eval.max = 1500,
                                               rel.tol = 1e-12)),
                  error = function(e) NULL)
    if (is.null(o) || !is.finite(o$objective) || o$objective >= 1e12) next
    if (is.null(best) || o$objective < best$objective) best <- o
  }
  # one Nelder-Mead polish of the best, since a marginal at low noise is
  # curved along a ridge and a port stop can sit 1e-3 short of its top
  o2 <- tryCatch(stats::optim(best$par, function(q) marginal(q * scale, y, D2, J, se),
                              method = "Nelder-Mead",
                              control = list(maxit = 2000, reltol = 1e-14)),
                 error = function(e) NULL)
  if (!is.null(o2) && o2$value < best$objective) {
    best$objective <- o2$value
    best$par <- o2$par
  }
  best$par <- best$par * scale
  best
}

rowout <- function(tag, y, x, ff, data, se = NULL, ...) {
  t0 <- Sys.time()
  w <- character()
  fit <- tryCatch(withCallingHandlers(frm(ff, data = data, family = gaussian(), ...),
                                      warning = function(cnd) {
                                        w <<- c(w, class(cnd)[1],
                                                substr(conditionMessage(cnd), 1, 90))
                                        invokeRestart("muffleWarning")
                                      }),
                  error = function(e) e)
  el <- as.numeric(Sys.time() - t0, units = "secs")
  if (inherits(fit, "error")) {
    cat(sprintf("NCGPT %-30s ERROR %s\n", tag, conditionMessage(fit)))
    return(invisible())
  }
  bk <- Filter(function(b) identical(b$covstruct, "gp"), fit$frame$re_blocks)[[1]]
  th <- fit$estimates$theta[bk$theta_idx]
  ll <- -fit$opt$objective
  sg <- if (is.null(se)) as.numeric(stats::sigma(fit))[1] else NA_real_
  # the same frm point read on my marginal at REL
  pf <- c(unname(fit$estimates$beta[1]), th[1], th[2],
          if (is.null(se)) log(sg))
  u <- sort(unique(x))
  atcf <- -marginal(pf, y, outer(u, u, "-")^2, match(x, u), se)
  ccf <- closed_max(y, x, se, extra_start = pf)
  cat(sprintf(paste0("NCGPT %-30s frm ll %.4f sigma %.5g sdgp %.5g rho %.5g conv %d",
                     " %.1fs | cf ll %.4f sigma %.5g sdgp %.5g rho %.5g",
                     " | gap(frm ll - cf ll) %.2e | frm point on cf %.4f",
                     " | warn %s\n"),
              tag, ll, sg, exp(th[1]), exp(th[2]), fit$opt$convergence, el,
              -ccf$objective,
              if (is.null(se)) exp(ccf$par[4]) else NA_real_,
              exp(ccf$par[2]), exp(ccf$par[3]), ll + ccf$objective, atcf,
              if (length(w)) paste(unique(w[seq(1, length(w), 2)]), collapse = ",")
              else "-"))
}

# ---- ordinary noise, curved fields (review M4), seeds 5701, 5702
for (seed in c(5701, 5702)) {
  set.seed(seed)
  n <- 100
  x <- sort(runif(n, 0, 10))
  e0 <- rnorm(n)
  for (nm in c("x^2", "exp(0.8x)", "sin(8x)")) {
    f <- switch(nm, "x^2" = x^2, "exp(0.8x)" = exp(0.8 * x), "sin(8x)" = sin(8 * x))
    y <- f + 0.3 * e0
    rowout(sprintf("%s s=0.3 seed %d", nm, seed), y, x, bf(y ~ gp(x)),
           data.frame(x = x, y = y))
  }
}
# ---- low noise, gpgap's design, seed 20261102
set.seed(20261102)
n <- 60
x <- round(runif(n, 0, 6), 3)
e0 <- rnorm(n)
for (s in c(1e-2, 3e-3, 1e-3, 3e-4, 1e-4)) {
  y <- sin(2 * x) + s * e0
  rowout(sprintf("sin(2x) s=%.0e", s), y, x, bf(y ~ gp(x)),
         data.frame(x = x, y = y))
}
# ---- replicates, 20 positions x 30, seed 131
set.seed(131)
u <- sort(runif(20, 0, 6))
x <- rep(u, each = 30)
e0 <- rnorm(length(x))
for (s in c(2e-2, 1e-2, 3e-4)) {
  y <- sin(2 * x) + s * e0
  rowout(sprintf("rep 20x30 s=%.0e", s), y, x, bf(y ~ gp(x)),
         data.frame(x = x, y = y))
}
# ---- se(): known noise, seed 5601
set.seed(5601)
n <- 60
x <- round(runif(n, 0, 6), 3)
e0 <- rnorm(n)
for (s in c(1e-3, 3e-4)) {
  y <- sin(2 * x) + s * e0
  d <- data.frame(x = x, y = y, s = s)
  rowout(sprintf("se(%.0e)", s), y, x, bf(y | se(s) ~ gp(x)), d, se = rep(s, n))
}
# ---- units, seed 911 (review M7)
set.seed(911)
n <- 60
x <- round(runif(n, 0, 6), 3)
e0 <- rnorm(n)
for (s in c(1e-1, 1e-3)) {
  y <- sin(2 * x) + s * e0
  rowout(sprintf("y s=%.0e", s), y, x, bf(y ~ gp(x)), data.frame(x = x, y = y))
  y1 <- 1000 * y
  rowout(sprintf("y*1000 s=%.0e", s), y1, x, bf(y1 ~ gp(x)),
         data.frame(x = x, y1 = y1))
  y2 <- y / 1000
  rowout(sprintf("y/1000 s=%.0e", s), y2, x, bf(y2 ~ gp(x)),
         data.frame(x = x, y2 = y2))
  x1 <- 1000 * x
  rowout(sprintf("x*1000 s=%.0e", s), y, x1, bf(y ~ gp(x1)),
         data.frame(x1 = x1, y = y))
}
