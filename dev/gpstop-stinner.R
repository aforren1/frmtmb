# Lane gpstop, item 2: inner modes of a student() exact gp() fit. At the
# default's outer optimum and at the best reference start's, the Laplace
# objective from several inner starts: the fit's own mode, zero, the other
# fit's mode, and the field that interpolates the data (the gaussian
# posterior mean of the field at the default's theta and a noise of 1e-4).
# Args: arm, rows as case:curve:n:seed:rho0:sd0 (the reference start).
a <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-gpstop/dev/gpstop-lib.R")
gs_libs(a[1])
rows <- strsplit(a[2], ",")[[1]]
curves <- list(mix = function(x) 2 * sin(3 * x) + 0.5 * x,
               wig = function(x) sin(7 * x) + 0.2 * x,
               two = function(x) sin(x) + 0.4 * sin(10 * x),
               smooth = function(x) sin(x), step = function(x) 2 * (x > 5) - 1)
at_inner <- function(obj, par, u) {
  env <- obj$env
  r <- env$random
  lp <- env$last.par
  lp[r] <- u
  env$last.par <- lp
  -obj$fn(par)
}
for (rw in rows) {
  p <- strsplit(rw, ":")[[1]]
  cs <- p[1]; cn <- p[2]; n <- as.integer(p[3]); seed <- as.integer(p[4])
  rho0 <- as.numeric(p[5]); sd0 <- as.numeric(p[6])
  set.seed(74000 + seed * 31 + n)
  x <- sort(round(runif(n, 0, 10), 4))
  d <- data.frame(x = x, y = curves[[cn]](x) + 0.02 * rt(n, 5))
  fd <- fe_fit(y ~ gp(x), data = d, family = student())$fit
  fr <- fe_fit(y ~ gp(x), data = d, family = student(),
               start = list(theta = c(log(sd0), log(rho0))))$fit
  od <- fd$obj; orf <- fr$obj
  pd <- fd$opt$par; pr <- fr$opt$par
  ud <- od$env$last.par.best[od$env$random]
  ur <- orf$env$last.par.best[orf$env$random]
  # the field that reads the data: a gaussian posterior mean in the tape's
  # own (whitened) coordinates, by least squares on Z at the default's theta
  Z <- as.matrix(fd$frame$linpreds[[1]]$Z)
  cat(sprintf("STI %s %s %d %d default ll %.4f ref ll %.4f nrandom %d\n", cs,
              cn, n, seed, as.numeric(logLik(fd)), as.numeric(logLik(fr)),
              length(ud)))
  for (wh in c("default", "ref")) {
    par <- if (wh == "default") pd else pr
    obj <- if (wh == "default") od else orf
    v <- c(own = at_inner(obj, par, if (wh == "default") ud else ur),
           zero = at_inner(obj, par, 0 * ud),
           other = at_inner(obj, par, if (wh == "default") ur else ud))
    cat(sprintf("STI %s %s %d %d at %s outer: own %.4f zero %.4f other %.4f\n",
                cs, cn, n, seed, wh, v[1], v[2], v[3]))
  }
  cat(sprintf("STIP %s %s %d %d default %s\nSTIP %s %s %d %d ref     %s\n",
              cs, cn, n, seed, paste(signif(pd, 4), collapse = " "), cs, cn, n,
              seed, paste(signif(pr, 4), collapse = " ")))
  flush.console()
}
cat("DONE\n")
