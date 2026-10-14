# Reviewer of lane optstart: would exgaussian's 0.70.0 divisor (sd / 2
# for each of sigma and beta, applied to the residual rms) keep the
# lane's gains and lose its rebig regression?
#   Rscript dev/optstart-rev-exgfix.R <lib>
# Designs as dev/optstart-battery.R (exgaussian, sigma = beta = 1.5,
# n 400 in 20 groups), seeds 1:10. Arms: the lane's default; "rms2", the
# lane's start with sigma = beta = rms / 2 given through `start`; "old",
# 0.70.0's start given through `start`.
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib", dirname(find.package("frmtmb")), "\n")
reg <- frmtmb:::family_registry
design <- function(des, s) {
  set.seed(s)
  n <- 400; G <- 20
  g <- factor(rep(seq_len(G), length.out = n))
  x <- rnorm(n)
  rs <- 1.5; b0 <- 0.3; b1 <- 0.5 * rs; sdg <- 0.5 * rs
  if (des == "hir2") b1 <- 1e3 * rs
  if (des == "offset") b0 <- 1e5 * rs
  if (des == "rebig") sdg <- 30 * rs
  bg <- rnorm(G, 0, sdg)
  o <- if (des != "offterm") rep(0, n) else 1e3 * rs * runif(n)
  mu <- o + b0 + b1 * x + bg[g]
  y <- reg$exgaussian()$sim(list(mu = mu, sigma = rep(rs, n),
                                 beta = rep(rs, n)), list(), n)
  d <- data.frame(y = y, x = x, g = g, o = o)
  fo <- y ~ x + (1 | g)
  if (des == "offterm") fo <- update(fo, . ~ . + offset(o))
  list(d = d, fo = fo)
}
one <- function(m, start) {
  w <- 0L
  f <- tryCatch(withCallingHandlers(suppressMessages(
    frm(m$fo, data = m$d, family = exgaussian(), start = start)),
    warning = function(x) {
      w <<- w + 1L
      invokeRestart("muffleWarning")
    }), error = function(e) NULL)
  c(ll = if (is.null(f)) NA else as.numeric(logLik(f)),
    code = if (is.null(f)) NA else f$opt$convergence, warned = w)
}
rows <- list()
for (des in c("ordinary", "hir2", "offset", "rebig", "offterm")) {
  for (s in 1:10) {
    m <- design(des, s)
    d <- m$d
    q <- qr(cbind(1, d$x))
    z <- d$y - d$o
    r <- qr.resid(q, z)
    bd <- function(v) c(`sigma_(Intercept)` = log(v),
                        `beta_(Intercept)` = log(v))
    arms <- list(lane = NULL,
                 rms2 = list(beta = unname(qr.coef(q, z)),
                             betad = bd(sqrt(mean(r^2)) / 2)),
                 old = list(beta = c(mean(d$y), 0),
                            betad = bd(sd(d$y) / 2)))
    for (nm in names(arms)) {
      v <- one(m, arms[[nm]])
      rows[[length(rows) + 1L]] <- data.frame(design = des, seed = s,
        arm = nm, ll = v[["ll"]], code = v[["code"]], warned = v[["warned"]])
    }
  }
  cat(des, "done\n")
}
X <- do.call(rbind, rows)
best <- ave(X$ll, X$design, X$seed,
            FUN = function(v) max(c(v, -Inf), na.rm = TRUE))
X$gap <- X$ll - best
print(aggregate(cbind(error = is.na(ll), short = is.na(ll) | gap < -1e-4,
                      nonzero = !is.na(code) & code != 0,
                      warned = warned > 0) ~ design + arm,
                data = X, FUN = sum, na.action = na.pass))
