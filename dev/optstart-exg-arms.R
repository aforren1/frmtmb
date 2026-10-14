# Lane optstart: exgaussian starts for sigma and beta from the
# least-squares residuals, compared through `start =` on the battery's
# exgaussian designs (dev/optstart-battery.R): moment (beta from the
# third moment, clipped to [0.2, 0.9] sqrt(m2), the lane's first rule)
# and equal (sigma = beta = sqrt(m2 / 2), the old rule's ratio at the
# residual spread). mu at least squares in both; the group sd at the
# unit, as the lane's start sets it.
#   Rscript dev/optstart-exg-arms.R <lib or ""> [seeds]
a <- commandArgs(TRUE)
.libPaths(c(if (nzchar(a[1])) a[1], "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
seeds <- if (length(a) > 1) eval(parse(text = a[2])) else 1:10
cat("lib", dirname(find.package("frmtmb")), "\n")
res <- list()
for (des in c("ordinary", "hir2", "offset", "rebig", "offterm")) {
  for (s in seeds) {
    set.seed(s)
    n <- 400
    g <- factor(rep(1:20, length.out = n))
    x <- rnorm(n)
    b0 <- 0.3
    b1 <- 0.5
    sdg <- 0.5
    if (des == "hir2") b1 <- 1e3
    if (des == "offset") b0 <- 1e5
    if (des == "rebig") sdg <- 30
    bg <- rnorm(20, 0, sdg)
    o <- if (des == "offterm") 1e3 * runif(n) else rep(0, n)
    mu <- o + b0 + b1 * x + bg[g]
    y <- frmtmb:::family_registry$exgaussian()$sim(
      list(mu = mu, sigma = rep(1, n), beta = rep(1.5, n)), list(), n)
    d <- data.frame(y = y, x = x, g = g, o = o)
    fo <- if (des == "offterm") y ~ x + offset(o) + (1 | g) else
      y ~ x + (1 | g)
    X <- cbind(1, x)
    q <- qr(X)
    r <- qr.resid(q, y - o)
    m2 <- mean(r^2)
    m3 <- mean(r^3)
    bm <- min(max((max(m3, 0) / 2)^(1 / 3), 0.2 * sqrt(m2)), 0.9 * sqrt(m2))
    th <- unclass(par_template(bf(fo), data = d, family = exgaussian()))$theta
    arms <- list(
      moment = list(beta = unname(qr.coef(q, y - o)),
                    betad = c(`sigma_(Intercept)` = log(sqrt(m2 - bm^2)),
                              `beta_(Intercept)` = log(bm)), theta = th),
      equal = list(beta = unname(qr.coef(q, y - o)),
                   betad = c(`sigma_(Intercept)` = log(sqrt(m2 / 2)),
                             `beta_(Intercept)` = log(sqrt(m2 / 2))),
                   theta = th))
    for (nm in names(arms)) {
      f <- tryCatch(suppressWarnings(frm(fo, data = d, family = exgaussian(),
                                         start = arms[[nm]])),
                    error = function(e) {
                      if (s == 1 && des == "ordinary") print(e)
                      NULL
                    })
      res[[length(res) + 1L]] <- data.frame(
        design = des, seed = s, arm = nm,
        code = if (is.null(f)) NA else f$opt$convergence,
        ll = if (is.null(f)) NA else as.numeric(logLik(f)))
    }
  }
}
X <- do.call(rbind, res)
best <- ave(X$ll, X$design, X$seed, FUN = function(v) max(v, na.rm = TRUE))
X$short <- is.na(X$ll) | X$ll < best - 1e-4
print(aggregate(cbind(error = is.na(X$ll), short = X$short,
                      code = !is.na(X$code) & X$code != 0) ~ design + arm,
                data = X, FUN = sum))
