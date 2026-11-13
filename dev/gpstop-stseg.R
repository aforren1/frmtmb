# Lane gpstop, item 2: the student() objective on the segment from the
# default's optimum to the best reference start's, 41 points, inner
# continuation from the default end. Args: arm, rows as
# case:curve:n:seed:rho0:sd0.
a <- commandArgs(TRUE)
source("C:/Users/adf44/source/r/frmtmb-wt-gpstop/dev/gpstop-lib.R")
gs_libs(a[1])
rows <- strsplit(a[2], ",")[[1]]
curves <- list(step = function(x) 2 * (x > 5) - 1)
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
  pd <- fd$opt$par; pr <- fr$opt$par
  obj <- fd$obj
  ts <- seq(0, 1, length.out = 41)
  v <- vapply(ts, function(t) -obj$fn(pd + t * (pr - pd)), 0)
  cat(sprintf("SEG %s %s %d %d t %s\n", cs, cn, n, seed,
              paste(sprintf("%.2f:%.3f", ts, v), collapse = " ")))
  # and back from the ref end
  obj2 <- fr$obj
  v2 <- vapply(rev(ts), function(t) -obj2$fn(pd + t * (pr - pd)), 0)
  cat(sprintf("SEGB %s %s %d %d t %s\n", cs, cn, n, seed,
              paste(sprintf("%.2f:%.3f", rev(ts), v2), collapse = " ")))
  # the gradient at each end, outer
  cat(sprintf("SEGG %s %s %d %d grad default %s | ref %s\n", cs, cn, n, seed,
              paste(signif(fd$obj$gr(pd), 3), collapse = " "),
              paste(signif(fr$obj$gr(pr), 3), collapse = " ")))
  flush.console()
}
cat("DONE\n")
