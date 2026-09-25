# Joint law of simulate(newdata = ) under a residual correlation term.
#   Rscript dev/simnewdata-review/rv-autocor.R > .../log/autocor.txt
# The exact covariance of two rows at newdata is sigma_i sigma_j R[t_i,t_j]
# within a group and 0 across groups. Each check reports the largest
# |z| of the empirical covariance/variance against that, over nsim draws.
source("dev/simnewdata-review/rv-prelude.R")
set.seed(8)
n_g <- 60
d <- expand.grid(time = 1:6, g = factor(seq_len(n_g)))
d$y <- 1 + unlist(lapply(seq_len(n_g), function(i) {
  as.vector(stats::arima.sim(list(ar = 0.7), 6))
}))
fit <- frm(bf(y ~ 1 + ar(time, gr = g, cov = TRUE)), data = d)
R <- autocor_matrix(fit)
sg <- sigma(fit)
cat(sprintf("fitted phi-lag1 %.4f sigma %.4f\n", R[1, 2], sg))

nsim <- 20000
check <- function(label, nd, fit, R, sg) {
  s <- tryCatch(simulate(fit, nsim = nsim, seed = 3, newdata = nd),
                error = function(e) e)
  if (inherits(s, "error")) {
    cat(label, ": ERROR:", conditionMessage(s), "\n"); return(invisible())
  }
  m <- as.matrix(s)
  n <- nrow(nd)
  tl <- as.numeric(rownames(R))
  ti <- match(as.numeric(nd$time), tl)
  S <- matrix(0, n, n)
  for (i in 1:n) for (j in 1:n) {
    if (nd$g[i] == nd$g[j]) S[i, j] <- sg^2 * R[ti[i], ti[j]]
  }
  E <- stats::cov(t(m))
  # sd of a sample covariance: sqrt((S_ii S_jj + S_ij^2) / nsim)
  sdc <- sqrt((outer(diag(S), diag(S)) + S^2) / nsim)
  z <- (E - S) / sdc
  mz <- (rowMeans(m) - stats::coef(fit)[[1]]) / (sqrt(diag(S)) / sqrt(nsim))
  cat(sprintf("%-40s rows %d  max|z| cov %.2f  max|z| mean %.2f\n",
              label, n, max(abs(z)), max(abs(mz))))
  invisible(list(E = E, S = S))
}

# (a) sorted, one group
check("a sorted 1 group t1..6",
      data.frame(time = 1:6, g = factor("A")), fit, R, sg)
# (b) scrambled time order within a group, interleaved with a second
# group, and a single-row group
nd <- data.frame(time = c(4, 1, 6, 2, 5, 2, 3, 3, 5),
                 g = factor(c("A", "A", "B", "A", "A", "B", "A", "C", "B")))
check("b scrambled+interleaved+single", nd, fit, R, sg)
# (c) a group with gaps: times 1 and 5 only (lag 4)
check("c gap lag 4", data.frame(time = c(5, 1), g = factor(c("A", "A"))),
      fit, R, sg)
# (d) new groups only; group labels that are existing fitted levels
check("d fitted group labels",
      data.frame(time = c(3, 2, 1), g = factor(c("5", "5", "5"))), fit, R, sg)
# (e) duplicated time within a group
cat("e dup time:", try_msg(simulate(fit, nsim = 2, newdata =
  data.frame(time = c(2, 2), g = factor(c("A", "A"))))), "\n")
# (f) no gr: whole newdata one group
fit0 <- frm(bf(y ~ 1 + ar(time, cov = TRUE)),
            data = d[d$g == "1", ])
cat("f fit without gr: d =", nrow(autocor_matrix(fit0)), "\n")
check("f no gr", data.frame(time = c(3, 1, 2), g = factor("x")),
      fit0, autocor_matrix(fit0), sigma(fit0))

# (g) the row order of the output is the row order of newdata:
# a row at time 1 and a row at time 6 of the same group, distinct
# means via a covariate, so a swap would show in the means
d$x <- stats::rnorm(nrow(d))
d$y2 <- d$y + 3 * d$x
fitx <- frm(bf(y2 ~ x + ar(time, gr = g, cov = TRUE)), data = d)
ndx <- data.frame(time = c(6, 1, 3), g = factor(c("A", "A", "B")),
                  x = c(5, -5, 0))
sx <- as.matrix(simulate(fitx, nsim = 4000, seed = 1, newdata = ndx))
cat("g row means", format(rowMeans(sx), digits = 4), " expected ",
    format(as.vector(frm_linpred(fitx, newdata = ndx)), digits = 4), "\n")

# (h) sigma varying by row: sigma ~ x
d$y3 <- d$y + (exp(0.5 * d$x) - 1) * stats::rnorm(nrow(d))
fits <- frm(bf(y3 ~ 1 + ar(time, gr = g, cov = TRUE), sigma ~ x),
            data = d)
nds <- data.frame(time = c(2, 3, 4), g = factor("A"), x = c(-1, 0, 1.5))
ss <- as.matrix(simulate(fits, nsim = nsim, seed = 5, newdata = nds))
sgr <- as.vector(frm_linpred(fits, newdata = nds, dpar = "sigma",
                             type = "response"))
Rs <- autocor_matrix(fits)
S <- outer(sgr, sgr) * Rs[2:4, 2:4]
cat("h sigma~x: emp cov\n"); print(round(stats::cov(t(ss)), 4))
cat("   exact\n"); print(round(S, 4))

# (i) rows of the fit dropped by NA: time 1 has NA response everywhere,
# so the fit's time levels are 2..6 (d = 5). unstr so that a level shift
# shows (ar is Toeplitz and would hide it).
set.seed(2)
du <- expand.grid(time = 1:4, g = factor(1:80))
Lc <- t(chol(matrix(c(1, .8, .1, .5,
                      .8, 1, .3, .2,
                      .1, .3, 1, -.4,
                      .5, .2, -.4, 1), 4)))
du$y <- as.vector(Lc %*% matrix(rnorm(4 * 80), 4))
du <- rbind(data.frame(time = 0, g = factor(1:80), y = NA), du)
fitu <- frm(bf(y ~ 1 + unstr(time, gr = g)), data = du)
Ru <- autocor_matrix(fitu)
cat("i unstr fitted time levels:", rownames(Ru), " nobs", nobs(fitu),
    " data_frame rows", nrow(fitu$frame$data_frame), "\n")
print(round(Ru, 3))
ndu <- data.frame(time = c(1, 2, 3, 4), g = factor("A"))
r <- tryCatch(as.matrix(simulate(fitu, nsim = nsim, seed = 2, newdata = ndu)),
              error = function(e) conditionMessage(e))
if (is.character(r)) cat("i ERROR:", r, "\n") else {
  cat("i emp cor at times 1..4\n"); print(round(stats::cor(t(r)), 3))
}
ndu0 <- data.frame(time = c(0, 1), g = factor("A"))
cat("i time 0 (only on NA rows):",
    try_msg(simulate(fitu, nsim = 2, newdata = ndu0)), "\n")
