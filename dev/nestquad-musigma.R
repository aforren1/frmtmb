# (1 | ga) in mu and sigma under the nested rule: convergence of the rule
# and of a plain-R reference as nodes rise. Seed 20261020 (the test data).
.libPaths(c("C:/Users/adf44/source/r/wt-nestquad-lib",
            "C:/Users/adf44/source/r/rellib-r18",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
set.seed(20261020)
ng <- 12; nt <- 6; n <- ng * nt
d <- data.frame(ga = factor(rep(seq_len(ng), each = nt)),
                gb = factor(rep(c("i", "ii"), length.out = n)),
                tt = factor(rep(seq_len(nt), ng)),
                x = rnorm(n), w = runif(n, 0.5, 2))
ab <- interaction(d$ga, d$gb, drop = TRUE)
d$eta <- 0.2 + 0.5 * d$x + rnorm(ng, 0, 0.7)[d$ga] +
  rnorm(nlevels(ab), 0, 1)[ab]
d$y <- rbinom(n, 1, plogis(d$eta))
d$yo <- 1
d$yg <- d$eta + rnorm(n, 0, 0.5)
fo <- bf(yg ~ 1 + x + (1 | ga), sigma ~ (1 | ga)) + gaussian()
fq <- frm(fo, data = d, quadrature = TRUE)
fl <- frm(fo, data = d)
par <- fq$opt$par
fr <- fq$frame
pl <- ns$imp_par_list(fr$par_template, par, fr$map)
print(pl)
print(vapply(fr$re_blocks, function(b) b$term_label, ""))
th <- pl$theta
s <- vapply(fr$re_blocks, function(bk) {
  sqrt(ns$cs_entry(bk$covstruct)$vcov(th[bk$theta_idx], bk)[1, 1])
}, 1)
print(s)
gh <- function(k) {
  i <- seq_len(k - 1); J <- matrix(0, k, k)
  J[cbind(i, i + 1)] <- J[cbind(i + 1, i)] <- sqrt(i)
  e <- eigen(J, symmetric = TRUE); list(x = e$values, lw = log(e$vectors[1, ]^2))
}
lse <- function(v) { m <- max(v); m + log(sum(exp(v - m))) }
b <- unname(pl$beta)
ls0 <- unname(pl$betad[1])
ref <- function(k) {
  g <- gh(k); tot <- 0
  for (a in levels(d$ga)) {
    r <- d$ga == a; mu0 <- b[1] + b[2] * d$x[r]
    lu <- vapply(s[[1]] * g$x, function(u) {
      lse(g$lw + vapply(s[[2]] * g$x, function(v)
        sum(dnorm(d$yg[r], mu0 + u, exp(ls0 + v), log = TRUE)), 1))
    }, 1)
    tot <- tot + lse(g$lw + lu)
  }
  tot
}
for (k in c()) cat(sprintf("REF k %d %.10f\n", k, ref(k)))
lap <- RTMB::MakeADFun(ns$build_objective(fr), fr$par_template, random = "b",
                       map = fr$map, silent = TRUE)
cat(sprintf("LAPLACE %.10f fit %.10f\n", -lap$fn(par),
            as.numeric(logLik(fq))))
plan <- ns$quad_check_blocks(fr)
cat("outer", plan$outer$term_label, plan$outer$dpar %||% "", "inner",
    plan$inner$term_label, "\n")
anc <- ns$quad_nest_anchor(lap, plan, par)
for (k in c()) {
  qo <- ns$quad_nest_objective(fr, plan, anc, c(k, k), pl)
  qo$prime(pl)
  cat(sprintf("NEST k %d %.10f\n", k, -qo$fn(pl)))
}
source("C:/Users/adf44/source/r/frmtmb-wt-nestquad/dev/nestquad-ref.R")
tot <- 0
for (a in levels(d$ga)) {
  r <- d$ga == a; mu0 <- b[1] + b[2] * d$x[r]
  fo_ <- function(u) dnorm(u, 0, s[[1]], log = TRUE) + ref_log_int(function(v)
    dnorm(v, 0, s[[2]], log = TRUE) +
      sum(dnorm(d$yg[r], mu0 + u, exp(ls0 + v), log = TRUE)), 1e-13)
  tot <- tot + ref_log_int(fo_, 1e-10)
}
cat(sprintf("REFINT %.10f\n", tot))
for (k in c()) {
  qo <- ns$quad_nest_objective(fr, plan, anc, c(k, k), pl)
  qo$prime(pl)
  cat(sprintf("NEST k %d %.10f\n", k, -qo$fn(pl)))
}
for (kk in list()) {
  qo <- ns$quad_nest_objective(fr, plan, anc, kk, pl)
  qo$prime(pl)
  cat(sprintf("NEST ko %d ki %d %.10f\n", kk[1], kk[2], -qo$fn(pl)))
}
plan2 <- plan
plan2$outer <- plan$inner
plan2$inner <- plan$outer
plan2$labels <- rev(plan$labels)
stopifnot(identical(plan$row_inner, plan$row_outer))
anc2 <- ns$quad_nest_anchor(lap, plan2, par)
for (kk in list(c(11, 11),
                c(31, 31))) {
  qo <- ns$quad_nest_objective(fr, plan2, anc2, kk, pl)
  qo$prime(pl)
  v <- -qo$fn(pl)
  cat(sprintf("SWAP ko %d ki %d %.10f err %.3e\n", kk[1], kk[2], v,
              abs(v - tot)))
}
ad <- new.env(parent = ns)
sys.source("C:/Users/adf44/source/r/frmtmb-wt-nestquad/dev/nestquad-adapt-version.R.txt", envir = ad)
for (pp in list(list("swap", plan2, anc2))) {
  for (kk in list(c(7, 7), c(11, 11), c(11, 15), c(15, 15), c(21, 21))) {
    qo <- ad$quad_nest_objective(fr, pp[[2]], pp[[3]], kk, pl)
    qo$prime(pl)
    v <- -qo$fn(pl)
    cat(sprintf("ADAPT %s %s ko %d ki %d %.10f err %.3e
", Sys.getenv("NQ_VARIANT"), pp[[1]], kk[1], kk[2],
                v, abs(v - tot)))
  }
}
