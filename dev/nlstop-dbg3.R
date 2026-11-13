# Lane nlstop: debug of the first build on items 4960 and 4705.
source("dev/nlstop-common.R")
ns <- asNamespace("frmtmb")
src <- parse("dev/nlstop-fresh.R")
keep <- vapply(src, function(e) {
  is.call(e) && identical(e[[1L]], as.name("<-")) &&
    as.character(e[[2L]])[1] %in% c("F", "pos", "B", "gen", "form")
}, NA)
for (e in src[keep]) eval(e)
d <- gen("bernoulli", "mm", 55 * 1000 + match("mm", names(B)))
r <- suppressMessages(frm(form("bernoulli", "mm"), family = bernoulli(),
                          data = d))
print(tryCatch(ns$nl_limit_message(r), error = function(e) e))
obj <- r$obj; p0 <- r$opt$par; j <- 1:2
g0 <- drop(obj$gr(p0))[j]
H <- matrix(NA, 2, 2)
for (i in 1:2) {
  h <- 1e-4 * abs(p0[j[i]]); up <- p0; dn <- p0
  up[j[i]] <- up[j[i]] + h; dn[j[i]] <- dn[j[i]] - h
  H[, i] <- (drop(obj$gr(up))[j] - drop(obj$gr(dn))[j]) / (2 * h)
}
D <- p0[j]
Hu <- (H + t(H)) / 2 * outer(D, D) + diag(D * g0, 2)
print(g0); print(H); print(Hu); print(eigen(Hu, symmetric = TRUE))
# 4705
set.seed(77)
dn <- data.frame(g = factor(rep(seq_len(8), each = 10)))
dn$x <- rnorm(nrow(dn))
u <- rnorm(8, 0, 0.6)[dn$g]
dn$yn <- 1 + 0.3 * dn$x + exp(0.3 + u)^0.8 + rnorm(nrow(dn), 0, 0.3)
form2 <- bf(yn ~ c0 + exp(a)^k, c0 ~ 1 + x, a ~ 1 + (1 | g), k ~ 1, nl = TRUE)
pr <- set_prior("normal(0.5, 1000)", nlpar = "k")
f <- suppressWarnings(suppressMessages(frm(form2, data = dn, prior = pr)))
nm <- ns$outer_par_names(f)
ia <- match("a_(Intercept)", nm); ik <- match("k_(Intercept)", nm)
it <- match("theta_1", nm)
kt <- 50
opt <- function(par, fn, gr, lower, upper, control) {
  o <- stats::nlminb(par, fn, gr, control = list(iter.max = 1000,
                                                 eval.max = 2000))
  p <- o$par
  s <- p[ik] / kt
  p[ia] <- p[ia] * s; p[ik] <- kt; p[it] <- p[it] + log(abs(s))
  list(par = p, objective = fn(p), convergence = 0L,
       message = "relative convergence (4)")
}
g <- suppressWarnings(suppressMessages(frm(form2, data = dn, prior = pr,
  control = frmtmb_control(optimizer = opt, restarts = 0))))
cat("short", g$opt$objective - f$opt$objective, "\n")
gr <- drop(g$obj$gr(g$opt$par))
cat("grad:", format(gr, digits = 3), "\n")
v <- ns$grad_verdict(g)
str(v[setdiff(names(v), "g")])
N <- ns$headroom_ritz_flat(g)
cat("ritz flat columns:", NCOL(N), "\n")
if (NCOL(N)) print(N)
if (NCOL(N)) print(ns$headroom_prior_gain(g, g$opt$par, N, gr, 1e-3))
