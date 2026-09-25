# Reviewer 2, item 2 (derivatives): tape gradient and tape Hessian of the
# lane's log distribution functions against finite differences, at the
# branch edges: the Taylor branch of ddm_lsinhc_s() (|x| 0.035 to 0.08),
# the tanh clamp at 40 in each blend (u = u0 exp(+-40 * 0.12)), the
# time hold at u = 1e-10, and zero drift.
# Parameters p = (v, a, w, t), t the decision time.
# Gradient check: Richardson central difference of the double function.
# Hessian check: Richardson central difference of the TAPE gradient.
source("dev/phase3b-review2/r2-prelude.R")
r2_lib("lane")
suppressPackageStartupMessages({library(RTMB); library(frmtmb.eam)})
ns <- asNamespace("frmtmb.eam")

outs <- list(
  lS = function(p) ns$ddm_rt_lcdf2(p[4], p[1], p[2], p[3])[["lS"]],
  lF = function(p) ns$ddm_rt_lcdf2(p[4], p[1], p[2], p[3])[["lF"]],
  lFl = function(p) ns$ddm_rt_lcdf_b(p[4], p[1], p[2], p[3], 0),
  lFu = function(p) ns$ddm_rt_lcdf_b(p[4], p[1], p[2], p[3], 1),
  lint = function(p) ns$ddm_rt_linterval_b(p[4], p[4] * 1.3, p[1], p[2],
                                           p[3], 0))
rich <- function(f, x, i, h) {
  e <- replace(numeric(length(x)), i, 1)
  D <- function(h) (f(x + h * e) - f(x - h * e)) / (2 * h)
  (4 * D(h / 2) - D(h)) / 3
}
check_one <- function(nm, p) {
  f <- outs[[nm]]
  tp <- tryCatch(MakeTape(function(q) f(q), p), error = function(e) NULL)
  if (is.null(tp)) return(data.frame(out = nm, err = "tape failed"))
  val <- f(p); tval <- tp(p)
  g <- as.numeric(tp$jacobian(p))
  gf <- tp$jacfun()
  H <- gf$jacobian(p)
  # step relative to each coordinate; v = 0 gets an absolute step
  h <- ifelse(p == 0, 1e-6, 1e-5 * abs(p))
  gfd <- vapply(1:4, function(i) rich(f, p, i, h[i]), 0)
  gfd2 <- vapply(1:4, function(i) rich(f, p, i, 10 * h[i]), 0)
  gg <- function(q) as.numeric(gf(q))
  Hfd <- sapply(1:4, function(i) rich(gg, p, i, h[i]))
  Hfd2 <- sapply(1:4, function(i) rich(gg, p, i, 10 * h[i]))
  # an error counts only past the finite difference's own noise, read
  # from the disagreement of two step sizes, and past 1e-10 of the scale
  rel <- function(a, b, b2) {
    noise <- abs(b - b2)
    ex <- pmax(abs(a - b) - 10 * noise, 0)
    # relative for entries above one, absolute below: an objective is a
    # sum of O(1) log terms, so an absolute 1e-10 error there is inert
    max(ex / pmax(abs(b), 1))
  }
  data.frame(out = nm, v = p[1], a = p[2], w = p[3], t = p[4], val = val,
             tape_val_diff = abs(tval - val),
             grad_rel = rel(g, gfd, gfd2), hess_rel = rel(H, Hfd, Hfd2),
             grad_nan = any(!is.finite(g)), hess_nan = any(!is.finite(H)),
             maxabsH = max(abs(H)))
}

P <- list()
add <- function(v, a, w, t, tag) P[[length(P) + 1L]] <<- list(p = c(v, a, w, t), tag = tag)
# Taylor edges on x = v a, with a = 1.4, u = 0.3 (large-time route live)
for (x in c(0, 1e-9, 0.02, 0.035, 0.045, 0.05, 0.055, 0.07, 0.0806, 0.1, -0.05)) {
  add(x / 1.4, 1.4, 0.5, 0.3 * 1.96, sprintf("taylor va=%g", x))
  add(x / 1.4, 1.4, 0.3, 2 * 1.96, sprintf("taylor va=%g w=.3 u=2", x))
}
# tanh clamp edges: F centre 0.1, S 0.05, B 0.2, scale 0.12, clamp 40
for (u0 in c(0.1, 0.05, 0.2)) for (sgn in c(-1, 1)) for (eps in c(-0.02, 0, 0.02)) {
  u <- u0 * exp(sgn * 40 * 0.12 + eps)
  add(0.7, 1.4, 0.5, u * 1.96, sprintf("clamp u0=%g %+d eps=%g", u0, sgn, eps))
  add(-2, 0.8, 0.2, u * 0.64, sprintf("clamp u0=%g %+d eps=%g b", u0, sgn, eps))
}
# time hold at u = 1e-10
for (m in c(0.5, 0.999, 1.001, 2, 1e3, 1e5)) {
  add(0.7, 1.4, 0.5, m * 1e-10 * 1.96, sprintf("hold u=%g e-10", m))
}
# near-zero drift in general positions
for (v in c(0, 1e-8, -1e-8, 1e-4)) for (u in c(0.01, 0.2, 1, 5)) {
  add(v, 1.4, 0.4, u * 1.96, sprintf("v0 v=%g u=%g", v, u))
}
set.seed(7)
for (i in 1:40) {
  a <- exp(runif(1, log(0.3), log(4))); u <- exp(runif(1, log(1e-3), log(10)))
  add(runif(1, -8, 8), a, runif(1, 0.05, 0.95), u * a^2, "random")
}

rows <- list()
for (pp in P) for (nm in names(outs)) {
  r <- check_one(nm, pp$p); r$tag <- pp$tag; rows[[length(rows) + 1L]] <- r
}
res <- do.call(rbind, lapply(rows, function(r) {
  if (is.null(r$grad_rel)) return(NULL); r
}))
saveRDS(res, "dev/phase3b-review2/deriv.rds")
options(width = 200)
cat("rows:", nrow(res), "\n")
cat("NaN gradient rows:", sum(res$grad_nan), " NaN Hessian rows:", sum(res$hess_nan), "\n")
res$grp <- sub(" .*", "", res$tag)
print(aggregate(cbind(grad_rel, hess_rel, tape_val_diff) ~ grp + out, res, max))
cat("\nworst 15 Hessian rows\n")
print(head(res[order(-res$hess_rel), c("tag", "out", "v", "a", "w", "t", "val",
                                       "grad_rel", "hess_rel", "maxabsH")], 15))
cat("\nworst 10 gradient rows\n")
print(head(res[order(-res$grad_rel), c("tag", "out", "v", "a", "w", "t", "val",
                                       "grad_rel", "hess_rel")], 10))
