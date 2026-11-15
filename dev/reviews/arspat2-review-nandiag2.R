# Review of lane arspat2: where the fcor() margin route's gradient turns
# NaN (n = 400, m = 300, the failing cell of arspat2-review-nan.R). Runs
# nlminb by hand on the dry-run objective and stops at the first
# non-finite fn or gradient, printing the outer point and the pieces.
# The dry run skips mi_ml_frame(); mi_sampling(, "laplace") applies it
# and retapes, which gives the ML route frm() fits.
source("dev/reviews/arspat2-review-common.R")
k <- 20; m <- 300
W <- rook(k); W <- W / rowSums(W)
n <- nrow(W)
A0 <- diag(n) - 0.4 * W
M <- solve(crossprod(A0))
set.seed(11)
d <- data.frame(x = rnorm(n), h = factor(rep(1:8, length.out = n)))
d$y <- as.vector(solve(A0, 1 + 0.5 * d$x + rnorm(8)[d$h] +
                         rnorm(n) * sqrt(4 / rchisq(1, 4))))
set.seed(m)
d$y[sort(sample(n, m))] <- NA
pr <- set_prior("gamma(2, 0.1)", class = "nu")
ob <- quiet(frm(y | mi() ~ x + (1 | h) + fcor(M), data = d,
                data2 = list(M = M), family = student(), prior = pr,
                dry_run = "objective"))
ob <- ns$mi_sampling(ob, "laplace")
tm <- ns$autocor_seam_t_missing(ob$frame)
cat("route margin:", tm$y$margin, "m", tm$y$m, "
")
obj <- ob$obj
# The joint objective (no Laplace) and the inner Hessian at the start,
# then after one fn(): is the NaN the joint density or the Laplace step?
env <- obj$env
rnd <- env$random
cat("joint nll at start:", format(env$f(env$par)), "\n")
pk <- function() sprintf("%.2f GB", ps::ps_memory_info()[["peak_wset"]] / 2^30)
cat("peak before Hessian:", pk(), "
")
H <- tryCatch(as.matrix(env$spHess(env$par, random = TRUE)), error = function(e) {
  cat("spHess error:", conditionMessage(e), "; peak", pk(), "
"); NULL })
if (is.null(H)) { cat("fn at start:", format(obj$fn(obj$par)), "; peak", pk(), "
"); quit(save = "no") }
ev <- eigen((H + t(H)) / 2, symmetric = TRUE, only.values = TRUE)$values
cat("inner Hessian at start: dim", nrow(H), "eigen range", format(range(ev), digits = 4),
    "finite", all(is.finite(H)), "\n")
gi <- env$f(env$par, order = 1)
cat("joint gradient finite:", all(is.finite(gi)), "max |g| over random",
    format(max(abs(gi[rnd]))), "\n")
v <- obj$fn(obj$par)
lp <- env$last.par
cat("fn:", format(v), "; joint nll at last.par:", format(env$f(lp)), "\n")
H2 <- as.matrix(env$spHess(lp, random = TRUE))
ev2 <- eigen((H2 + t(H2)) / 2, symmetric = TRUE, only.values = TRUE)$values
cat("inner Hessian at last.par: eigen range", format(range(ev2), digits = 4),
    "finite", all(is.finite(H2)), "n negative", sum(ev2 < 0), "\n")
# the y_m block and the b block separately
ib <- which(names(env$par)[rnd] == "b")
im <- which(names(env$par)[rnd] == "miss")
cat("y_m block diag range", format(range(diag(H2)[im]), digits = 4),
    "; b block eigen", format(range(eigen(H2[ib, ib], only.values = TRUE)$values), digits = 4),
    "\n")
cat("y_m block minus I, max abs:", format(max(abs(H2[im, im] - diag(length(im))))), "\n")
