# Reviewer: is quadrature = TRUE over two scalar blocks the exact
# iterated integral? An independent Gauss-Hermite integral written from
# the model alone (own density, own nodes), its own optimum, and its own
# numerical Hessian, against frmtmb's quadrature fit.
#   Rscript dev/fuzzfix-rev-quad-exact.R <lib|base> <model>
# model: bern | pois | cum | slope | dbar | musig
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r10"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (a[1] == "base") c(base, ul) else c(a[1], base, ul))
suppressMessages(library(frmtmb))
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
model <- a[2]
`%||%` <- function(a, b) if (is.null(a)) b else a

gh <- function(k) {
  i <- seq_len(k - 1)
  J <- matrix(0, k, k)
  J[cbind(i, i + 1)] <- J[cbind(i + 1, i)] <- sqrt(i / 2)
  e <- eigen(J, symmetric = TRUE)
  # nodes and log weights for E[f(Z)], Z ~ N(0, 1)
  list(z = sqrt(2) * e$values, lw = log(e$vectors[1, ]^2))
}
lse <- function(v) { m <- max(v); m + log(sum(exp(v - m))) }
lse_rows <- function(M) { m <- apply(M, 1, max); m + log(rowSums(exp(M - m))) }

seed <- c(bern = 101, pois = 102, cum = 103, slope = 104, dbar = 105,
          musig = 106)[[model]]
set.seed(seed)
ng <- 14; nb <- 2; m <- 4
d <- expand.grid(r = seq_len(m), gb = factor(c("i", "ii")),
                 ga = factor(seq_len(ng)))
n <- nrow(d)
d$x <- rnorm(n)
d$w <- runif(n, 0.5, 2)
d$cell <- interaction(d$ga, d$gb, drop = TRUE)
ua <- rnorm(ng, 0, 0.8)[d$ga]
vab <- rnorm(nlevels(d$cell), 0, 0.7)[d$cell]
eta <- 0.2 + 0.5 * d$x + ua
if (model == "slope") {
  eta <- eta + d$x * vab * 1.2
} else if (model == "dbar") {
  eta <- eta + d$x * rnorm(ng, 0, 0.9)[d$ga]
} else if (model != "musig") {
  eta <- eta + vab
}

# outer group (u) and inner cell (v) per model, and the multiplier of v
outer <- d$ga
inner <- if (model %in% c("dbar", "musig")) d$ga else d$cell
vmult <- if (model %in% c("slope", "dbar")) d$x else rep(1, n)
wts <- if (model %in% c("pois", "cum", "slope")) d$w else rep(1, n)

if (model == "bern" || model == "slope" || model == "dbar") {
  d$y <- rbinom(n, 1, plogis(eta))
  fam <- bernoulli()
  ldens <- function(eta, p) dbinom(d$y, 1, plogis(eta), log = TRUE)
} else if (model == "pois") {
  d$y <- rpois(n, exp(eta))
  fam <- poisson()
  ldens <- function(eta, p) dpois(d$y, exp(eta), log = TRUE)
} else if (model == "cum") {
  d$y <- factor(cut(eta + rlogis(n), c(-Inf, -0.5, 0.7, Inf),
                    labels = FALSE), ordered = TRUE)
  fam <- cumulative()
  yi <- as.integer(d$y)
  ldens <- function(eta, p) {
    th <- c(-Inf, p$thr, Inf)
    log(plogis(th[yi + 1] - eta) - plogis(th[yi] - eta))
  }
} else if (model == "musig") {
  # (1 | g) in mu and (1 | g) in sigma: two blocks over one factor
  ls <- -0.3 + rnorm(ng, 0, 0.5)[d$ga]
  d$y <- eta + rnorm(n, 0, exp(ls))
  fam <- gaussian()
}

re <- switch(model, bern = "(1 | ga/gb)", pois = "(1 | ga/gb)",
             cum = "(1 | ga/gb)", slope = "(1 | ga) + (0 + x | ga:gb)",
             dbar = "(x || ga)", musig = "(1 | ga)")
lhs <- if (all(wts == 1)) "y" else "y | weights(w)"
fo <- if (model == "musig") {
  bf(as.formula(paste(lhs, "~ x +", re)), sigma ~ 1 + (1 | ga))
} else bf(as.formula(paste(lhs, "~ x +", re)))
fo0 <- bf(as.formula(paste(lhs, "~ x")))

# parameter unpacking for the hand-written model
unpack <- function(par) {
  if (model == "cum") {
    list(thr = par[1:2], beta = c(0, par[3]), lsd = par[4:5])
  } else if (model == "musig") {
    list(beta = par[1:2], s0 = par[3], lsd = par[4:5])
  } else list(beta = par[1:2], lsd = par[3:4])
}
row_ll <- function(p, u, v) {
  e0 <- p$beta[1] + p$beta[2] * d$x
  if (model == "musig") {
    return(dnorm(d$y, e0 + u, exp(p$s0 + v), log = TRUE))
  }
  ldens(e0 + u + vmult * v, p)
}
# log marginal: prod over outer of int N(u) prod over inner-in-outer of
# int N(v) prod rows f^w, one grid of nodes per level
exact_ll <- function(par, k = 60) {
  p <- unpack(par)
  q <- gh(k)
  su <- exp(p$lsd[1]); sv <- exp(p$lsd[2])
  e0 <- p$beta[1] + p$beta[2] * d$x
  tot <- 0
  for (g in levels(outer)) {
    rg <- which(outer == g)
    acc <- q$lw  # over u nodes
    for (c in unique(as.character(inner[rg]))) {
      rc <- rg[as.character(inner[rg]) == c]
      # K_u x K_v matrix of sum over rows of w * log f
      M <- matrix(0, k, k)
      for (i in rc) {
        if (model == "musig") {
          L <- outer(su * q$z, sv * q$z, function(u, v)
            dnorm(d$y[i], e0[i] + u, exp(p$s0 + v), log = TRUE))
        } else {
          et <- outer(su * q$z, vmult[i] * sv * q$z, `+`) + e0[i]
          L <- matrix(ldens_one(i, et, p), k, k)
        }
        M <- M + wts[i] * L
      }
      acc <- acc + lse_rows(sweep(M, 2, q$lw, `+`))
    }
    tot <- tot + lse(acc)
  }
  tot
}
ldens_one <- function(i, et, p) {
  if (model %in% c("bern", "slope", "dbar")) {
    dbinom(d$y[i], 1, plogis(et), log = TRUE)
  } else if (model == "pois") {
    dpois(d$y[i], exp(et), log = TRUE)
  } else {
    th <- c(-Inf, p$thr, Inf); yy <- as.integer(d$y[i])
    log(plogis(th[yy + 1] - et) - plogis(th[yy] - et))
  }
}

# 1. my density and parameterization against frmtmb's, no random effects
f0 <- frm(fo0, family = fam, data = d)
b0 <- f0$estimates$beta
cat("beta names:", names(b0), "\n")
p0 <- if (model == "cum") {
  list(thr = unname(b0[1:2]), beta = c(0, b0[[3]]))
} else list(beta = unname(b0[1:2]), s0 = unname(f0$estimates$betad[1]))
my0 <- if (model == "musig") {
  sum(dnorm(d$y, p0$beta[1] + p0$beta[2] * d$x, exp(p0$s0), log = TRUE))
} else sum(wts * ldens(p0$beta[1] + p0$beta[2] * d$x, p0))
cat(sprintf("DENSITY CHECK glm logLik frmtmb %.10f mine %.10f diff %.3e\n",
            as.numeric(logLik(f0)), my0, as.numeric(logLik(f0)) - my0))

# 2. the quadrature fit
t0 <- proc.time()[[3]]
fq <- tryCatch(suppressMessages(frm(fo, family = fam, data = d,
                                    quadrature = TRUE)),
               error = function(e) e)
tq <- proc.time()[[3]] - t0
if (inherits(fq, "error")) {
  cat("QUAD ERROR:", conditionMessage(fq), "\n"); quit(status = 0)
}
fl <- suppressMessages(frm(fo, family = fam, data = d))
cat(sprintf("QUAD time %.1fs left_random %d nb %d\n", tq,
            length(fq$obj$env$random), length(fq$estimates$b)))
sds <- vapply(fq$frame$re_blocks, function(bk) {
  th <- fq$estimates$theta[bk$theta_idx]
  sqrt(frmtmb:::covstruct_registry[[bk$covstruct]]$vcov(th, bk)[1, 1])
}, 1)
lab <- vapply(fq$frame$re_blocks, function(bk)
  paste(bk$term_label, bk$group_name, bk$dpar %||% ""), "")
cat("blocks:", paste(lab, sprintf("%.6f", sds), collapse = " | "), "\n")
`%||%` <- function(a, b) if (is.null(a)) b else a
# block order: outer first (u), then inner (v)
iu <- 1L; iv <- 2L
if (model == "slope") { iu <- 1L; iv <- 2L }
bq <- fq$estimates$beta
par_q <- if (model == "cum") {
  c(bq[1:2], bq[3], log(sds[c(iu, iv)]))
} else if (model == "musig") {
  c(bq[1:2], fq$estimates$betad[1], log(sds[c(iu, iv)]))
} else c(bq[1:2], log(sds[c(iu, iv)]))
par_q <- unname(par_q)
llq <- as.numeric(logLik(fq))
ex60 <- exact_ll(par_q, 60); ex100 <- exact_ll(par_q, 100)
lap_at <- -fl$obj$fn(fq$opt$par)
cat(sprintf(paste0("AT QUAD OPT logLik_quad %.10f exact60 %.10f ",
                   "exact100 %.10f | quad-exact %.3e  gh60-gh100 %.3e ",
                   "laplace_at-exact %.3e\n"),
            llq, ex60, ex100, llq - ex100, ex60 - ex100, lap_at - ex100))

# 3. my own optimum of the exact integral, from the quadrature optimum
op <- optim(par_q, function(p) -exact_ll(p, 40), method = "BFGS",
            control = list(reltol = 1e-12, maxit = 500))
ex_opt <- exact_ll(op$par, 100)
cat(sprintf("EXACT OPT logLik %.10f (quad %.10f, diff %.3e) conv %d\n",
            ex_opt, llq, llq - ex_opt, op$convergence))
cat("par quad :", sprintf("%.6f", par_q), "\n")
cat("par exact:", sprintf("%.6f", op$par), "\n")

# 4. SEs of the fixed effects: numerical Hessian of the exact integral
H <- numDeriv::hessian(function(p) -exact_ll(p, 60), op$par)
se_ex <- sqrt(diag(solve(H)))
V <- vcov(fq)
cat("vcov names:", rownames(V), "\n")
nb <- if (model == "cum") 3 else 2
se_q <- sqrt(diag(V))[seq_len(nb)]
cat("SE quad :", sprintf("%.6f", se_q), "\n")
cat("SE exact:", sprintf("%.6f", se_ex[seq_len(nb)]), "\n")
cat(sprintf("SE ratio quad/exact: %s\n",
            paste(sprintf("%.5f", se_q / se_ex[seq_len(nb)]),
                  collapse = " ")))
if (model == "musig") {
  cat(sprintf("SE sigma intercept quad %.6f exact %.6f\n",
              sqrt(diag(V))[3], se_ex[3]))
}

# 5. the post-fit methods
ws <- character()
res <- withCallingHandlers({
  sm <- tryCatch(summary(fq), error = function(e) e)
  vc <- tryCatch(VarCorr(fq), error = function(e) e)
  rf <- tryCatch(ranef(fq), error = function(e) e)
  list(sm, vc, rf)
}, warning = function(w) { ws <<- c(ws, conditionMessage(w))
  invokeRestart("muffleWarning") },
   message = function(m) invokeRestart("muffleMessage"))
cat("summary:", if (inherits(res[[1]], "error")) conditionMessage(res[[1]])
    else "ok", "| VarCorr:", if (inherits(res[[2]], "error"))
    conditionMessage(res[[2]]) else "ok", "| ranef:",
    if (inherits(res[[3]], "error")) conditionMessage(res[[3]]) else
      paste("ok, finite", all(is.finite(unlist(lapply(res[[3]], unlist))))),
    "\n")
if (length(ws)) cat("WARNINGS:", unique(ws), sep = "\n  ")
if (!inherits(res[[1]], "error")) {
  print(res[[1]]$fixed %||% res[[1]])
}
if (!inherits(res[[2]], "error")) print(res[[2]])
cat("DONE\n")
