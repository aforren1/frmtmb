# Reviewer of lane loglik, fix round 2: when the classed warning
# frmtmb_autocor_jacobian_unchecked fires. Seed 121.
# Output: dev/review-loglik-log/fr2-warn.txt
.libPaths(c("C:/Users/adf44/source/r/wt-loglik-lib",
            "C:/Users/adf44/source/r/rellib-r13",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({
  library(frmtmb)
  library(frmtmb.spatial)
  library(frmtmb.sample)
})
cat("lib:", find.package("frmtmb"), find.package("frmtmb.spatial"), "\n")
count_w <- function(expr) {
  nu <- 0L
  other <- character()
  val <- withCallingHandlers(expr, warning = function(w) {
    if (inherits(w, "frmtmb_autocor_jacobian_unchecked")) nu <<- nu + 1L
    else other <<- c(other, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(val = val, n = nu, other = other)
}
set.seed(121)
n <- 40
d <- data.frame(x = rnorm(n))
d$y <- 1 + 0.5 * d$x + rnorm(n)
d$y2 <- d$y

# (a) a map defined at its start only, with the right slot
P <- matrix(0, n, n)
P[cbind(1:(n - 1), 2:n)] <- 1
P <- P + t(P)
Wr <- P / rowSums(P)
only_start <- function(theta) {
  if (!inherits(theta, "advector") && abs(theta[1] - 0.2) > 1e-12) {
    stop("defined at its start only")
  }
  theta[1]
}
frmtmb_register_autocor(
  "rvonly", args = "W",
  build = function(spec, rows, n_all, lookup) {
    W <- lookup(spec$W, "rvonly()")[rows, rows, drop = FALSE]
    list(W = W, ev = Re(eigen(W, only.values = TRUE)$values))
  },
  whiten = function(y, mu, sigma, theta, aux) {
    A <- diag(length(y)) - only_start(theta) * aux$W
    list(u = as.vector(A %*% (y - mu)) / sigma,
         logjac = sum(log1p(-only_start(theta) * aux$ev)) -
           length(y) * log(sigma))
  },
  unwhiten = function(u, mu, sigma, theta, aux) {
    A <- diag(nrow(aux$W)) - only_start(theta) * aux$W
    mu + sigma * solve(A, u)
  },
  npar = 1L, start = function(aux) 0.2,
  natural = function(theta, aux) c(rho = theta[1]), sigma_varies = FALSE,
  jacobian = function(sigma, theta, aux) {
    (diag(nrow(aux$W)) - only_start(theta) * aux$W) / sigma
  })
r <- count_w(frm(y ~ x + rvonly(W), data = d, data2 = list(W = Wr),
                 dry_run = "objective"))
cat("rvonly: warnings during frm(dry_run):", r$n, "\n")
fit <- r$val
lab <- c(brms_par_labels(fit), "lp__")
m <- matrix(0, 5, length(lab), dimnames = list(NULL, lab))
m[, "b_Intercept"] <- 1
m[, "sigma"] <- 1
m[, "thetaac_1"] <- 0.2
ds <- structure(list(stanfit = NULL, draws = m, fit = fit),
                class = "frmtmb_draws")
r2 <- count_w(log_lik(ds))
cat("rvonly: warnings during log_lik() of 5 draws:", r2$n,
    " finite:", all(is.finite(r2$val)), "\n")
r3 <- count_w(frm(bf(y ~ x + rvonly(W)) + bf(y2 ~ x + rvonly(W)) +
                    set_rescor(FALSE), data = d, data2 = list(W = Wr),
                  dry_run = "objective"))
cat("rvonly in two responses: warnings during frm(dry_run):", r3$n, "\n")

# (b) correct terms at their own starts never warn: fcor, sar lag and
# error over several W, fitted, and assembled for sampling
V <- 0.6^abs(outer(1:n, 1:n, "-")) + diag(0.3, n)
dg <- 1 / sqrt(rowSums(P))
Ws <- P * outer(dg, dg)
Ps <- matrix(0, n, n)
Ps[1, -1] <- 1
Ps[-1, 1] <- 1
Pstar <- Ps + P * 0
Wstar <- Pstar / rowSums(Pstar)          # star graph: eigenvalues -1, 0, 1
Wbin <- P                                 # unnormalized: rho in about (-0.5, 0.5)
k <- 8
id <- matrix(seq_len(k * k), k, k)
Pl <- matrix(0, k * k, k * k)
Pl[cbind(as.vector(id[-k, ]), as.vector(id[-1, ]))] <- 1
Pl[cbind(as.vector(id[, -k]), as.vector(id[, -1]))] <- 1
Pl <- Pl + t(Pl)
Wlat <- Pl / rowSums(Pl)
# eigenvalues in (-0.05, 1]: the lower bound 1 / min is about -20, so
# rho's interval is far from symmetric about the start rho = 0
Wskew <- 0.95 * diag(0, n) + Wr
Wskew <- (Wskew + diag(0.9, n)) / 1.9
diag(Wskew) <- 0
Wskew <- Wskew / rowSums(Wskew)
Wtiny <- Wr * 40                          # rho's interval about (-0.025, 0.025)
dl <- data.frame(x = rnorm(k * k))
dl$y <- 1 + rnorm(k * k)
cases <- list(
  list("fcor(V)", y ~ x + fcor(V), d, list(V = V)),
  list("fcor(V), sigma ~ x", bf(y ~ x + fcor(V), sigma ~ x), d, list(V = V)))
for (wn in c("Wr", "Ws", "Wstar", "Wbin", "Wskew", "Wtiny")) {
  W <- get(wn)
  for (type in c("lag", "error")) {
    cases[[length(cases) + 1L]] <- list(
      paste0("sar(", wn, ", ", type, ")"),
      stats::as.formula(paste0("y ~ x + sar(W, type = '", type, "')")),
      d, list(W = W))
  }
}
for (type in c("lag", "error")) {
  cases[[length(cases) + 1L]] <- list(
    paste0("sar(Wlat 8x8, ", type, ")"),
    stats::as.formula(paste0("y ~ x + sar(W, type = '", type, "')")),
    dl, list(W = Wlat))
}
for (cs in cases) {
  for (dr in c("objective", "fit")) {
    r <- tryCatch(count_w(if (dr == "fit") {
      frm(cs[[2]], data = cs[[3]], data2 = cs[[4]])
    } else {
      frm(cs[[2]], data = cs[[3]], data2 = cs[[4]], dry_run = "objective")
    }), error = function(e) list(n = NA, other = conditionMessage(e)))
    lims <- ""
    if (!is.null(r$val) && grepl("^sar", cs[[1]])) {
      aux <- r$val$frame$autocor$y$aux
      lims <- sprintf("rho in (%.3g, %.3g)", aux$lo, aux$hi)
    }
    cat(sprintf("%-26s %-9s unchecked warnings %s  other: %s  %s\n",
                cs[[1]], dr, r$n,
                if (length(r$other)) substr(paste(unique(r$other),
                                                  collapse = " | "), 1, 80)
                else "none", lims))
  }
}
cat("DONE\n")
