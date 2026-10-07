# Reviewer of lane setier, re-check: why the loss-ratio window loses the
# identified near-collinear pair of dev/setier-rev2-window.R, and three
# probe rules compared on the same directions:
#   r1   round 1: both full steps lose at least grad_tol
#   lane punch round 1: also each side's full/half loss ratio in (3, 5.5)
#   sym  symmetric second differences, which cancel the linear term a
#        fit short of exact stationarity leaves: c(t) = f(+t) + f(-t) -
#        2 f0; keep when c(full) >= 2 grad_tol and c(full) / c(half) is
#        in (3, 5.5)
# Each rule is patched into the lane's namespace in turn (this process
# only), and the designs are refitted.
#   Rscript dev/setier-rev2-probe3.R
.libPaths(c("C:/Users/adf44/source/r/wt-setier-lib",
            "C:/Users/adf44/source/r/rellib-r6",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages({library(frmtmb); library(lme4)})
ns <- asNamespace("frmtmb")
lane_rule <- get("se_curvature_real", ns)
mk <- function(kind) {
  f <- function(fit, p, free, dir_free, lambda) {
    tol <- fit$control$grad_tol %||% 1e-3
    obj <- fit$obj
    saved <- obj_state_save(obj)
    on.exit(obj_state_restore(obj, saved), add = TRUE)
    d <- numeric(length(p))
    d[free] <- dir_free * sqrt(4 * tol / lambda)
    ev <- function(x) tryCatch(obj$fn(x), error = function(e) NA_real_)
    f0 <- ev(p)
    fu <- ev(p + d); fd <- ev(p - d)
    hu <- ev(p + d / 2); hd <- ev(p - d / 2)
    if (!all(is.finite(c(f0, fu, fd, hu, hd)))) return(FALSE)
    if (kind == "r1") return(fu - f0 >= tol && fd - f0 >= tol)
    cf <- fu + fd - 2 * f0
    ch <- hu + hd - 2 * f0
    cf >= 2 * tol && ch > 0 && cf / ch > 3 && cf / ch < 5.5
  }
  en <- new.env(parent = ns)
  en$kind <- kind
  environment(f) <- en
  f
}
rules <- list(r1 = mk("r1"), lane = lane_rule, sym = mk("sym"))
unlockBinding("se_curvature_real", ns)
fits <- list(
  pois = function(s) {
    set.seed(s); n <- 400
    d <- data.frame(x1 = rnorm(n)); d$x2 <- d$x1 + rnorm(n) * 1e-5
    d$y <- rpois(n, exp(-1 + 0.3 * d$x1))
    list(f = frm(y ~ x1 + x2, family = poisson(), data = d),
         ref = sqrt(diag(vcov(glm(y ~ x1 + x2, poisson, data = d)))))
  },
  bern = function(s) {
    set.seed(s); n <- 400
    d <- data.frame(x1 = rnorm(n)); d$x2 <- d$x1 + rnorm(n) * 1e-5
    d$y <- rbinom(n, 1, plogis(-3 + 0.5 * d$x1))
    list(f = frm(y ~ x1 + x2, family = bernoulli(), data = d),
         ref = sqrt(diag(vcov(glm(y ~ x1 + x2, binomial, data = d)))))
  },
  ppoly = function(s) {
    set.seed(s); n <- 400
    d <- data.frame(x = runif(n, 1, 2))
    for (k in 1:5) d[[paste0("x", k)]] <- d$x^k
    d$y <- rpois(n, exp(0.5 + sin(3 * d$x)))
    fo <- y ~ x1 + x2 + x3 + x4 + x5
    list(f = frm(fo, family = poisson(), data = d),
         ref = sqrt(diag(vcov(glm(fo, poisson, data = d)))))
  },
  gpoly_re = function(s) {
    set.seed(s)
    d <- data.frame(x = runif(200, 1, 2), g = factor(rep(1:20, 10)))
    d$y <- sin(3 * d$x) + rnorm(20, 0, 0.3)[d$g] + rnorm(200, 0, 0.2)
    for (k in 1:6) d[[paste0("x", k)]] <- d$x^k
    fo <- y ~ x1 + x2 + x3 + x4 + x5 + x6 + (1 | g)
    list(f = frm(fo, data = d),
         ref = sqrt(diag(as.matrix(vcov(suppressMessages(suppressWarnings(
           lmer(fo, data = d, REML = FALSE))))))))
  },
  c0k = function(s) {
    set.seed(77)
    G <- 8
    dn <- data.frame(g = factor(rep(seq_len(G), each = 10)))
    dn$x <- rnorm(nrow(dn))
    u <- rnorm(G, 0, 0.6)[dn$g]
    dn$yn <- 1 + 0.3 * dn$x + exp(0.3 + u)^0.8 + rnorm(nrow(dn), 0, 0.3)
    list(f = frm(bf(yn ~ c0 + exp(a)^k, c0 ~ 1 + x, a ~ 1 + (1 | g),
                    k ~ 1, nl = TRUE), data = dn), ref = NULL)
  })
seeds <- list(pois = 1:10, bern = 1:10, ppoly = 1:10, gpoly_re = 1:5,
              c0k = 77)
for (rn in names(rules)) {
  assign("se_curvature_real", rules[[rn]], envir = ns)
  cat("\n== rule", rn, "\n")
  for (des in names(fits)) {
    nl <- 0; rels <- numeric()
    for (s in seeds[[des]]) {
      r <- suppressMessages(suppressWarnings(fits[[des]](s)))
      lost <- ns$sdr_of(r$f)$se_lost
      if (length(lost)) nl <- nl + 1 else if (!is.null(r$ref)) {
        se <- suppressWarnings(fixef(r$f))[, "Est.Error"]
        rels <- c(rels, max(abs(unname(se) / unname(r$ref) - 1)))
      }
      if (des == "c0k") cat("  c0k lost:", paste(names(lost), collapse = ","),
                            "\n")
    }
    if (des != "c0k") {
      cat(sprintf("  %-9s fits %d | losing an SE %d | max |SE/ref-1| of the rest %s\n",
                  des, length(seeds[[des]]), nl,
                  if (length(rels)) signif(max(rels), 3) else "-"))
    }
  }
}
# the mechanism on pois seed 1: the linear term along the probed direction
assign("se_curvature_real", lane_rule, envir = ns)
r <- suppressWarnings(fits$pois(1))
f <- r$f
H <- f$obj$he(f$opt$par)
D <- sqrt(abs(diag(H)))
e <- eigen(H / outer(D, D), symmetric = TRUE)
k <- which.min(e$values)
d <- e$vectors[, k] / D * sqrt(4e-3 / e$values[k])
p <- f$opt$par
g <- f$obj$gr(p)
f0 <- f$obj$fn(p)
cat("\npois seed 1: eigenvalue", signif(e$values[k], 3), "| gradient . step",
    signif(sum(g * d), 3), "\n  losses full +", signif(f$obj$fn(p + d) - f0, 3),
    "-", signif(f$obj$fn(p - d) - f0, 3), "| half +",
    signif(f$obj$fn(p + d / 2) - f0, 3), "-",
    signif(f$obj$fn(p - d / 2) - f0, 3), "\n")
