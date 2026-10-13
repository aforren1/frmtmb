# Reviewer of lane optrap, re-check: attacks on newton_polish(), on the
# lane's build, in one process. Parts:
#   saddle  a direct call at a point whose Hessian is indefinite
#   bound   degree-6 raw polynomials (dev/optrap-p1-polystudy.R's data)
#           with an upper bound on the x6 coefficient below lm()'s value
#   cost    the time newton_polish() adds to fixed-effect fits of 2e5
#           rows: gaussian with 99 covariates, poisson with 50, Gamma
#           with 20; interleaved with the step patched out, 3 rounds
#   bign    how often a step is taken on y ~ x1 + ... + x9 at n = 1e2,
#           1e4 and 1e6 (seeds 1 to 10), and its gain
#   Rscript dev/optrap-rev2-polish.R <part>
a <- commandArgs(TRUE)
.libPaths(c("C:/Users/adf44/source/r/wt-optrap-lib",
            "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
np_orig <- get("newton_polish", ns)
set_polish <- function(on) {
  utils::assignInNamespace("newton_polish",
                           if (on) np_orig else function(obj, opt, ...) opt,
                           "frmtmb")
}
part <- a[1]
cat("part", part, "lib", dirname(find.package("frmtmb")), "\n")
if (part == "saddle") {
  # f = x^2 - y^2 + 0.1 y: no minimum; the Newton step from (1, 0.1)
  # goes to the stationary point (0, 0.05), a saddle, with f lower
  obj <- RTMB::MakeADFun(function(p) p$x^2 - p$y^2 + 0.1 * p$y,
                         list(x = 1, y = 0.1), silent = TRUE)
  opt <- list(par = c(x = 1, y = 0.1), objective = obj$fn(c(1, 0.1)),
              convergence = 0L)
  out <- np_orig(obj, opt, list(lower = -Inf, upper = Inf))
  cat("start", opt$par, "f", opt$objective, "\n")
  cat("after", out$par, "f", out$objective, "steps",
      out$polish[["steps"]] %||% 0, "\n")
  cat("Hessian eigenvalues at the end:",
      eigen(obj$he(out$par))$values, "\n")
}
if (part == "bound") {
  for (s in 1:10) {
    set.seed(6000 + s)
    d <- data.frame(x = runif(200, 1, 2))
    d$y <- sin(3 * d$x) + rnorm(200, 0, 0.2)
    for (k in 1:6) d[[paste0("x", k)]] <- d$x^k
    fo <- y ~ x1 + x2 + x3 + x4 + x5 + x6
    b6 <- coef(lm(fo, data = d))[["x6"]]
    ub <- b6 - 0.5 * abs(b6)
    pr <- set_prior("", class = "b", coef = "x6", ub = ub)
    res <- list()
    for (on in c(TRUE, FALSE)) {
      set_polish(on)
      f <- suppressWarnings(frm(fo, data = d, prior = pr))
      res[[if (on) "on" else "off"]] <- f
    }
    set_polish(TRUE)
    x6 <- function(f) f$estimates$beta[7]
    cat(sprintf("seed %d ub %.4g | polish: obj %.10f x6 %.6g steps %s code %d | off: obj %.10f x6 %.6g | x6 <= ub %s\n",
                s, ub, res$on$opt$objective, x6(res$on),
                res$on$opt$polish[["steps"]] %||% 0, res$on$opt$convergence,
                res$off$opt$objective, x6(res$off), x6(res$on) <= ub))
  }
}
if (part == "cost") {
  set.seed(1)
  n <- 2e5
  mk <- function(p) {
    Z <- matrix(rnorm(n * p), n, p)
    colnames(Z) <- paste0("z", seq_len(p))
    as.data.frame(Z)
  }
  cases <- list(
    gauss99 = function() {
      d <- mk(99)
      d$y <- 1 + 0.5 * d$z1 + rnorm(n)
      list(d = d, fam = gaussian())
    },
    pois50 = function() {
      d <- mk(50)
      d$y <- rpois(n, exp(0.3 + 0.2 * d$z1))
      list(d = d, fam = poisson())
    },
    gamma20 = function() {
      d <- mk(20)
      d$y <- rgamma(n, 2, 2 / exp(0.3 + 0.2 * d$z1))
      list(d = d, fam = Gamma(link = "log"))
    })
  for (nm in names(cases)) {
    cs <- cases[[nm]]()
    fo <- as.formula(paste("y ~",
                           paste(setdiff(names(cs$d), "y"), collapse = "+")))
    t <- matrix(NA, 3, 2, dimnames = list(NULL, c("on", "off")))
    st <- NA
    for (r in 1:3) for (on in c(TRUE, FALSE)) {
      set_polish(on)
      t0 <- proc.time()[[3]]
      f <- suppressWarnings(frm(fo, family = cs$fam, data = cs$d))
      t[r, if (on) "on" else "off"] <- proc.time()[[3]] - t0
      if (on) st <- f$opt$polish[["steps"]] %||% 0
    }
    set_polish(TRUE)
    cat(sprintf("%s: fit secs min with steps %.2f, without %.2f (ratio %.2f); steps taken %s\n",
                nm, min(t[, "on"]), min(t[, "off"]),
                min(t[, "on"]) / min(t[, "off"]), st))
  }
}
if (part == "bign") {
  for (n in c(1e2, 1e4, 1e6)) {
    stepped <- 0L
    gains <- numeric(0)
    for (s in 1:10) {
      set.seed(s)
      Z <- matrix(rnorm(n * 9), n, 9)
      d <- as.data.frame(Z)
      names(d) <- paste0("x", 1:9)
      d$y <- 1 + Z %*% seq(0.1, 0.9, by = 0.1) + rnorm(n)
      f <- suppressWarnings(frm(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 +
                                  x9, data = d))
      if (!is.null(f$opt$polish)) {
        stepped <- stepped + 1L
        gains <- c(gains, f$opt$polish[["gain"]])
      }
    }
    cat(sprintf("n %g: fits with a step %d of 10 | gains %s\n", n, stepped,
                paste(signif(gains, 3), collapse = " ")))
  }
}
