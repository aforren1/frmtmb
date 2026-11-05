# Short check of fix round 2, item 1: chol_L's `is.numeric(K)` on the tape. Is it
# TRUE on an advector here? Can the tape hold a jitter at the sampler's 1e-12: at
# which starts does the factor fail while taping? If it does, how far is that tape
# from the plain one, value and gradient, at a theta where both factor? Installed
# build only; the wrappers live in this process. Seeds 8801 up.
#
#   bash dev/review-gpgap-run.sh dev/review-gpgap-tape.R tape.txt
.libPaths(c("C:/Users/adf44/source/r/wt-gpgap-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
ns <- asNamespace("frmtmb")
nss <- asNamespace("frmtmb.sample")
tst <- RTMB::MakeADFun(function(p) {
  K <- diag(2) * p$a
  cat(sprintf("TAPE is.numeric(advector) %s, class %s\n", is.numeric(K), class(K)[1]))
  sum(K)
}, list(a = 1), silent = TRUE)
orig <- get("gp_chol_num", envir = ns)
LOG <- new.env()
LOG$rows <- list()
logged <- function(K) {
  adv <- inherits(K, "advector")
  R <- tryCatch(chol(K), error = function(e) NULL)
  jit <- 1e-12
  used <- 0
  while (is.null(R) && jit <= frmtmb:::gp_nugget) {
    Kj <- K
    diag(Kj) <- diag(Kj) + jit
    R <- tryCatch(chol(Kj), error = function(e) NULL)
    used <- jit
    jit <- jit * 10
  }
  LOG$rows[[length(LOG$rows) + 1L]] <- c(adv = adv, used = used)
  if (is.null(R)) chol(K) else R
}
plain <- function(K) chol(K)
put <- function(f) {
  unlockBinding("gp_chol_num", ns); assign("gp_chol_num", f, envir = ns)
  lockBinding("gp_chol_num", ns)
}
designs <- list(
  sin15 = function() { x <- round(runif(15, 0, 6), 3); data.frame(x, y = sin(2 * x) + 1e-3 * rnorm(15)) },
  sin60 = function() { x <- round(runif(60, 0, 6), 3); data.frame(x, y = sin(2 * x) + 0.1 * rnorm(60)) },
  sin200 = function() { x <- round(runif(200, 0, 6), 3); data.frame(x, y = sin(2 * x) + 0.1 * rnorm(200)) },
  sin480 = function() { x <- round(runif(480, 0, 6), 4); data.frame(x, y = sin(2 * x) + 0.1 * rnorm(480)) },
  trend100 = function() { x <- round(runif(100, 0, 10), 3); data.frame(x, y = x^2 + 0.3 * rnorm(100)) },
  dense200 = function() { x <- seq(0, 0.2, length.out = 200); data.frame(x, y = sin(2 * x) + 0.01 * rnorm(200)) }
)
sd0 <- 8800L
for (nm in names(designs)) {
  sd0 <- sd0 + 1L
  set.seed(sd0)
  d <- designs[[nm]]()
  f <- suppressWarnings(frm(bf(y ~ gp(x)), family = gaussian(), data = d))
  f <- nss$gp_sample_nugget(f, 1L)
  th0 <- f$estimates$theta
  for (mult in c(1, 10, 100, 1e3, 1e4, 1e6)) {
    fi <- f
    fi$estimates$theta[2] <- th0[2] + log(mult)
    put(logged)
    LOG$rows <- list()
    r <- tryCatch({
      o <- nss$ncp_objective(fi, 1L, list(), NULL)
      c(ok = TRUE)
    }, error = function(e) {
      cat("   tape error:", substr(conditionMessage(e), 1, 90), "\n")
      c(ok = FALSE)
    })
    put(orig)
    rows <- do.call(rbind, LOG$rows)
    tp <- if (is.null(rows)) NULL else rows[rows[, "adv"] == 1, , drop = FALSE]
    kn <- exp(-Reduce(`+`, f$frame$re_blocks[[1]]$aux_D2) / (2 * exp(2 * fi$estimates$theta[2])))
    diag(kn) <- diag(kn) + 1e-12
    p <- f$frame$re_blocks[[1]]$gp_order
    ch <- !inherits(try(chol(kn[p, p]), silent = TRUE), "try-error")
    cat(sprintf(paste0("TAPE %-9s rho x %-6g (rho %.3g) | R's chol at the start %s | taping: built %s, calls on advector %d, ",
                       "jitter on the tape %s\n"),
                nm, mult, exp(fi$estimates$theta[2]), ch, r[["ok"]],
                if (is.null(tp)) 0L else nrow(tp),
                if (is.null(tp) || !nrow(tp)) "none" else paste(unique(tp[, "used"]), collapse = ",")))
  }
}
# what a taped jitter does: force the failure by a smaller nugget (as ncgp did),
# and read the tape at a theta where both factor
set.seed(8810)
d <- designs$sin60()
f <- suppressWarnings(frm(bf(y ~ gp(x)), family = gaussian(), data = d))
f <- nss$gp_sample_nugget(f, 1L)
for (nu in c(1e-12, 1e-14, 1e-16, 1e-17)) {
  fn <- f
  fn$frame$re_blocks[[1]]$gp_nugget <- nu
  thg <- f$estimates$theta
  fb <- fn
  fb$estimates$theta[2] <- thg[2] + log(1e4)
  put(logged); LOG$rows <- list()
  o_bad <- tryCatch(nss$ncp_objective(fb, 1L, list(), NULL), error = function(e) NULL)
  rows <- do.call(rbind, LOG$rows); tp <- rows[rows[, "adv"] == 1, , drop = FALSE]
  LOG$rows <- list()
  o_ok <- tryCatch(nss$ncp_objective(fn, 1L, list(), NULL), error = function(e) NULL)
  rows2 <- do.call(rbind, LOG$rows); tp2 <- rows2[rows2[, "adv"] == 1, , drop = FALSE]
  put(orig)
  if (is.null(o_bad) || is.null(o_ok)) {
    cat(sprintf("JIT nu %g: a tape did not build\n", nu))
    next
  }
  # one parameter point, a typical z at the good theta
  par <- o_ok$env$last.par.best
  nm <- names(par)
  set.seed(1)
  par[nm == "b"] <- rnorm(sum(nm == "b"))
  par[nm == "theta"] <- thg
  v1 <- o_ok$env$f(par); v2 <- o_bad$env$f(par)
  g1 <- as.numeric(o_ok$env$f(par, order = 1)); g2 <- as.numeric(o_bad$env$f(par, order = 1))
  cat(sprintf(paste0("JIT nu %g | tape built at rho x 1e4: jitter %s; at the fit's rho: jitter %s | at the fit's theta: ",
                     "value %.6g against %.6g, diff %.3e; max |grad| %.3e, max |grad diff| %.3e\n"),
              nu, paste(unique(tp[, "used"]), collapse = ","), paste(unique(tp2[, "used"]), collapse = ","),
              v1, v2, v2 - v1, max(abs(g1)), max(abs(g2 - g1))))
}
