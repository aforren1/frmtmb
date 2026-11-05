# Short check of fix round 2, item 1, second part. Harder designs for a jitter on the
# tape at 1e-12: tied rows under gr = FALSE, positions 1e-9 apart, 2-D, by. What a taped
# jitter does when it is there, against the same nugget without one. And the numeric
# route: how many draws of a real run does gp_chol_num() jitter? Seed 8821 up.
#
#   bash dev/review-gpgap-run.sh dev/review-gpgap-tape2.R tape2.txt
.libPaths(c("C:/Users/adf44/source/r/wt-gpgap-lib",
            "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
ns <- asNamespace("frmtmb")
nss <- asNamespace("frmtmb.sample")
orig <- get("gp_chol_num", envir = ns)
LOG <- new.env(); LOG$rows <- list()
logged <- function(K) {
  adv <- inherits(K, "advector")
  R <- tryCatch(chol(K), error = function(e) NULL)
  jit <- 1e-12; used <- 0
  while (is.null(R) && jit <= frmtmb:::gp_nugget) {
    Kj <- K; diag(Kj) <- diag(Kj) + jit
    R <- tryCatch(chol(Kj), error = function(e) NULL)
    used <- jit; jit <- jit * 10
  }
  LOG$rows[[length(LOG$rows) + 1L]] <- c(adv = adv, used = used)
  if (is.null(R)) chol(K) else R
}
put <- function(f) {
  unlockBinding("gp_chol_num", ns); assign("gp_chol_num", f, envir = ns)
  lockBinding("gp_chol_num", ns)
}
designs <- list(
  tied_grF = function() {
    x <- sample(seq(0, 6, by = 0.5), 300, TRUE)
    frm(bf(y ~ gp(x, gr = FALSE)), data = data.frame(x, y = sin(2 * x) + 0.1 * rnorm(300)))
  },
  near1e9 = function() {
    x <- rep(runif(60, 0, 6), each = 3) + rep(c(0, 1e-9, 2e-9), 60)
    frm(bf(y ~ gp(x)), data = data.frame(x, y = sin(2 * x) + 0.1 * rnorm(180)))
  },
  twoD = function() {
    d <- data.frame(x = runif(200, 0, 6), z = runif(200, 0, 6))
    d$y <- sin(d$x) * cos(d$z) + 0.05 * rnorm(200)
    frm(bf(y ~ gp(x, z, iso = FALSE)), data = d)
  },
  by = function() {
    d <- data.frame(x = runif(200, 0, 6), f = factor(rep(1:2, 100)))
    d$y <- sin(2 * d$x) + 0.05 * rnorm(200)
    frm(bf(y ~ gp(x, by = f)), data = d)
  })
s0 <- 8820L
for (nm in names(designs)) {
  s0 <- s0 + 1L
  set.seed(s0)
  f <- suppressWarnings(designs[[nm]]())
  bks <- which(vapply(f$frame$re_blocks, function(b) identical(b$covstruct, "gp"), NA))
  f <- nss$gp_sample_nugget(f, bks)
  th0 <- f$estimates$theta
  for (mult in c(1, 1e2, 1e4, 1e6)) {
    fi <- f
    for (i in bks) {
      idx <- f$frame$re_blocks[[i]]$theta_idx
      fi$estimates$theta[idx[-1]] <- th0[idx[-1]] + log(mult)
    }
    put(logged); LOG$rows <- list()
    r <- tryCatch({
      nss$ncp_objective(fi, bks, list(), NULL)
      TRUE
    }, error = function(e) substr(conditionMessage(e), 1, 70))
    put(orig)
    rows <- do.call(rbind, LOG$rows)
    tp <- if (is.null(rows)) NULL else rows[rows[, "adv"] == 1, , drop = FALSE]
    cat(sprintf("TAPE2 %-9s rho x %-6g built %s, advector calls %d, jitter on the tape: %s\n", nm, mult, r,
                if (is.null(tp)) 0L else nrow(tp),
                if (is.null(tp) || !nrow(tp)) "none" else paste(unique(tp[, "used"]), collapse = ",")))
  }
}
# a taped jitter against the same nugget without it: nugget 1e-16 tapes 1e-12 (tape.txt)
set.seed(8810)
x <- round(runif(60, 0, 6), 3)
d <- data.frame(x, y = sin(2 * x) + 1e-3 * rnorm(60))
f <- nss$gp_sample_nugget(suppressWarnings(frm(bf(y ~ gp(x)), data = d)), 1L)
fj <- f; fj$frame$re_blocks[[1]]$gp_nugget <- 1e-16
fp <- f; fp$frame$re_blocks[[1]]$gp_nugget <- 2e-12 # what the jittered tape really is, nu + 1e-12...
ft <- f
o12 <- nss$ncp_objective(ft, 1L, list(), NULL)
oj <- nss$ncp_objective(fj, 1L, list(), NULL)
op <- nss$ncp_objective(fp, 1L, list(), NULL)
par <- o12$env$last.par.best
nm <- names(par)
set.seed(1)
par[nm == "b"] <- rnorm(sum(nm == "b"))
par[nm == "theta"] <- f$estimates$theta
par[nm == "betad"] <- log(1e-3)
rd <- function(o) list(v = o$env$f(par), g = as.numeric(o$env$f(par, order = 1)))
a <- rd(o12); b <- rd(oj); c1 <- rd(op)
cat(sprintf(paste0("JITX nu 1e-16 tape (holds jitter 1e-12) against a tape at nu 1e-12: value diff %.3e, max grad diff %.3e, of max grad %.3e\n",
                   "JITX tape at nu 2e-12 against nu 1e-12 (what a second jitter step would be): value diff %.3e, max grad diff %.3e\n"),
            b$v - a$v, max(abs(b$g - a$g)), max(abs(a$g)), c1$v - a$v, max(abs(c1$g - a$g))))
# the numeric route: a real run, counted
for (arm in list(list(tag = "low n 40", s = 1e-3, A = 1, n = 40L), list(tag = "units n 60", s = 0.1, A = 1000, n = 60L))) {
  set.seed(6601)
  x <- round(runif(arm$n, 0, 6), 3)
  d <- data.frame(x, y = arm$A * (sin(2 * x) + arm$s * rnorm(arm$n)))
  put(logged); LOG$rows <- list()
  fs <- suppressMessages(suppressWarnings(frm_sample(bf(y ~ gp(x)), data = d, family = gaussian(), chains = 2,
                                                     iter = 300, warmup = 150, seed = 3, refresh = 0, cores = 1,
                                                     control = list(max_treedepth = 8))))
  invisible(brms::posterior_epred(fs, newdata = data.frame(x = c(0.5, 2.2))))
  put(orig)
  rows <- do.call(rbind, LOG$rows)
  cat(sprintf("NUMR %s: numeric factor calls %d, jittered %d (sizes %s); advector calls %d with jitter %d\n", arm$tag,
              sum(rows[, "adv"] == 0), sum(rows[, "adv"] == 0 & rows[, "used"] > 0),
              paste(unique(rows[rows[, "adv"] == 0, "used"]), collapse = ","),
              sum(rows[, "adv"] == 1), sum(rows[, "adv"] == 1 & rows[, "used"] > 0)))
}
