# Short check 2, item 2. The box. (a) The two L1 cells with no number: unstr() on
# Beta() and Gamma(log), seed 167, the test's data, with the box and without it,
# and the exact value of each fit's point by importance sampling on each group's
# 5 latents. (b) Non-latent fits carry no box. (c) A predicted dispersion. (d)
# After nlfix's density: is the box still needed on beta_binomial()?
source("C:/Users/adf44/source/r/frmtmb-wt-latres/dev/review-latres-common.R")
options(warn = -1)
ns <- asNamespace("frmtmb")
assignInNamespace("latres_refuse_margin", function(fit) invisible(NULL), ns = "frmtmb")
setbox <- function(v) assignInNamespace("latres_disp_box", v, ns = "frmtmb")
setcap <- function(v) assignInNamespace("latres_disp_cap", v, ns = "frmtmb")
box0 <- get("latres_disp_box", ns); cap0 <- get("latres_disp_cap", ns)
capt <- function(expr) {
  w <- character(0)
  v <- withCallingHandlers(tryCatch(expr, error = function(e) e),
    warning = function(c) { w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning") },
    message = function(c) { w <<- c(w, conditionMessage(c)); invokeRestart("muffleMessage") })
  list(v = v, w = w)
}
set.seed(167)
d <- expand.grid(t = 1:5, g = factor(1:30))
d <- d[order(d$g, d$t), ]
n <- nrow(d)
d$x <- rnorm(n)
e <- as.vector(replicate(30, arima.sim(list(ar = 0.5), 5, sd = 0.5)))
eta <- 0.8 + 0.4 * d$x + e
d$count <- rpois(n, exp(eta))
d$bin <- rbinom(n, 1, plogis(eta - 0.8))
d$trials <- rbinom(n, 10, plogis(eta - 0.8))
d$tr <- 10L
d$pos <- rgamma(n, 4, 4 / exp(eta))
d$unit <- rbeta(n, 8 * plogis(eta - 0.8), 8 * (1 - plogis(eta - 0.8)))

# exact log Z at the fit's outer point, by importance sampling group by group on
# u, proposal a t5 at the fit's inner mode with the inner Hessian's scale
is_value <- function(f, fam_d, nd = 40000) {
  bk <- fi("frame_latres")(f$frame)[[1]]
  a <- bk$aux
  th <- f$estimates$theta[bk$theta_idx]
  s <- exp(th[1])
  L <- as.matrix(fi("us_chol_L")(th[-1], a$d))
  b <- f$estimates$b[bk$b_idx]
  env <- f$obj$env
  par <- f$opt$par
  f$obj$fn(par)
  H <- as.matrix(env$spHess(env$last.par, random = TRUE))
  bm <- as.numeric(env$last.par[env$random])
  fixed <- as.numeric(model.matrix(f) %*% f$estimates$beta)
  y <- f$frame$y[[1]]
  G <- a$G; dd <- a$d
  tot <- 0; se2 <- 0
  set.seed(5)
  for (g in seq_len(G)) {
    idx <- (g - 1) * dd + seq_len(dd)
    Hg <- H[idx, idx, drop = FALSE]
    R <- chol(solve(Hg))
    rows <- which(a$g_of == g)
    Lr <- L[a$J[rows], , drop = FALSE]
    z <- matrix(rnorm(nd * dd), nd, dd)
    w <- sqrt(5 / rchisq(nd, 5))
    u <- sweep(z * w, 2, 0, "+") %*% R
    u <- sweep(u, 2, bm[idx], "+")
    # log t5 density
    qf <- rowSums((z * w)^2)
    lq <- lgamma((5 + dd) / 2) - lgamma(5 / 2) - dd / 2 * log(5 * pi) - sum(log(diag(R))) - (5 + dd) / 2 * log1p(qf / 5)
    lp <- -0.5 * rowSums(u^2) - dd / 2 * log(2 * pi)
    errs <- s * u %*% t(Lr)
    et <- sweep(errs, 2, fixed[rows], "+")
    for (k in seq_along(rows)) lp <- lp + fam_d(y[rows[k]], et[, k])
    lw <- lp - lq
    mx <- max(lw)
    wv <- exp(lw - mx)
    mw <- mean(wv)
    tot <- tot + mx + log(mw)
    se2 <- se2 + var(wv) / nd / mw^2
  }
  c(value = tot, se = sqrt(se2))
}
cells <- list(
  list("Beta unstr", unit ~ x + unstr(t, g), Beta(), function(f) { phi <- exp(f$opt$par[names(f$opt$par) == "betad"][1]); function(y, e) { m <- plogis(e); dbeta(y, m * phi, (1 - m) * phi, log = TRUE) } }),
  list("Gamma() unstr", pos ~ x + unstr(t, g), Gamma(), function(f) { sh <- exp(f$opt$par[names(f$opt$par) == "betad"][1]); function(y, e) { v <- rep(-Inf, length(e)); k <- e > 0; v[k] <- dgamma(y, sh, sh * e[k], log = TRUE); v } }))
for (cc in cells) {
  for (arm in c("box", "nobox")) {
    setbox(if (arm == "box") box0 else Inf)
    r <- capt(frm(cc[[2]], family = cc[[3]], data = d))
    f <- r$v
    bk <- fi("frame_latres")(f$frame)[[1]]
    th <- f$estimates$theta[bk$theta_idx]
    C <- as.matrix(fi("us_chol_cor")(th[-1], bk$aux$d))
    ch <- !inherits(try(chol(C), silent = TRUE), "try-error")
    ld <- f$opt$par[names(f$opt$par) == "betad"][1]
    v <- is_value(f, cc[[4]](f))
    cat(sprintf("%-16s %-6s code %d ll %.4f log disp %.3f sderr %.4f max|theta cor| %.3g min eig %.3g chol %s max|cor| %.6f exact %.4f (se %.4f) ll-exact %.4f\n",
                cc[[1]], arm, f$opt$convergence, as.numeric(logLik(f)), ld, exp(th[1]), max(abs(th[-1])), min(eigen(C, TRUE, TRUE)$values), ch, max(abs(C[lower.tri(C)])), v[1], v[2], as.numeric(logLik(f)) - v[1]))
    cat("   ", substr(paste(r$w, collapse = " | "), 1, 200), "\n")
  }
  setbox(box0)
}

cat("\n== (b) non-latent fits: the upper bound\n")
fp <- frm(unit ~ x + (1 | g), family = Beta(), data = d)
cat("upper on betad, no latent:", fp$upper[names(fp$upper) == "betad"], " length frame_latres", length(fi("frame_latres")(fp$frame)), "\n")
print(body(fi("latres_disp_dpars"))[[3]])

if (nzchar(Sys.getenv("REV_PRED"))) {
cat("\n== (c) a predicted dispersion, seeds 1 to 30\n")
pred <- function(lab, sim, form, fam) {
  rows <- list()
  for (s in 1:30) {
    set.seed(s)
    dd <- expand.grid(t = 1:6, g = factor(1:25))
    dd <- dd[order(dd$g, dd$t), ]
    m <- nrow(dd)
    dd$x <- rnorm(m); dd$z <- rnorm(m)
    dd$e <- as.vector(replicate(25, arima.sim(list(ar = .5), 6, sd = .5)))
    dd$tr <- 10L
    dd$y <- sim(dd)
    for (arm in c("cap", "nocap")) {
      setcap(if (arm == "cap") cap0 else Inf)
      r <- capt(frm(form, data = dd, family = fam))
      if (inherits(r$v, "error")) { rows[[length(rows) + 1]] <- data.frame(s = s, arm = arm, ll = NA, code = NA, mx = NA, nan = NA); next }
      f <- r$v
      dp <- fi("eval_dpars")(f)[[1]]
      nmd <- setdiff(names(dp), "mu")[1]
      rows[[length(rows) + 1]] <- data.frame(s = s, arm = arm, ll = as.numeric(logLik(f)), code = f$opt$convergence, mx = max(log(dp[[nmd]])), nan = any(!is.finite(sqrt(diag(as.matrix(vcov(f, full = TRUE)))))))
    }
    setcap(cap0)
  }
  o <- do.call(rbind, rows)
  for (arm in c("cap", "nocap")) {
    q <- o[o$arm == arm, ]
    cat(sprintf("%-26s %-6s fits %d code!=0 %d largest ll %.4g median ll %.2f largest log disp %.1f NaN SE %d\n", lab, arm, sum(!is.na(q$ll)), sum(q$code != 0, na.rm = TRUE), max(q$ll, na.rm = TRUE), median(q$ll, na.rm = TRUE), max(q$mx, na.rm = TRUE), sum(q$nan, na.rm = TRUE)))
  }
  w <- reshape(o[, c("s", "arm", "ll")], idvar = "s", timevar = "arm", direction = "wide")
  print(round(summary(w$ll.cap - w$ll.nocap), 4))
}
pred("Beta, phi 1000, phi ~ z", function(d) { m <- plogis(.3 * d$x + d$e); pmin(pmax(rbeta(nrow(d), 1000 * m, 1000 * (1 - m)), 1e-4), 1 - 1e-4) }, frmtmb::bf(y ~ x + ar(t, g, cov = TRUE), phi ~ z), Beta())
pred("Beta, phi 8, phi ~ z", function(d) { m <- plogis(.3 * d$x + d$e); pmin(pmax(rbeta(nrow(d), 8 * m, 8 * (1 - m)), 1e-4), 1 - 1e-4) }, frmtmb::bf(y ~ x + ar(t, g, cov = TRUE), phi ~ z), Beta())
pred("bb, binomial data, phi ~ z", function(d) rbinom(nrow(d), 10, plogis(.3 * d$x + d$e)), frmtmb::bf(y | trials(tr) ~ x + ar(t, g, cov = TRUE), phi ~ z), beta_binomial())

}
cat("\n== (d) nlfix's density on beta_binomial(), the box off\n")
src <- readLines("C:/Users/adf44/source/r/frmtmb-wt-nlfix/R/families.R")
i1 <- grep("^lgamma_count_shift <- ", src); i2 <- grep("^beta_binomial_lpmf <- ", src)
txt <- c(src[i1:(i1 + which(src[(i1 + 1):length(src)] == "}")[1])], src[i2:(i2 + which(src[(i2 + 1):length(src)] == "}")[1])])
enl <- new.env(parent = ns)
eval(parse(text = txt), enl)
fx <- get("beta_binomial_lpmf", enl)
real_db <- RTMBdist::dbetabinom
patched <- function(x, size, shape1, shape2, log = FALSE) {
  if (inherits(x, c("osa", "simref"))) return(real_db(x, size, shape1, shape2, log = log))
  fx(x, size, shape1, shape2, shape1 + shape2)
}
environment(patched) <- enl
assignInNamespace("dbetabinom", patched, ns = "RTMBdist")
for (arm in c("nlfix, box 15", "nlfix, no box")) {
  setbox(if (grepl("no box", arm)) Inf else box0)
  rows <- list()
  for (s in 1:30) {
    set.seed(s)
    dd <- expand.grid(t = 1:6, g = factor(1:25))
    dd <- dd[order(dd$g, dd$t), ]
    m <- nrow(dd)
    dd$x <- rnorm(m)
    dd$e <- as.vector(replicate(25, arima.sim(list(ar = .5), 6, sd = .5)))
    dd$tr <- 10L
    dd$y <- rbinom(m, 10, plogis(.3 * dd$x + dd$e))
    r <- capt(frm(y | trials(tr) ~ x + ar(t, g, cov = TRUE), data = dd, family = beta_binomial()))
    f <- r$v
    if (inherits(f, "error")) next
    rows[[length(rows) + 1]] <- data.frame(s = s, ll = as.numeric(logLik(f)), code = f$opt$convergence, ld = f$opt$par[names(f$opt$par) == "betad"][1])
  }
  o <- do.call(rbind, rows)
  cat(sprintf("%-16s fits %d code!=0 %d largest ll %.4g median ll %.3f largest log phi %.1f\n", arm, nrow(o), sum(o$code != 0), max(o$ll), median(o$ll), max(o$ld)))
  assign(paste0("o_", sub(".*, ", "", sub(" ", "_", arm))), o)
}
setbox(box0)
