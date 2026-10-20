# Review of lane illse2, item 1: REML fixed-effect SEs on ordinary
# (well-conditioned) mixed models, frmtmb against lme4 and the lane
# against 0.71.0.
#   ILLSE2_ARM=lane|base Rscript dev/illse2-rev-reml.R
source("C:/Users/adf44/source/r/frmtmb-wt-illse2/dev/illse2-common.R")
out_rds <- sprintf("C:/Users/adf44/source/r/frmtmb-wt-illse2/dev/illse2-rev-log/reml-%s.rds", arm)
data("sleepstudy", package = "lme4")
ss <- sleepstudy
set.seed(1)
ssu <- ss[sort(sample(nrow(ss), 120)), ]
set.seed(2)
G <- 30
ng <- sample(2:15, G, replace = TRUE)
g <- factor(rep(seq_len(G), ng))
sim <- data.frame(g = g, x = rnorm(length(g)),
                  z = rnorm(G)[g])
sim$y <- 1 + 0.5 * sim$x + 0.8 * sim$z + rnorm(G, 0, 0.7)[g] +
  rnorm(nrow(sim), 0, 1)
ssb <- ssu
ssb$t <- ssb$Days + 1e5
cases <- list(
  sleep_bal = list(fo = Reaction ~ Days + (Days | Subject), d = ss),
  sleep_unb = list(fo = Reaction ~ Days + (1 | Subject), d = ssu),
  sim_unb   = list(fo = y ~ x + z + (1 | g), d = sim),
  sleep_t1e5 = list(fo = Reaction ~ t + (1 | Subject), d = ssb)
)
# GLS of a (1 | g) model at t2 = sd^2 / sigma^2, double precision on a
# well-conditioned design
gls1 <- function(X, y, gg, t2) {
  Xt <- X; yt <- y
  for (ix in split(seq_along(y), gg)) {
    k <- length(ix)
    a <- (1 - 1 / sqrt(1 + k * t2)) / k
    yt[ix] <- y[ix] - a * sum(y[ix])
    Xt[ix, ] <- X[ix, , drop = FALSE] -
      a * rep(colSums(X[ix, , drop = FALSE]), each = k)
  }
  q <- qr(Xt)
  Ri <- backsolve(qr.R(q), diag(ncol(X)))
  Ri[q$pivot, ] <- Ri
  list(beta = qr.coef(q, yt), XtVX_inv = tcrossprod(Ri))
}
res <- list()
for (nm in names(cases)) {
  cs <- cases[[nm]]
  r <- cap(frm(cs$fo, data = cs$d, REML = TRUE))
  f <- r$v
  fe <- fixef(f)
  V <- suppressWarnings(vcov(f))
  Vf <- suppressWarnings(vcov(f, full = TRUE))
  l <- cap(lme4::lmer(cs$fo, data = cs$d, REML = TRUE))
  lv <- as.matrix(vcov(l$v))
  se_f <- unname(fe[, "Est.Error"])
  se_l <- unname(sqrt(diag(lv)))
  cat(sprintf("== %s (%s) code %d | warn %s | whiten %s | joint_root %s\n",
              nm, arm, f$opt$convergence, wtag(r$w),
              !is.null(f$beta_whiten), !is.null(ns$sdr_of(f)$joint_root)))
  cat(sprintf("  REML crit frm %.10f lme4 %.10f | lme4 warn %s\n",
              as.numeric(logLik(f)), as.numeric(logLik(l$v)), wtag(l$w)))
  cat("  est  frm:", format(unname(fe[, "Estimate"]), digits = 12), "\n")
  cat("  est lme4:", format(unname(lme4::fixef(l$v)), digits = 12), "\n")
  cat("  SE   frm:", format(se_f, digits = 12), "\n")
  cat("  SE  lme4:", format(se_l, digits = 12), "\n")
  cat("  SE frm/lme4 - 1:", format(se_f / se_l - 1, digits = 3), "\n")
  cat("  sqrt(diag vcov) / fixef SE - 1:",
      format(sqrt(diag(V)) / se_f - 1, digits = 3), "\n")
  cond <- NA
  if (nm != "sleep_bal") {
    X <- stats::model.matrix(lme4::nobars(cs$fo), cs$d)
    y <- cs$d[[all.vars(cs$fo)[1]]]
    gg <- cs$d[[if (nm == "sim_unb") "g" else "Subject"]]
    p <- f$opt$par
    lsig <- p[names(p) == "betad"]
    lsd <- p[names(p) == "theta"]
    t2 <- exp(2 * (lsd - lsig))
    gl <- gls1(X, y, gg, t2)
    Vc <- exp(2 * lsig) * gl$XtVX_inv
    h <- 1e-4
    db <- (gls1(X, y, gg, exp(log(t2) + 2 * h))$beta -
             gls1(X, y, gg, exp(log(t2) - 2 * h))$beta) / (2 * h)
    J <- cbind(-db, db)
    k <- nrow(Vf)
    Vt <- Vf[k - 1:0, k - 1:0]
    Vp <- J %*% Vt %*% t(J)
    cond <- sqrt(diag(Vc))
    cat("  SE frm / conditional (lme4 def at frm's point) - 1:",
        format(se_f / cond - 1, digits = 3), "\n")
    cat("  SE frm / sqrt(cond + propagated) - 1:",
        format(se_f / sqrt(diag(Vc + Vp)) - 1, digits = 3), "\n")
    cat("  propagated / conditional variance:",
        format(diag(Vp) / diag(Vc), digits = 3), "\n")
  }
  res[[nm]] <- list(est = unname(fe[, "Estimate"]), se = se_f, V = V,
                    Vf = Vf, se_l = se_l, cond = cond,
                    code = f$opt$convergence, w = r$w)
}
dir.create(dirname(out_rds), showWarnings = FALSE)
saveRDS(res, out_rds)
