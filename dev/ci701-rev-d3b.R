# Review of lane ci701: se_edge_cor() on us blocks of dimension 3 and 4,
# where a correlation parameter is not t / sqrt(1 + t^2) and pulling one
# parameter to 70.7 need not put any correlation at 0.9999. Each fit is
# made on the lane build twice in one process: as installed, and with
# se_edge_cor() stubbed to FALSE (0.70.0's check plus D1).
#
#   Rscript dev/ci701-rev-d3.R
.libPaths(c("C:/Users/adf44/source/r/wt-ci701-lib",
            "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
cat("lib:", find.package("frmtmb"), "\n")
orig <- get("se_edge_cor", ns)
stub <- function(fit, np) logical(np)

corr_at <- function(fit, th) {
  out <- list()
  for (bk in fit$frame$re_blocks) {
    d <- bk$dim
    if (!identical(bk$covstruct, "us") || d < 2) next
    tc <- th[bk$theta_idx[-seq_len(d)]]
    C <- ns$us_chol_cor(tc, d)
    out[[length(out) + 1]] <- C[lower.tri(C)]
  }
  unlist(out)
}

run <- function(tag, expr_fit) {
  for (a in c("lane", "noD2")) {
    assignInNamespace("se_edge_cor", if (a == "lane") orig else stub, ns)
    w <- character(0); bm <- character(0)
    fit <- withCallingHandlers(expr_fit(), warning = function(x) {
      w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning")
    }, message = function(m) {
      if (inherits(m, "frmtmb_boundary_fit")) bm <<- c(bm, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
    p <- fit$opt$par
    th <- fit$estimates$theta
    np <- length(p)
    ec <- orig(fit, np)
    pulled <- ns$se_pull_back(p, ec)
    thp <- pulled[names(p) == "theta"]
    g <- max(abs(drop(fit$obj$gr(p))))
    lost <- ns$sdr_of(fit)$se_lost
    an <- ns$se_fit_analysis(fit, force = TRUE)
    se <- suppressWarnings(sqrt(diag(an$shown)))
    thi <- which(names(p) == "theta")
    cat(sprintf("RES %s %s code=%d g=%.3g nwarn=%d bmsg=%d\n", tag, a,
                fit$opt$convergence, g, length(w), length(bm)))
    cat("  theta   :", sprintf("%.5g", th), "\n")
    cat("  edge_cor:", which(ec), "\n")
    cat("  corr fit:", sprintf("%.5f", corr_at(fit, th)), "\n")
    cat("  corr pul:", sprintf("%.5f", corr_at(fit, thp)), "\n")
    cat("  fit lost:", if (length(lost)) paste(names(lost), lost, sep = ":",
                                               collapse = ",") else "none",
        "\n")
    cat("  an  lost:", if (length(an$lost)) paste(names(an$lost), an$lost,
                                                  sep = ":", collapse = ",")
        else "none", "\n")
    cat("  se all  :", sprintf("%.5g", se), "\n")
    cat("  loglik  :", sprintf("%.8f", as.numeric(logLik(fit))), "\n")
    if (length(w)) cat("  warn:", substr(w, 1, 120), sep = "\n   ")
    if (length(bm)) cat("  bmsg:", substr(bm, 1, 160), sep = "\n   ")
    cat("\n")
  }
  assignInNamespace("se_edge_cor", orig, ns)
}

# test-v16.R: us(spp + 0 | site) on 3 species of a rank-2 covariance
set.seed(61)
n_site <- 60; S <- 6
L <- matrix(0, S, 2)
L[, 1] <- c(1.2, 0.8, 0.5, -0.4, 0.9, 0.2)
L[2:S, 2] <- c(0.7, -0.6, 0.4, 0.3, -0.5)
f <- matrix(rnorm(n_site * 2), n_site, 2)
u <- f %*% t(L)
dd <- data.frame(y = rpois(n_site * S, exp(0.5 + as.vector(t(u)))),
                 spp = factor(rep(seq_len(S), n_site)),
                 site = factor(rep(seq_len(n_site), each = S)))
d2 <- dd[dd$spp %in% c("1", "2", "3"), , drop = FALSE]
d2$spp <- droplevels(d2$spp)
v16 <- function(st) run(sprintf("v16us st=%g", st[5]), function()
  frm(bf(y ~ 1 + us(spp + 0 | site)) + poisson(), data = d2, start = list(theta = st)))

v16(c(0.091939, 0.069813, -0.15009, 0.26129, 26368, -14159))
v16(c(0.091939, 0.069813, -0.15009, 0.26129, 1318.4 * 100, -707.94 * 100))
mkB <- function(s) {
  set.seed(s)
  ng <- 30; m <- 10
  g <- factor(rep(seq_len(ng), each = m))
  d <- data.frame(g = g, x = rnorm(ng * m), z = rnorm(ng * m))
  uu <- rnorm(ng); vv <- rnorm(ng)
  d$y <- 1 + 0.5 * d$x - 0.3 * d$z + 0.8 * uu[g] + 0.8 * vv[g] * d$x +
    0.8 * (uu[g] + vv[g]) / sqrt(2) * d$z + rnorm(ng * m, 0, 0.5)
  d
}
for (s in 1:2) {
  dB <- mkB(s)
  for (T in c(2e4, 1e5)) {
    run(sprintf("caseB s=%d T=%g", s, T), function()
      frm(bf(y ~ x + z + (1 + x + z | g)), family = gaussian(), data = dB,
          start = list(theta = c(-0.2, -0.2, -0.1, 0.2, 0.75 * T,
                                 0.5 * T))))
  }
}
