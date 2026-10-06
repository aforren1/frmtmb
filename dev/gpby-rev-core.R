# Reviewer checks on the lane build, no Stan: data_gp parity, class
# "sd" on a gp(), extra_cov identities, PSD and cost, summary()$gp
# against variables().
.libPaths(c("C:/Users/adf44/source/r/wt-gpby-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb)})
cat("lib:", find.package("frmtmb"), "\n")
wt <- "C:/Users/adf44/source/r/frmtmb-wt-gpby"
src <- parse(file.path(wt, "tests/testthat/test-gp-by.R"))
for (e in src) if (is.call(e) && identical(e[[1]], as.name("<-"))) eval(e)
d <- gpby_data()
d$z <- { set.seed(77); stats::runif(nrow(d), 0, 3) }
ulps <- function(a, b) {
  a <- as.numeric(a); b <- as.numeric(b)
  max(abs(a - b) / pmax(.Machine$double.eps * pmax(abs(a), abs(b)),
                        .Machine$double.xmin))
}

cat("\n### B1 data_gp parity\n")
sd_of <- function(f, data = d) {
  suppressMessages(brms::standata(brms::bf(f), data = data,
                                  family = gaussian()))
}
par_exact <- function(fs, data = d) {
  f <- stats::as.formula(fs)
  sdat <- sd_of(f, data)
  fr <- frm(bf(f), data = data, dry_run = "frame")
  gis <- fr$linpreds[["y.mu"]]$gps
  kgp <- as.integer(sdat$Kgp_1)
  byfac <- !is.null(sdat[["Igp_1_1"]])
  for (j in seq_len(kgp)) {
    sfx <- if (byfac) paste0("_1_", j) else "_1"
    bk <- fr$re_blocks[[gis[[j]]$block_id]]
    Xb <- as.matrix(sdat[[paste0("Xgp", sfx)]])
    pos <- gis[[j]]$positions
    div <- bk$gp_lscale_div
    dmb <- sdat[[paste0("dmax", sfx)]]
    # brms orders distinct positions by first appearance; sort both
    ob <- do.call(order, as.data.frame(Xb))
    Xs <- Xb[ob, , drop = FALSE]
    Ig <- if (byfac) as.integer(sdat[[paste0("Igp", sfx)]]) else
      seq_len(nrow(data))
    Z <- as.matrix(fr$linpreds[["y.mu"]]$Z)[, bk$c_idx, drop = FALSE]
    rows <- which(rowSums(abs(Z)) > 0)
    cat(sprintf(paste0("EXACT %-40s sub %d | dmax identical %s (ulps %.1f)",
                       " | pos/div vs Xgp identical %s ulps %.1f | ",
                       "npos %d vs %d | Igp identical %s\n"),
                fs, j, identical(unname(div), unname(as.numeric(dmb))),
                ulps(div, dmb),
                identical(unname(pos / div), unname(Xs)),
                ulps(pos / div, Xs), nrow(pos), nrow(Xb),
                identical(rows, Ig)))
  }
}
for (fs in c("y ~ gp(x, by = f)", "y ~ gp(x, by = f, cmc = FALSE)",
             "y ~ gp(x, by = w)", "y ~ gp(x, z, by = f)",
             "y ~ gp(x, by = f, scale = FALSE)",
             "y ~ gp(x, by = f, gr = FALSE)")) {
  r <- try(par_exact(fs))
}
par_hsgp <- function(fs, data = d) {
  f <- stats::as.formula(fs)
  sdat <- sd_of(f, data)
  fr <- frm(bf(f), data = data, dry_run = "frame")
  gis <- fr$linpreds[["y.mu"]]$gps
  Zall <- as.matrix(fr$linpreds[["y.mu"]]$Z)
  kgp <- as.integer(sdat$Kgp_1)
  byfac <- !is.null(sdat[["Igp_1_1"]])
  for (j in seq_len(kgp)) {
    sfx <- if (byfac) paste0("_1_", j) else "_1"
    gi <- gis[[j]]
    bk <- fr$re_blocks[[gi$block_id]]
    Ig <- if (byfac) as.integer(sdat[[paste0("Igp", sfx)]]) else
      seq_len(nrow(data))
    Jg <- sdat[[paste0("Jgp", sfx)]]
    Xb <- as.matrix(sdat[[paste0("Xgp", sfx)]])
    if (!is.null(Jg)) Xb <- Xb[Jg, , drop = FALSE]
    Cg <- sdat[[paste0("Cgp", sfx)]]
    if (!is.null(Cg)) Xb <- Xb * as.numeric(Cg)
    Zf <- Zall[Ig, bk$c_idx, drop = FALSE]
    sl <- as.matrix(sdat[[paste0("slambda", sfx)]])
    cat(sprintf(paste0("HSGP  %-40s sub %d | basis identical %s ulps %.1f ",
                       "maxabs %.2e | slambda identical %s ulps %.1f | ",
                       "NBgp %d vs %d\n"),
                fs, j, identical(unname(Zf), unname(Xb)), ulps(Zf, Xb),
                max(abs(Zf - Xb)),
                identical(unname(as.vector(gi$omega)), as.vector(sl)),
                ulps(gi$omega, sl), as.integer(sdat$NBgp_1),
                nrow(gi$omega)))
  }
}
for (fs in c("y ~ gp(x, by = f, k = 8)", "y ~ gp(x, by = f, k = 8, cmc = FALSE)",
             "y ~ gp(x, by = w, k = 8)", "y ~ gp(x, z, by = f, k = 5)",
             "y ~ gp(x, by = f, k = 8, scale = FALSE)",
             "y ~ gp(x, by = f, k = 8, gr = FALSE)",
             "y ~ gp(x, by = f, k = 8, c = 2)")) {
  r <- try(par_hsgp(fs))
}
cat("standata names exact:", paste(names(sd_of(y ~ gp(x, by = f))),
                                    collapse = " "), "\n")

cat("\n### B2 class sd on gp\n")
d$g <- factor(rep(1:6, length.out = nrow(d)))
msg <- function(expr) tryCatch({expr; "NO ERROR"},
                               error = function(e) conditionMessage(e))
cat("sd, gp only:", msg(frm(bf(y ~ gp(x)), data = d,
                            prior = set_prior("exponential(1)",
                                              class = "sd"))), "\n")
cat("sd, group gp(x):", msg(frm(bf(y ~ gp(x)), data = d,
                                prior = set_prior("exponential(1)",
                                                  class = "sd",
                                                  group = "gp(x)"))), "\n")
cat("sdgp with group:", msg(frm(bf(y ~ gp(x)), data = d,
                                prior = set_prior("exponential(1)",
                                                  class = "sdgp",
                                                  group = "g"))), "\n")
f0 <- frm(bf(y ~ gp(x) + (1 | g)), data = d)
f1 <- frm(bf(y ~ gp(x) + (1 | g)), data = d,
          prior = set_prior("exponential(50)", class = "sd"))
f2 <- frm(bf(y ~ gp(x) + (1 | g)), data = d,
          prior = set_prior("exponential(50)", class = "sdgp"))
cat("gp + (1|g): sdgp flat", exp(f0$estimates$theta[1]),
    "| class sd exp(50)", exp(f1$estimates$theta[1]),
    "| class sdgp exp(50)", exp(f2$estimates$theta[1]), "\n")
print(variables(f0))
gpf <- get_prior(bf(y ~ gp(x, by = f) + (1 | g)), data = d)
print(gpf[gpf$class %in% c("sd", "sdgp", "lscale"),
          c("prior", "class", "coef", "group", "dpar")])
bp <- brms::default_prior(brms::bf(y ~ gp(x, by = f) + (1 | g)), data = d)
print(bp[bp$class %in% c("sd", "sdgp", "lscale"),
         c("prior", "class", "coef", "group", "dpar")])
gps <- get_prior(bf(y ~ x, sigma ~ gp(x, by = f)), data = d)
print(gps[gps$class %in% c("sdgp", "lscale"),
          c("prior", "class", "coef", "dpar")])
bps <- brms::default_prior(brms::bf(y ~ x, sigma ~ gp(x, by = f)), data = d)
print(bps[bps$class %in% c("sdgp", "lscale"),
          c("prior", "class", "coef", "dpar")])

cat("\n### B3 extra_cov\n")
fit <- frm(bf(y ~ gp(x)), data = d)
set.seed(3)
grids <- list(
  extrap = data.frame(x = seq(5, 9, length.out = 30)),
  repeated = data.frame(x = rep(c(6.3, 7.1, 8.4), each = 4)),
  near = data.frame(x = c(7, 7 + 1e-12, 7 + 1e-9, 7 + 1e-6, 0.1 + 0.2,
                          0.3, 8, 8)),
  mixed_obs = data.frame(x = c(d$x[1:5], 6.55, 6.55, 7.3)))
for (nm in names(grids)) {
  nd <- grids[[nm]]
  lb <- frm_lp_basis(fit, newdata = nd, extra_cov = TRUE)
  E <- lb$extra_cov
  Es <- E - lb$extra_white
  ev <- eigen(E, symmetric = TRUE, only.values = TRUE)$values
  evs <- eigen(Es, symmetric = TRUE, only.values = TRUE)$values
  lb0 <- frm_lp_basis(fit, newdata = nd)
  cat(sprintf(paste0("GRID %-9s n %d | identical(diag, extra_var) %s | ",
                     "identical extra_var w/o extra_cov %s | symmetric %s |",
                     " min eig E %.3e (max %.3e) | min eig E-white %.3e |",
                     " full Sigma min eig %.3e\n"),
              nm, nrow(nd), identical(diag(E), lb$extra_var),
              identical(lb$extra_var, lb0$extra_var), isSymmetric(E),
              min(ev), max(ev), min(evs),
              min(eigen(as.matrix(lb$A) %*% lb$V %*% t(as.matrix(lb$A)) + E,
                        symmetric = TRUE, only.values = TRUE)$values)))
}
# closed form of the kriging covariance at the fitted kernel
th <- fit$estimates$theta
bk <- Filter(function(b) b$covstruct == "gp", fit$frame$re_blocks)[[1]]
pos <- fit$frame$linpreds[["y.mu"]]$gps[[1]]$positions[, 1]
kf <- function(a, b) {
  exp(2 * th[1]) * (exp(-outer(a, b, "-")^2 / (2 * exp(2 * th[2]))) +
                      1e-6 * (outer(a, b, "-") == 0))
}
xg <- grids$extrap$x
Sref <- kf(xg, xg) - kf(xg, pos) %*% solve(kf(pos, pos), kf(pos, xg))
lb <- frm_lp_basis(fit, newdata = grids$extrap, extra_cov = TRUE)
cat(sprintf("closed-form kriging cov: max rel diff %.3e\n",
            max(abs(lb$extra_cov - Sref)) / max(abs(Sref))))

cat("\n### B3b cost on a 2000-row grid\n")
big <- data.frame(x = seq(-1, 9, length.out = 2000))
tm <- function(expr, k = 3) {
  min(vapply(seq_len(k), function(i) {
    system.time(expr)[["elapsed"]]
  }, 0))
}
t0 <- tm(frm_lp_basis(fit, newdata = big))
t1 <- tm(frm_lp_basis(fit, newdata = big, extra_cov = TRUE))
g0 <- gc(reset = TRUE)
lb <- frm_lp_basis(fit, newdata = big, extra_cov = TRUE)
g1 <- gc()
cat(sprintf(paste0("2000 rows: frm_lp_basis %.3f s, extra_cov = TRUE %.3f ",
                   "s; max Vcells used during call %.1f MB; object %.1f MB",
                   "\n"),
            t0, t1, g1[2, 6] , as.numeric(object.size(lb)) / 2^20))
Rprof(tmp <- tempfile(), memory.profiling = TRUE)
lb <- frm_lp_basis(fit, newdata = big, extra_cov = TRUE)
Rprof(NULL)
sp <- summaryRprof(tmp, memory = "both")
print(utils::head(sp$by.total[, c("total.time", "mem.total")], 15))
fitb <- frm(bf(y ~ gp(x, by = f)), data = d)
bigf <- data.frame(x = rep(seq(-1, 9, length.out = 700), 3),
                   f = factor(rep(c("a", "b", "c"), each = 700),
                              levels = levels(d$f)))
t2 <- tm(frm_lp_basis(fitb, newdata = bigf, extra_cov = TRUE), 2)
cat(sprintf("2100 rows, 3 sub-GPs: extra_cov = TRUE %.3f s\n", t2))
lbf <- frm_lp_basis(fitb, newdata = bigf, extra_cov = TRUE)
cat("cross-level block all zero:",
    all(lbf$extra_cov[1:700, 701:2100] == 0), "\n")

cat("\n### B4 summary()$gp vs gp_brms_values\n")
for (fs in c("y ~ gp(x, by = f)", "y ~ gp(x, by = f, k = 8)",
             "y ~ gp(x, z, by = f, iso = FALSE)",
             "y ~ gp(x, z, by = f, k = 5, iso = FALSE)",
             "y ~ gp(x, by = f, scale = FALSE)")) {
  ff <- frm(bf(stats::as.formula(fs)), data = d)
  sg <- summary(ff)$gp
  vals <- frmtmb:::gp_brms_values(ff, ff$estimates$theta)
  est <- sg[, 1]
  nmv <- sub("^(sdgp|lscale)_", "\\1(", names(vals))
  nmv <- paste0(nmv, ")")
  cat(sprintf("%-42s names equal %s | max rel diff %.2e | %s\n", fs,
              identical(rownames(sg), nmv),
              max(abs(est / vals - 1)),
              paste(rownames(sg), collapse = " ")))
}
cat("DONE\n")
