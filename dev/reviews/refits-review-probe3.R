# Review of lane refits: 3434 influence() smooths (s, by, fs, t2, te, gp),
# cost; 3188 attribute size, serialization, check against nnet::multinom;
# 2939 names and the exact marginal at the estimates.
source("C:/Users/adf44/source/r/frmtmb-wt-refits/dev/reviews/refits-review-lib.R")

cat("\n==== 3434 influence() smooths ====\n")
set.seed(8)
n <- 60
dl <- data.frame(x = runif(n, 0, 3), z = rnorm(n), v = runif(n),
                 f = factor(rep(c("a", "b"), n / 2)),
                 g = factor(rep(1:5, 12)), blk = factor(rep(1:6, 10)))
dl$y <- sin(dl$x) + 0.3 * dl$z + dl$v^2 + rnorm(n, 0, 0.3)
forms <- list(
  s = y ~ z + s(x, k = 5),
  s_cr = y ~ z + s(x, bs = "cr", k = 5),
  s_by = y ~ z + f + s(x, by = f, k = 5),
  s_fs = y ~ z + s(x, g, bs = "fs", k = 4),
  t2 = y ~ z + t2(x, v, k = 4),
  te = y ~ z + te(x, v, k = 4),
  gp_exact = y ~ z + gp(x),
  gp_approx = y ~ z + gp(x, k = 8)
)
res <- list()
for (nm in names(forms)) {
  fit <- go(paste(nm, "fit"), {
    f <- suppressMessages(frm(forms[[nm]], data = dl)); class(f)[1]; f })
  if (is.character(fit)) next
  t0 <- proc.time()[["elapsed"]]
  inf <- go(paste(nm, "influence(groups = blk)"),
            { i <- suppressWarnings(influence(fit, groups = "blk", data = dl))
              round(i$fixed[, 1:min(3, ncol(i$fixed)), drop = FALSE], 6) })
  t1 <- proc.time()[["elapsed"]]
  cat(nm, "influence elapsed s:", round(t1 - t0, 2), "\n")
  i2 <- tryCatch(suppressWarnings(influence(fit, groups = "blk", data = dl)),
                 error = function(e) NULL)
  if (!is.null(i2)) {
    res[[nm]] <- list(fixed = i2$fixed, theta = i2$theta)
    cat(nm, "theta table row 1:", sprintf("%.6f", i2$theta[1, ]), "\n")
    cat(nm, "theta full:", sprintf("%.6f", i2$theta_full), "\n")
  }
}
saveRDS(res, paste0("C:/Users/adf44/source/r/frmtmb-wt-refits/dev/reviews/refits-review-infl-", arm, ".rds"))

cat("\n== gp exact: delete the row with the largest x (dmax moves) ==\n")
fg <- suppressMessages(frm(y ~ z + gp(x), data = dl))
imax <- which.max(dl$x)
dl$one <- factor(seq_len(n) == imax)
ig <- suppressWarnings(influence(fg, groups = "one", data = dl))
man <- suppressMessages(frm(y ~ z + gp(x), data = dl[-imax, ]))
go("gp exact deletion: table vs manual refit on n-1 rows",
   list(table_fixed = ig$fixed["TRUE", ], manual_fixed = fixef(man, flatten = TRUE)[1:2],
        table_theta = ig$theta["TRUE", ], full_theta = ig$theta_full))
go("gp exact VarCorr-scale lscale of manual vs full",
   list(manual = VarCorr(man), full = VarCorr(fg)))

cat("\n==== 3188 attribute ====\n")
set.seed(5)
dm <- data.frame(x = rnorm(200), N = 10L)
pr <- cbind(1, exp(0.3 + 0.5 * dm$x), exp(-0.2 * dm$x))
pr <- pr / rowSums(pr)
dm$Y <- t(sapply(1:200, function(i) rmultinom(1, 10, pr[i, ])))
colnames(dm$Y) <- c("a", "b", "c")
fm <- frm(Y | trials(N) ~ x, family = multinomial(), data = dm)
md <- attr(fm, "marginaleffects_modeldata")
go("attribute columns", if (is.null(md)) "NULL" else names(md))
go("sizes (bytes): fit serialized, attribute serialized, frame data_frame",
   c(fit = length(serialize(fm, NULL)),
     attr = length(serialize(md, NULL)),
     data_frame = length(serialize(fm$frame$data_frame, NULL))))
go("saveRDS/readRDS keeps the attribute",
   { tf <- tempfile(fileext = ".rds"); saveRDS(fm, tf); r <- readRDS(tf)
     identical(attr(r, "marginaleffects_modeldata"), md) })
go("refit(): attribute response is the old one",
   { rf <- refit(fm, fm$frame$y[[1]][sample(200), ])
     md2 <- attr(rf, "marginaleffects_modeldata")
     c(identical(md2$Y.a, md$Y.a), identical(rf$frame$y[[1]][, 1], md$Y.a)) })
go("predictions vs nnet::multinom (probabilities x N)", {
  p <- marginaleffects::predictions(fm)
  pd <- as.data.frame(p)
  mn <- nnet::multinom(Y ~ x, data = dm, trace = FALSE)
  pn <- as.data.frame(marginaleffects::predictions(mn, type = "probs"))
  ord_f <- order(pd$group, pd$rowid); ord_n <- order(pn$group, pn$rowid)
  c(max_rel_est = max(abs(pd$estimate[ord_f] / 10 - pn$estimate[ord_n]) /
                        pn$estimate[ord_n]),
    max_rel_se = max(abs(pd$std.error[ord_f] / 10 - pn$std.error[ord_n]) /
                       pn$std.error[ord_n]),
    groups_f = paste(unique(pd$group), collapse = ","),
    groups_n = paste(unique(pn$group), collapse = ","))
})
go("avg_predictions vs nnet::multinom", {
  a <- as.data.frame(marginaleffects::avg_predictions(fm))
  b <- as.data.frame(marginaleffects::avg_predictions(nnet::multinom(
    Y ~ x, data = dm, trace = FALSE), type = "probs"))
  cbind(a[, c("group", "estimate", "std.error")],
        nn_est_x10 = b$estimate * 10, nn_se_x10 = b$std.error * 10)
})
go("a large fit: bytes of the attribute against the fit", {
  set.seed(6)
  nb <- 20000
  db <- data.frame(x = rnorm(nb), N = 10L)
  prb <- cbind(1, exp(0.3 + 0.5 * db$x), exp(-0.2 * db$x)); prb <- prb / rowSums(prb)
  db$Y <- t(vapply(seq_len(nb), function(i) rmultinom(1, 10, prb[i, ])[, 1], numeric(3)))
  colnames(db$Y) <- c("a", "b", "c")
  fb <- frm(Y | trials(N) ~ x, family = multinomial(), data = db)
  c(fit_MB = length(serialize(fb, NULL)) / 1e6,
    attr_MB = length(serialize(attr(fb, "marginaleffects_modeldata"), NULL)) / 1e6)
})

cat("\n==== 2939 names and the exact marginal ====\n")
set.seed(11)
ng <- 12
A <- 0.5 ^ abs(outer(1:ng, 1:ng, "-"))
dimnames(A) <- list(paste0("g", 1:ng), paste0("g", 1:ng))
d <- data.frame(g = factor(rep(paste0("g", 1:ng), each = 10),
                           levels = paste0("g", 1:ng)),
                f = factor(rep(c("a", "b", "c"), 40)),
                f4 = factor(rep(c("a", "b", "c", "d"), 30)),
                x = rnorm(120), z = rnorm(120))
u <- MASS::mvrnorm(3, rep(0, ng), A)
d$y <- 1 + 0.5 * d$x + u[1, as.integer(d$g)] +
  (d$f == "b") * u[2, as.integer(d$g)] + 0.3 * d$z * u[3, as.integer(d$g)] +
  rnorm(120, 0, 0.7)
exact_ll <- function(fit, cols) {
  vc <- VarCorr(fit)
  sds <- vc$g$sd[, "Estimate"]
  sig <- vc$residual__$sd[, "Estimate"]
  V <- diag(sig^2, nrow(d))
  gi <- as.integer(d$g)
  for (k in seq_along(cols)) {
    Z <- matrix(0, nrow(d), ng); Z[cbind(seq_len(nrow(d)), gi)] <- cols[[k]]
    V <- V + sds[k]^2 * Z %*% A %*% t(Z)
  }
  X <- model.matrix(~ x, d)
  mu <- drop(X %*% fixef(fit)[, 1])
  R <- chol(V)
  e <- backsolve(R, d$y - mu, transpose = TRUE)
  -sum(log(diag(R))) - 0.5 * sum(e^2) - nrow(d) / 2 * log(2 * pi)
}
fc <- frm(y ~ x + (1 + z + f || gr(g, cov = A)), data = d, data2 = list(A = A))
go("mix: VarCorr names", rownames(VarCorr(fc)$g$sd))
go("mix: sd variables", grep("^sd_|^cor_", variables(fc), value = TRUE))
go("mix: ranef dims", lapply(ranef(fc), dim))
go("mix: ranef coef names", dimnames(ranef(fc)$g)[[3]])
go("mix: coef(fit)$g columns", dimnames(coef(fc)$g)[[3]])
go("mix: r_ variables head/tail", { v <- grep("^r_", variables(fc), value = TRUE); c(head(v, 3), tail(v, 3), length(v)) })
go("mix: logLik - exact at the estimates",
   as.numeric(logLik(fc)) - exact_ll(fc, list(1, d$z, d$f == "b", d$f == "c")))
go("mix: optim from the estimates cannot raise the exact marginal", {
  vc <- VarCorr(fc)
  p0 <- c(fixef(fc)[1:2, 1], log(vc$g$sd[, 1]), log(vc$residual__$sd[, 1]))
  cols <- list(1, d$z, d$f == "b", d$f == "c")
  gi <- as.integer(d$g); X <- model.matrix(~ x, d)
  nll <- function(p) {
    V <- diag(exp(2 * p[7]), nrow(d))
    for (k in 1:4) { Z <- matrix(0, nrow(d), ng); Z[cbind(seq_len(nrow(d)), gi)] <- cols[[k]]
      V <- V + exp(2 * p[2 + k]) * Z %*% A %*% t(Z) }
    R <- chol(V); e <- backsolve(R, d$y - drop(X %*% p[1:2]), transpose = TRUE)
    sum(log(diag(R))) + 0.5 * sum(e^2) + nrow(d) / 2 * log(2 * pi)
  }
  o <- optim(p0, nll, method = "BFGS", control = list(reltol = 1e-14, maxit = 500))
  c(at_fit = -nll(p0), optim = -o$value, gain = nll(p0) - o$value)
})
f4 <- frm(y ~ x + (0 + f4 || gr(g, cov = A)), data = d, data2 = list(A = A))
go("f4: logLik - exact at the estimates",
   as.numeric(logLik(f4)) - exact_ll(f4, lapply(levels(d$f4), function(l) d$f4 == l)))
go("f4: VarCorr names", rownames(VarCorr(f4)$g$sd))
go("f4: summary prints", { capture.output(summary(f4)) |> grep(pattern = "sd\\(", value = TRUE) })
go("f4: confint_varcorr rows", rownames(confint_varcorr(f4)))
go("f4: predict new level allow_new_levels",
   { nd <- data.frame(g = "gnew", f4 = factor("b", levels = levels(d$f4)), x = 0)
     predict(f4, newdata = nd, allow_new_levels = TRUE) })
go("f4: ranef(f4) condVar present", { r <- ranef(f4); dim(r$g) })
go("f4: re_formula = ~ (0 + f4 || gr(g, cov = A))",
   { a <- fitted(f4, newdata = d[1:5, ], re_formula = ~ (0 + f4 || gr(g, cov = A)))[, 1]
     b <- fitted(f4, newdata = d[1:5, ])[, 1]; max(abs(a - b)) })
