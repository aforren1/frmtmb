# REVIEW re-check after punch round 1.
#
#   Rscript dev/arcovsample-rev-14-recheck.R
#
# (a) integrity of the REBUILT lane library for the files punch 1 moved;
# (b) the corrected wording, RENDERED, including the two-argument \eqn;
# (c) the worker's rule measured the way the rule states it: not only
#     "is position k shifted" but "does it carry lag i", which is tested
#     by changing coefficient i and seeing whether the row moves;
# (d) a 1-row group and a 3-row group at p = 3 through log_lik();
# (e) where the OLD Student-t form actually becomes NaN, against the
#     3.6e305 the new NEWS entry names.

LANE <- "C:/Users/adf44/source/r/wt-arcovsample-lib"
.libPaths(c(LANE, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
suppressMessages(library(frmtmb.sample))
stopifnot("arma_cond_resp" %in% getNamespaceExports("frmtmb"))
WT <- "C:/Users/adf44/source/r/frmtmb-wt-arcovsample"
`%||%` <- function(a, b) if (is.null(a)) b else a

cat("=== (a) rebuilt install against the worktree source\n")
src_fun <- function(path, name) {
  for (e in parse(file = path)) {
    if (is.call(e) && length(e) >= 3L &&
        as.character(e[[1L]])[1L] %in% c("<-", "=") &&
        identical(as.character(e[[2L]]), name)) {
      return(eval(e[[3L]], envir = new.env()))
    }
  }
  NULL
}
cmp <- function(pkg, path, name) {
  inst <- get0(name, envir = asNamespace(pkg), inherits = FALSE)
  s <- src_fun(file.path(WT, path), name)
  ok <- !is.null(inst) && !is.null(s) &&
    identical(deparse(inst), deparse(s))
  cat("  ", if (ok) "same" else "DIFFERS", "  ", pkg, "::", name, "\n",
      sep = "")
  ok
}
ok <- c(cmp("frmtmb", "R/autocor.R", "arma_cond_resp"),
        cmp("frmtmb", "R/autocor.R", "arma_cond_dpars"),
        cmp("frmtmb", "R/predict-brms.R", "rescor_row_loglik"),
        cmp("frmtmb.sample", "extensions/frmtmb.sample/R/loo.R",
            "draws_row_loglik"),
        cmp("frmtmb.sample", "extensions/frmtmb.sample/R/loo.R",
            "draws_loglik_factors"))
cat("  all same: ", all(ok), "\n", sep = "")
# R/compat.R carries the sentence as DATA, so compare the string
compat_src <- paste(readLines(file.path(WT, "R/compat.R"), warn = FALSE),
                    collapse = "\n")
cat("  R/compat.R source has the new rule: ",
    grepl("lag i first reaches within-group position i + 1",
          compat_src, fixed = TRUE), "\n", sep = "")

cat("\n=== (b) the corrected wording, rendered\n")
outdir <- file.path(WT, "dev/arcovsample-rev-log/rd2")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
for (f in c("extensions/frmtmb.sample/man/sample-log_lik.Rd",
            "man/frmtmb-autocor.Rd")) {
  nm <- basename(f)
  txt <- capture.output(tools::Rd2txt(file.path(WT, f), out = stdout(),
                                      package = "x"))
  writeLines(txt, file.path(outdir, sub("\\.Rd$", ".txt", nm)))
  one <- paste(txt, collapse = "\n")
  cat("\n---- ", nm, " (", length(txt), " rendered lines)\n", sep = "")
  # the retracted claim must be GONE, in every spelling
  gone <- c("first max(p, q) rows", "first rows get no lagged",
            "first rows of a group get no lagged",
            "the first rows of a group get no lagged")
  for (g in gone) {
    cat("    retracted '", g, "': ", sum(grepl(g, txt, fixed = TRUE)),
        "\n", sep = "")
  }
  keep <- c("lag", "position", "FIRST row", "max(p, q) + 1")
  for (k in keep) {
    cat("    present  '", k, "': ", sum(grepl(k, txt, fixed = TRUE)),
        "\n", sep = "")
  }
  # a two-argument \eqn renders its SECOND argument in text; a broken
  # one leaves a brace or the latex
  cat("    stray brace or backslash in rendered text: ",
      sum(grepl("[\\\\]|[{}]", txt)), "\n", sep = "")
  cat("    em dashes ", sum(grepl("\u2014", txt, fixed = TRUE)),
      "  lines over 80 cols ", sum(nchar(txt) > 80), "\n", sep = "")
}
cat("\n-- the rule as rendered, frmtmb-autocor:\n")
t2 <- readLines(file.path(outdir, "frmtmb-autocor.txt"), warn = FALSE)
i <- grep("lag", t2, fixed = TRUE)
i <- i[i > 100][1L]
cat(paste(t2[grep("CONDITIONAL likelihood", t2)[1L] +
                seq(0, 14) - 1L], collapse = "\n"), "\n")
cat("\n-- the rule as rendered, sample-log_lik:\n")
t1 <- readLines(file.path(outdir, "sample-log_lik.txt"), warn = FALSE)
j <- grep("one-step conditional mean", t1, fixed = TRUE)[1L]
cat(paste(t1[j:(j + 11L)], collapse = "\n"), "\n")

cat("\n=== (c) does position k carry ONLY the lags up to k - 1?\n")
# Falsification: change coefficient i alone. If the rule is right, a row
# at within-group position k must NOT move when coefficient i > k - 1
# changes, and MUST move when i <= k - 1.
set.seed(6161L)
ng <- 4L
dd <- do.call(rbind, lapply(seq_len(ng), function(g) {
  data.frame(g = factor(g, levels = seq_len(ng)), t = seq_len(7L))
}))
dd$x <- rnorm(nrow(dd))
dd$y <- 0.5 + 0.4 * dd$x + as.numeric(arima.sim(list(ar = 0.5),
                                                nrow(dd)))
dd <- dd[sample(nrow(dd)), ]
rownames(dd) <- NULL
ord <- order(dd$g, dd$t)
pos <- unlist(lapply(split(seq_along(ord), dd$g[ord]), seq_along))

fit3 <- frm(bf(y ~ x + ar(t, g, p = 3)), family = gaussian(), data = dd,
            dry_run = "objective")
mu_at <- function(fit, a) {
  lab <- c(frmtmb::brms_par_labels(fit), "lp__")
  m <- matrix(0, 1L, length(lab), dimnames = list(NULL, lab))
  m[, "b_Intercept"] <- 0.5; m[, "b_x"] <- 0.4; m[, "sigma"] <- 0.8
  for (k in seq_along(a)) m[, paste0("thetaac_", k)] <- a[k]
  ds <- structure(list(stanfit = NULL, draws = m, fit = fit),
                  class = "frmtmb_draws")
  sh <- frmtmb.sample:::draws_fit_at(
    ds, 1L, frmtmb.sample:::draws_par_index(fit))
  as.numeric(frmtmb::arma_cond_dpars(sh, frmtmb::eval_dpars(sh))$y$mu)
}
base <- c(0.30, 0.20, 0.15)
m0 <- mu_at(fit3, base)
cat("  ar(p = 3): rows that MOVE when coefficient i alone changes\n")
cat(sprintf("  %-14s %s\n", "coefficient",
            paste(sprintf("%3d", 1:7), collapse = " ")))
cat(sprintf("  %-14s %s\n", "position",
            paste(sprintf("%3d", pos[1:7]), collapse = " ")))
for (i in 1:3) {
  a <- base; a[i] <- a[i] + 0.11
  d <- abs(mu_at(fit3, a) - m0)[ord]
  cat(sprintf("  lag %-10d %s\n", i,
              paste(sprintf("%3s", ifelse(d[1:7] > 1e-12, "Y", ".")),
                    collapse = " ")))
  # the rule, as a single assertion over ALL rows
  want <- pos >= i + 1L
  cat("    moves exactly at positions >= ", i + 1L, ": ",
      identical(unname(d > 1e-12), unname(want)), "\n", sep = "")
}

cat("\n=== (d) a 1-row group and a 3-row group at p = 3, through log_lik()\n")
re_prior_of <- function(sh) {
  tot <- 0
  for (bk in sh$frame[["re_blocks"]] %||% list()) {
    f <- frmtmb:::covstruct_registry[[bk$covstruct]]$nll
    tot <- tot + as.numeric(f(sh$estimates[["b"]][bk$b_idx],
                              sh$estimates[["theta"]][bk$theta_idx], bk))
  }
  tot
}
set.seed(6162L)
lens <- c(9L, 7L, 3L, 1L, 6L)
de <- do.call(rbind, lapply(seq_along(lens), function(i) {
  data.frame(g = factor(i, levels = seq_along(lens)),
             t = seq_len(lens[i]))
}))
de$x <- rnorm(nrow(de))
de$y <- 0.5 + 0.4 * de$x + rnorm(nrow(de), 0, 0.8)
de <- de[sample(nrow(de)), ]
rownames(de) <- NULL
cat("  N = ", nrow(de), "  group sizes ", paste(lens, collapse = " "),
    "  seed 6162\n", sep = "")
run <- function(lab, form, np) {
  f <- frm(form, family = gaussian(), data = de, dry_run = "objective")
  lab2 <- c(frmtmb::brms_par_labels(f), "lp__")
  m <- matrix(NA_real_, 3L, length(lab2), dimnames = list(NULL, lab2))
  m[, "b_Intercept"] <- c(0.5, 0.6, 0.45)
  m[, "b_x"] <- c(0.4, 0.35, 0.5)
  m[, "sigma"] <- c(0.8, 0.9, 0.7)
  for (k in seq_len(np)) m[, paste0("thetaac_", k)] <- c(0.3, -0.2, 0.15) /
    k
  m[, "lp__"] <- 0
  stopifnot(!anyNA(m))
  ds <- structure(list(stanfit = NULL, draws = m, fit = f),
                  class = "frmtmb_draws")
  ll <- log_lik(ds)
  idx <- frmtmb.sample:::draws_par_index(f)
  r <- numeric(3L)
  for (k in 1:3) {
    sh <- frmtmb.sample:::draws_fit_at(ds, k, idx)
    nll <- as.numeric(frmtmb::build_objective(sh$frame)(sh$estimates))
    r[k] <- sum(ll[k, ]) - (-nll - re_prior_of(sh))
  }
  cat("  ", lab, ": dim ", paste(dim(ll), collapse = "x"),
      "  max|rowsum - objective| ", format(max(abs(r)), digits = 8),
      "  rel ", format(max(abs(r) / abs(rowSums(ll))), digits = 4),
      "  finite ", all(is.finite(ll)), "\n", sep = "")
}
run("ar(p = 3)", bf(y ~ x + ar(t, g, p = 3)), 3L)
run("arma(p = 3, q = 2)", bf(y ~ x + arma(t, g, p = 3, q = 2)), 5L)
run("ma(q = 3)", bf(y ~ x + ma(t, g, q = 3)), 3L)
run("ar(p = 1)", bf(y ~ x + ar(t, g, p = 1)), 1L)

cat("\n=== (e) where the OLD Student-t head becomes NaN, K = 2\n")
old_head <- function(nu, K = 2) lgamma((nu + K) / 2) - lgamma(nu / 2)
for (nu in c(1e305, 2e305, 3e305, 3.6e305, 4e305, 5e305, 5.1e305,
             5.2e305, 6e305)) {
  cat(sprintf("  nu %8.2e  old head %14s  log(nu*pi) %14s\n", nu,
              format(old_head(nu)), format(log(nu * pi))))
}
cat("\nDONE\n")
