## Reviewer, claim 1 and claim 8: does anything catch a MISLABEL? Apply
## one mutation in memory (assignInNamespace, nothing installed), then
## run the lane's pinning test files and an independent semantic check.
##   Rscript dev/brmsnames-rev-mutants.R none|rlev|rcoef|corord <testfile>...
## Mutations:
##   rlev   brms_r_labels() writes the levels in reverse order
##   rcoef  brms_r_labels() writes the coefficient names in reverse order
##   corord draws_cor_index() returns R's column-major lower triangle
## Semantic check: data seed 5, sampler seed 8, y ~ x + (1 + x + z + w | g).
av <- commandArgs(trailingOnly = TRUE)
mut <- av[1]; files <- av[-1]
.libPaths(c("C:/Users/adf44/source/r/brmsnames-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true",
           FRMTMB_STAN_CACHE = normalizePath("dev/stan-cache"))
q <- function(e) suppressWarnings(suppressMessages(e))
q(library(testthat)); q(library(frmtmb)); q(library(frmtmb.sample))
ns <- asNamespace("frmtmb"); nss <- asNamespace("frmtmb.sample")
if (mut %in% c("rlev", "rcoef")) {
  f <- get("brms_r_labels", ns)
  src <- deparse(f)
  if (mut == "rlev") {
    i <- grep("lev <- bk[[\"levels\"]]", src, fixed = TRUE)
    src[i] <- sub("lev <- bk[[\"levels\"]]", "lev <- rev(bk[[\"levels\"]])",
                  src[i], fixed = TRUE)
  } else {
    i <- grep("nm <- paste0(\"r_\"", src, fixed = TRUE)
    src <- append(src, "    cf <- rev(cf)", after = i - 1L)
  }
  stopifnot(length(i) == 1L)
  g <- eval(parse(text = src)); environment(g) <- ns
  assignInNamespace("brms_r_labels", g, "frmtmb")
}
if (mut == "corord") {
  g <- function(K) which(lower.tri(diag(K)))
  environment(g) <- nss
  assignInNamespace("draws_cor_index", g, "frmtmb.sample")
}
cat("mutation:", mut, "\n")

# independent semantic check
set.seed(5)
G <- 10; n <- 300
d <- data.frame(x = rnorm(n), z = rnorm(n), w = rnorm(n),
                g = factor(rep(1:G, 30)))
S <- matrix(c(1, .5, -.3, .2, .5, 1, .3, -.2, -.3, .3, 1, .1, .2, -.2, .1, 1), 4)
L <- t(chol(S * 0.4))
U <- t(L %*% matrix(rnorm(4 * G), 4))
d$y <- 1 + d$x + U[d$g, 1] + U[d$g, 2] * d$x + U[d$g, 3] * d$z +
  U[d$g, 4] * d$w + rnorm(n)
fit <- q(frm(bf(y ~ x + (1 + x + z + w | g)), family = gaussian(), data = d))
ds <- q(frm_sample(fit, chains = 1, iter = 100, refresh = 0, seed = 8))
idx <- frmtmb.sample:::draws_par_index(fit)
vc <- VarCorr(ds, summary = FALSE)$g$cor
re <- ranef(ds, summary = FALSE)$g
bad_r <- 0; bad_re <- 0; bad_c <- 0; nr <- 0; nc <- 0
for (i in c(1, 25, 50)) {
  sh <- frmtmb.sample:::draws_fit_at(ds, i, idx)
  M <- unclass(ranef(sh))[[1]]
  C <- stats::cov2cor(varcorr_matrices(sh)[[1]])
  cn <- gsub("[()]", "", colnames(M))
  for (lv in rownames(M)) for (k in seq_along(cn)) {
    nr <- nr + 1
    col <- paste0("r_g[", lv, ",", cn[k], "]")
    if (!identical(unname(ds$draws[i, col]), unname(M[lv, k]))) bad_r <- bad_r + 1
    if (!identical(unname(re[i, lv, cn[k]]), unname(M[lv, k]))) bad_re <- bad_re + 1
  }
  for (a in 1:4) for (b in 1:4) if (a != b) {
    nc <- nc + 1
    if (abs(vc[i, cn[a], cn[b]] - C[a, b]) > 1e-12) bad_c <- bad_c + 1
  }
}
cat(sprintf("semantic: r_ columns %d checked, %d wrong; ranef(ds) cells %d wrong; cor cells %d checked, %d wrong\n",
            nr, bad_r, bad_re, nc, bad_c))

for (file in files) {
  pkg <- if (grepl("frmtmb.sample", file)) "frmtmb.sample" else "frmtmb"
  res <- q(test_file(file, reporter = "silent", package = pkg))
  df <- as.data.frame(res)
  cat(sprintf("%s: blocks %d PASS %d FAIL %d ERROR %d SKIP %d\n", basename(file),
              nrow(df), sum(df$passed), sum(df$failed), sum(df$error), sum(df$skipped)))
}
cat("DONE\n")
