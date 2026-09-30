# Reviewer of lane defects: frmtmb's default_prior() table against brms
# 2.23.0's on a set of models: the rows (class, coef, group, resp, dpar,
# nlpar) and their ORDER. Run once per arm:
#   Rscript dev/defects-rev-priors.R lane|base
arm <- commandArgs(trailingOnly = TRUE)[1]
libs <- c("C:/Users/adf44/source/r/rellib-r3",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
if (identical(arm, "lane")) libs <- c("C:/Users/adf44/source/r/wt-defects-lib", libs)
.libPaths(libs)
suppressMessages({loadNamespace("frmtmb"); loadNamespace("brms")})
set.seed(20260929)
n <- 60
d <- data.frame(y = rnorm(n), x = rnorm(n), z = rnorm(n), w = rnorm(n),
                g = factor(rep(1:6, 10)), h = factor(rep(1:5, 12)),
                f = factor(sample(c("a", "b", "c"), n, TRUE)),
                cnt = rpois(n, 3), yb = rbinom(n, 1, 0.5),
                yo = factor(sample(1:4, n, TRUE), ordered = TRUE),
                y2 = rnorm(n), t = 1:n, yp = runif(n, 0.05, 0.95))
models <- list(
  glm_order = quote(list(y ~ z + x + w, gaussian())),
  factor = quote(list(y ~ f + x, gaussian())),
  smooth = quote(list(y ~ z + s(x), gaussian())),
  smooth2 = quote(list(y ~ s(w) + s(x) + z, gaussian())),
  smooth_by = quote(list(y ~ s(x, by = f) + f, gaussian())),
  re = quote(list(y ~ x + (1 + x | g), gaussian())),
  re2 = quote(list(y ~ x + (1 | g) + (1 + z | h), gaussian())),
  dist = quote(list(bf(y ~ x, sigma ~ z + x), gaussian())),
  noint = quote(list(y ~ 0 + x + z, gaussian())),
  noint_f = quote(list(y ~ 0 + f + x, gaussian())),
  pois = quote(list(cnt ~ x + z, poisson())),
  ord = quote(list(yo ~ x + z, cumulative())),
  ord_cs = quote(list(yo ~ x + cs(z), sratio())),
  student = quote(list(y ~ x, student())),
  nl = quote(list(bf(y ~ b1 * exp(b2 * x), b1 ~ 1, b2 ~ z, nl = TRUE),
                  gaussian())),
  nl2 = quote(list(bf(y ~ b1 + b2 * x, b1 ~ 1 + z, b2 ~ 1 + w, nl = TRUE),
                   gaussian())),
  nl_smooth = quote(list(bf(y ~ lp, lp ~ z + s(x) + (1 | g), nl = TRUE),
                         gaussian())),
  mv = quote(list(bf(y ~ x + z) + bf(y2 ~ w) + set_rescor(FALSE),
                  gaussian())),
  zi = quote(list(bf(cnt ~ x, zi ~ z), zero_inflated_poisson())),
  mo = quote(list(y ~ mo(yo) + x, gaussian())),
  beta = quote(list(yp ~ x, Beta())),
  beta_bad = quote(list(y ~ x, Beta())),
  bern_bad = quote(list(cnt ~ x, bernoulli())),
  pois_bad = quote(list(y ~ x, poisson())),
  cat_num = quote(list(cnt ~ x, categorical()))
)
key <- function(t) {
  t <- as.data.frame(t)
  cols <- c("class", "coef", "group", "resp", "dpar", "nlpar")
  for (k in cols) if (is.null(t[[k]])) t[[k]] <- ""
  apply(t[cols], 1, function(r) paste(ifelse(is.na(r), "", r),
                                      collapse = "|"))
}
env_for <- function(pkg) {
  e <- new.env(parent = asNamespace(pkg))
  e$gaussian <- stats::gaussian; e$poisson <- stats::poisson
  e
}
for (nm in names(models)) {
  mb <- eval(models[[nm]], env_for("brms"))
  mf <- eval(models[[nm]], env_for("frmtmb"))
  b <- tryCatch(key(suppressWarnings(suppressMessages(
    brms::default_prior(mb[[1]], data = d, family = mb[[2]])))),
    error = function(e) paste("BRMS ERROR:", conditionMessage(e)))
  f <- tryCatch(key(suppressWarnings(suppressMessages(
    frmtmb::default_prior(mf[[1]], data = d, family = mf[[2]])))),
    error = function(e) paste("FRMTMB ERROR:", substr(conditionMessage(e),
                                                       1, 200)))
  common_b <- b[b %in% f]
  common_f <- f[f %in% b]
  cat(sprintf("== %s: brms %d rows, frmtmb %d rows, common %d, order of common %s\n",
              nm, length(b), length(f), length(common_b),
              if (identical(common_b, common_f)) "SAME" else "DIFFERENT"))
  # the order within each (resp, dpar, nlpar, class, group) block
  blk <- function(k) sub("[|][^|]*[|]([^|]*[|][^|]*[|][^|]*[|][^|]*)$",
                         "|\\1", k)
  wb <- split(common_b, blk(common_b)); wf <- split(common_f, blk(common_f))
  wdiff <- names(wb)[!vapply(names(wb), function(k) identical(wb[[k]], wf[[k]]), NA)]
  cat("   within-block order differs in:", if (length(wdiff)) paste(wdiff, collapse = ", ") else "none", "\n")
  if (!identical(common_b, common_f)) {
    cat("  brms  :", paste(common_b, collapse = "  ;  "), "\n")
    cat("  frmtmb:", paste(common_f, collapse = "  ;  "), "\n")
  }
  ob <- setdiff(b, f); of <- setdiff(f, b)
  if (length(ob)) cat("  only brms  :", paste(ob, collapse = "  ;  "), "\n")
  if (length(of)) cat("  only frmtmb:", paste(of, collapse = "  ;  "), "\n")
}
