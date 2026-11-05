# Prior tables of frmtmb against brms 2.23.0 for smooths and mm() terms:
# class, coef, group, resp, dpar, nlpar and lb, row by row, in order.
# Seed 43. Fit route in core; run with REFUSE_LIB=<lib> for the changed
# build and with REFUSE_LIB empty for rellib-r16.
lib <- Sys.getenv("REFUSE_LIB", "C:/Users/adf44/source/r/wt-refuse-lib")
.libPaths(c(if (nzchar(lib)) lib, "C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb)); attach(list(mmc = brms::mmc), name = "refuse_brms_mmc", warn.conflicts = FALSE)
cat("frmtmb", format(packageVersion("frmtmb")), "from",
    find.package("frmtmb"), "\n")
set.seed(43)
n <- 200
dd <- data.frame(x = runif(n), z = runif(n), w = rnorm(n),
                 f = factor(rep(c("a", "b"), n / 2)),
                 g = factor(rep(1:10, each = 20)))
dd$g1 <- factor(sample(letters[1:8], n, TRUE), levels = letters[1:8])
dd$g2 <- factor(sample(letters[1:8], n, TRUE), levels = letters[1:8])
dd$x1 <- rnorm(n)
dd$x2 <- rnorm(n)
dd$y <- sin(2 * pi * dd$x) + dd$z + rnorm(n, 0, 0.3)
dd$y2 <- dd$x + rnorm(n, 0, 0.3)
A <- 0.5^abs(outer(1:8, 1:8, "-"))
dimnames(A) <- list(letters[1:8], letters[1:8])
specs <- list(
  s = quote(bf(y ~ s(x))),
  s2 = quote(bf(y ~ s(x) + s(z, k = 5))),
  tight = quote(bf(y ~ s(z,k=5))),
  cr = quote(bf(y ~ s(z, bs = "cr", k = 5))),
  t2 = quote(bf(y ~ t2(x, z))),
  by = quote(bf(y ~ f + s(x, by = f))),
  sre = quote(bf(y ~ s(x) + (1 + w | g))),
  sig = quote(bf(y ~ s(x), sigma ~ s(z))),
  nl = quote(bf(y ~ a * exp(b * x), a ~ s(z), b ~ 1, nl = TRUE)),
  mv = quote(mvbf(bf(y ~ s(x)), bf(y2 ~ s(z)), rescor = FALSE)),
  mmid = quote(bf(y ~ x + (1 + mmc(x1, x2) | p | mm(g1, g2)),
                  sigma ~ 1 + (1 | p | mm(g1, g2)))),
  mmidarg = quote(bf(y ~ x + (1 | mm(g1, g2, id = "p")),
                     sigma ~ 1 + (1 | mm(g1, g2, id = "p")))),
  mmcov = quote(bf(y ~ x + (1 + mmc(x1, x2) | mm(g1, g2, cov = A)))),
  mmcov1 = quote(bf(y ~ x + (1 | mm(g1, g2, cov = A)))),
  mmidcov = quote(bf(y ~ x + (1 | p | mm(g1, g2, cov = A)),
                     sigma ~ 1 + (1 | p | mm(g1, g2, cov = A))))
)
cols <- c("class", "coef", "group", "resp", "dpar", "nlpar")
norm <- function(tab) {
  tab <- as.data.frame(tab)
  for (c in cols) tab[[c]] <- as.character(tab[[c]])
  out <- tab[, cols]
  out$lb <- ifelse(is.na(tab$lb) | tab$lb == "", "", as.character(tab$lb))
  out <- out[out$class != "theta", , drop = FALSE]
  rownames(out) <- NULL
  out
}
for (nm in names(specs)) {
  cat("\n=====", nm, ":", deparse1(specs[[nm]]), "\n")
  fp <- tryCatch(norm(frmtmb::default_prior(
    eval(specs[[nm]], asNamespace("frmtmb")), data = dd,
    data2 = list(A = A))), error = function(e) conditionMessage(e))
  sp <- specs[[nm]]
  if (identical(sp[[1L]], as.name("mvbf"))) {
    bform <- brms::bf(y ~ s(x)) + brms::bf(y2 ~ s(z)) + brms::set_rescor(FALSE)
  } else {
    sp[[1L]] <- quote(brms::bf)
    bform <- eval(sp)
  }
  bp <- tryCatch(norm(brms::default_prior(bform, data = dd,
                                          data2 = list(A = A))),
                 error = function(e) conditionMessage(e))
  if (is.character(fp) || is.character(bp)) {
    cat("frmtmb:", if (is.character(fp)) substr(fp, 1, 200) else "table", "\n")
    cat("brms:", if (is.character(bp)) substr(bp, 1, 200) else "table", "\n")
    next
  }
  same <- isTRUE(all.equal(fp, bp, check.attributes = FALSE))
  cat("identical rows:", same, "\n")
  if (!same) {
    key <- function(t) do.call(paste, c(t, sep = "|"))
    cat("only frmtmb:\n"); cat(setdiff(key(fp), key(bp)), sep = "\n")
    cat("only brms:\n"); cat(setdiff(key(bp), key(fp)), sep = "\n")
    if (setequal(key(fp), key(bp))) {
      cat("same rows, order differs\n")
      print(fp); print(bp)
    }
  }
}
