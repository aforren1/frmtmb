# Lane wcheck2, the asym_laplace refusal: what the frame records for each
# term the Laplace route integrates (block fields), and frm_compat()'s
# rows for asym_laplace against the structures and specials.
#   Rscript dev/wcheck2-ald-blocks.R <base|lane|lib>
a <- commandArgs(TRUE)
lib_arg <- a[1]
source("dev/wcheck2-common.R")
`%||%` <- function(a, b) if (is.null(a)) b else a
set.seed(1)
n <- 80
d <- data.frame(x = rnorm(n), z = runif(n), g = factor(rep(1:20, each = 4)),
                tt = factor(rep(1:4, 20)))
d$y <- 1 + d$x + rnorm(n)
fos <- list(
  "(1 | g)" = bf(y ~ x + (1 | g)),
  "(1 + x | g)" = bf(y ~ x + (1 + x | g)),
  "(1 + x || g)" = bf(y ~ x + (1 + x || g)),
  "s(z)" = bf(y ~ x + s(z)),
  "gp(z)" = bf(y ~ x + gp(z)),
  "gp(z, k = 5)" = bf(y ~ x + gp(z, k = 5, c = 5 / 4)),
  "ar1(tt + 0 | g)" = bf(y ~ x + ar1(tt + 0 | g)),
  "sigma ~ (1 | g)" = bf(y ~ x, sigma ~ (1 | g)),
  "(1 | p | g) mu and sigma" = bf(y ~ x + (1 | p | g), sigma ~ (1 | p | g)))
for (nm in names(fos)) {
  fr <- tryCatch(q(frm(fos[[nm]], family = asym_laplace(), data = d,
                       dry_run = "frame")), error = function(e) e)
  if (inherits(fr, "error")) {
    cat(nm, ": ERROR", conditionMessage(fr), "\n")
    next
  }
  cat(nm, "\n")
  for (bk in fr$re_blocks) {
    cps <- vapply(bk$components %||% list(), function(cp) cp$lp_key %||% "",
                  "")
    cat(sprintf("  covstruct %s, dpar %s, lp %s, theta %s, ncp %s, label %s, fields %s\n",
                bk$covstruct %||% "?", bk$dpar %||% "?",
                paste(cps, collapse = ","),
                paste(bk$theta_idx, collapse = ","),
                format(isTRUE(bk$ncp)), bk$term_label %||% "",
                paste(names(bk), collapse = " ")))
  }
}

cm <- frm_compat("asym_laplace")
cat("\nfrm_compat(\"asym_laplace\"):", nrow(cm), "rows\n")
print(cm[cm$status != "works" | grepl("cov|special", cm$kind_b %||% ""),
         intersect(names(cm), c("a", "b", "status", "note"))][1:5, ])
cat(names(cm), "\n")
DONE <- TRUE
cat("DONE\n")
