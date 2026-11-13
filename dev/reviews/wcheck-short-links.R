# Short check of lane wcheck: the natural-scale prior rule adds one
# power for the log link's Jacobian. Other sigma links carry another
# Jacobian: identity none (sigma itself), squareplus about sigma^2
# (sigma ~ 1 / |eta| as eta -> -Inf). Verdict vs the profile slope per
# unit of -log sigma.
#   Rscript dev/reviews/wcheck-short-links.R <lib|base>
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r20"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (a[1] == "base") c(base, ul) else c(a[1], base, ul))
suppressMessages(library(frmtmb))
options(frmtmb.notices = FALSE)
cat("lib:", format(packageVersion("frmtmb")), find.package("frmtmb"), "\n")
q <- function(e) suppressWarnings(suppressMessages(e))
outcome <- function(expr) {
  r <- tryCatch(q(expr), error = function(e) e)
  if (inherits(r, "error")) {
    m <- gsub("\\s+", " ", conditionMessage(r))
    paste(if (grepl("has no maximum|is improper", m)) "REFUSED" else "ERROR",
          substr(m, 1, 120))
  } else sprintf("fit logLik %.6g betad %s conv %s", as.numeric(logLik(r)),
                 paste(format(r$estimates$betad, digits = 4), collapse = " "),
                 r$opt$convergence)
}
set.seed(20261111)
n <- 80
d <- data.frame(obs = factor(seq_len(n)), x = rnorm(n))
d$y <- 0.5 + 0.4 * d$x + rnorm(n)
fo <- bf(y | weights(w) ~ x + (1 | obs))
pr <- set_prior("normal(0, 1)", class = "sigma")
prof <- function(fam, eta) {
  ob <- q(frm(fo, family = fam, data = d, prior = pr, dry_run = "objective"))
  ob <- if (!is.null(ob$obj)) ob$obj else ob
  p <- ob$par
  k <- which(names(p) == "betad")[1L]
  free <- setdiff(seq_along(p), k)
  vapply(eta, function(t) {
    p[k] <- t
    fr <- function(r) {q2 <- p; q2[free] <- r
      v <- tryCatch(ob$fn(q2), error = function(e) NA_real_)
      if (!is.finite(v)) 1e10 else v}
    op <- nlminb(p[free], fr)
    if (op$objective < 1e9) p[free] <<- op$par
    -op$objective
  }, 0)
}
links <- list(
  log = list(eta = -c(2, 4, 6, 8, 10), sig = function(e) exp(e)),
  identity = list(eta = exp(-c(2, 4, 6, 8, 10)), sig = function(e) e),
  squareplus = list(eta = -exp(c(2, 4, 6, 8, 10)),
                    sig = function(e) (e + sqrt(e^2 + 4)) / 2))
for (E in c(0.5, 1.5)) {
  d$w <- 1 + E / n
  for (ln in names(links)) {
    fam <- brmsfamily("gaussian", link_sigma = ln)
    L <- links[[ln]]
    v <- tryCatch(prof(fam, L$eta), error = function(e) {
      cat("  profile ERROR", conditionMessage(e), "\n"); NULL})
    ls <- -log(L$sig(L$eta))
    cat(sprintf("E %g, link_sigma %s, normal(0, 1) on sigma | check: %s\n",
                E, ln, outcome(frm(fo, family = fam, data = d, prior = pr))))
    if (!is.null(v)) {
      cat("  -log sigma", format(ls, digits = 4), "\n  -fn",
          format(v, digits = 7), "\n  slope per unit -log sigma",
          format(diff(v) / diff(ls), digits = 3), "\n")
    }
  }
}
cat("DONE\n")
