# Short check of lane wcheck: a power-tail prior on Gamma's shape under
# links whose Jacobian has no power of the shape as it grows (identity,
# softplus, squareplus all go like eta). The lane's rule adds one power
# for every link. E = 12 (growth 6 per unit log shape), student_t(5.5)
# (p = 6.5): the lane says p - 1 < 6, refused; with no Jacobian power
# the objective falls like shape^(6 - 6.5).
#   Rscript dev/reviews/wcheck-short-shapelinks.R <lib|base>
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
mu <- exp(0.5 + 0.3 * d$x + rnorm(n, 0, 0.3))
d$yg <- rgamma(n, shape = 4, rate = 4 / mu)
d$w <- 1 + 12 / n
fo <- bf(yg | weights(w) ~ x + (1 | obs))
for (nu in c(5.5, 6.5)) {
  pr <- set_prior(sprintf("student_t(%g, 0, 2.5)", nu), class = "shape")
  for (ln in c("log", "identity", "softplus")) {
    fam <- brmsfamily("Gamma", link = "log", link_shape = ln)
    eta <- if (ln == "log") log(exp(c(4, 6, 8, 10, 12))) else exp(c(4, 6, 8, 10, 12))
    ob <- tryCatch(q(frm(fo, family = fam, data = d, prior = pr,
                         dry_run = "objective")), error = function(e) e)
    cat(sprintf("student_t(%g) on shape, link_shape %s | check: %s\n", nu, ln,
                outcome(frm(fo, family = fam, data = d, prior = pr))))
    if (inherits(ob, "error")) {cat("  ERROR", conditionMessage(ob), "\n"); next}
    ob <- if (!is.null(ob$obj)) ob$obj else ob
    p <- ob$par
    k <- which(names(p) == "betad")[1L]
    free <- setdiff(seq_along(p), k)
    v <- vapply(eta, function(t) {
      p[k] <- t
      fr <- function(r) {q2 <- p; q2[free] <- r
        vv <- tryCatch(ob$fn(q2), error = function(e) NA_real_)
        if (!is.finite(vv)) 1e10 else vv}
      op <- nlminb(p[free], fr)
      if (op$objective < 1e9) p[free] <<- op$par
      -op$objective
    }, 0)
    lk <- c(4, 6, 8, 10, 12)
    cat("  log shape", lk, "\n  -fn", format(v, digits = 7),
        "\n  slope per unit log shape", format(diff(v) / diff(lk), digits = 3),
        "\n")
  }
}
cat("DONE\n")
