# Does the response's unit in imp_outer_steps() matter? no_units passes every test.
# Same data on the lane's build and on the mutant: y in large units, identity link.
source("C:/Users/adf44/source/r/frmtmb-wt-impsamp/dev/review-impsamp-common.R")
set_libs(Sys.getenv("REV_ARM", "after"))
suppressPackageStartupMessages(library(frmtmb))
cat(find.package("frmtmb"), "\n")
ses <- function(f) tryCatch(sqrt(diag(frmtmb:::sdr_of(f)$cov.fixed)), error = function(e) NA)
set.seed(1)
G <- 15L
g <- factor(rep(seq_len(G), each = 10))
z <- rnorm(150)
u <- rnorm(G, 0, 0.7)[g]
for (fam in c("student", "gaussian", "skew")) for (s in c(1e2, 1e4, 1e5)) {
  yy <- (0.3 + 0.5 * z + u + 0.5 * rt(150, 5)) * s
  d <- data.frame(y = yy, z, g)
  fo <- if (fam == "student") bf(y ~ z + (1 | g), nu = 4) else bf(y ~ z + (1 | g))
  fy <- switch(fam, student = student(), gaussian = gaussian(), skew = exgaussian())
  w <- character(0)
  a <- withCallingHandlers(tryCatch(frm(fo, family = fy, data = d, importance = 300L), error = function(e) e),
                           warning = function(x) {
                             w <<- c(w, substr(conditionMessage(x), 1, 60))
                             invokeRestart("muffleWarning")
                           })
  d1 <- transform(d, y = y / s)
  r <- suppressWarnings(tryCatch(frm(fo, family = fy, data = d1, importance = 300L), error = function(e) e))
  if (inherits(a, "error") || inherits(r, "error")) {
    cat(sprintf("%s y*%g STOP: %s\n", fam, s, substr(conditionMessage(if (inherits(a, "error")) a else r), 1, 130)))
    next
  }
  sa <- ses(a); sr <- ses(r)
  cat(sprintf("%s y*%g rounds %d conv %d grad %.3g slope z %+.4f SE ratio %s ll-unit %+.5f mcse %.4f warns %d %s\n", fam, s, a$importance$rounds, a$opt$convergence, a$importance$grad,
              (a$opt$par[[2]] / s - r$opt$par[[2]]) / sr[[2]],
              format(sa[[2]] / s / sr[[2]], digits = 5), as.numeric(logLik(a)) - (as.numeric(logLik(r)) - 150 * log(s)),
              a$importance$mcse, length(w), paste(w, collapse = ";")))
}
