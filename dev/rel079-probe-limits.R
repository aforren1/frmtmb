# Finding 13 and finding 1's limit with real rows, not moved options: a latent residual
# at 3000, 3001 and 20001 rows, alone and beside (1 | s). Each route says what it does:
# values, NA with a message, or a refusal. One process per arm, for the peak memory.
#   Rscript dev/rel079-probe-limits.R <rows> <alone|block>
a <- commandArgs(TRUE)
n <- as.integer(a[1]); arm <- a[2]
.libPaths(c("C:/Users/adf44/source/r/rellib-r17", "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), as.character(packageVersion("frmtmb")),
    format(file.mtime(system.file("DESCRIPTION", package = "frmtmb"))), "\n")
set.seed(7500 + n)
nt <- 50L
ng <- ceiling(n / nt)
d <- expand.grid(t = seq_len(nt), s = factor(seq_len(ng)))[seq_len(n), ]
d$x <- rnorm(n)
e <- ave(rnorm(n), d$s, FUN = function(z) stats::filter(z, 0.5, "recursive"))
d$y <- rpois(n, exp(1 + 0.3 * d$x + rnorm(ng, 0, 0.4)[d$s] + 0.4 * e))
form <- if (arm == "block") {
  y ~ x + (1 | s) + ar(time = t, gr = s, cov = TRUE)
} else {
  y ~ x + ar(time = t, gr = s, cov = TRUE)
}
t0 <- proc.time()[[3]]
f <- suppressWarnings(suppressMessages(frm(form, data = d, family = poisson())))
cat(sprintf("n %d %s fit %.1f s conv %d\n", n, arm, proc.time()[[3]] - t0, f$opt$convergence))
nd <- data.frame(t = nt + 1L, s = factor(c(1, 2), levels = levels(d$s)), x = c(0, 1))
route <- function(lab, expr) {
  said <- character()
  t0 <- proc.time()[[3]]
  gc(reset = TRUE, full = FALSE)
  v <- tryCatch(withCallingHandlers(expr,
      message = function(m) { said <<- c(said, class(m)[1L]); invokeRestart("muffleMessage") },
      warning = function(w) invokeRestart("muffleWarning")),
    error = function(e) structure(list(), err = conditionMessage(e), cls = class(e)[1L]))
  tt <- proc.time()[[3]] - t0
  peak <- sum(gc()[, 6L])
  if (!is.null(attr(v, "err"))) {
    cat(sprintf("%-24s REFUSED in %.1f s, %s: %s\n", lab, tt, attr(v, "cls"), substr(attr(v, "err"), 1, 110)))
    return(invisible())
  }
  se <- if (is.list(v) && !is.null(v$se.fit)) v$se.fit else if (is.matrix(v) && "Est.Error" %in% colnames(v)) v[, "Est.Error"] else NA
  cat(sprintf("%-24s %.1f s, peak %.0f MB, se finite %d of %d, first %s, messages: %s\n", lab, tt, peak,
              sum(is.finite(se)), length(se), format(se[1L], digits = 3), paste(unique(said), collapse = ",")))
}
route("fitted()", fitted(f))
route("residuals()", residuals(f))
route("fitted(new)", fitted(f, newdata = nd))
route("fitted(new, linear)", fitted(f, newdata = nd, scale = "linear"))
route("frm_linpred(se)", frm_linpred(f, se.fit = TRUE))
route("frm_linpred(new, se)", frm_linpred(f, newdata = nd, se.fit = TRUE))
route("conditional_effects", {
  ce <- conditional_effects(f, "x")
  matrix(ce[[1L]]$se__, ncol = 1L, dimnames = list(NULL, "Est.Error"))
})
cat("LIMITS END\n")
