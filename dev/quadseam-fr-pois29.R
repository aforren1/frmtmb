# Fix round: the review's poissd2 seed 29 (gen() of
# dev/reviews/quadseam-review-acc.R) at the default grad_tol. What the
# gradient verdict says, and whether the final tape's AD Hessian is
# finite; then the verbose stages.
#   Rscript dev/quadseam-fr-pois29.R <lib|base>
a <- commandArgs(TRUE)
base <- "C:/Users/adf44/source/r/rellib-r20"
ul <- "C:/Users/adf44/AppData/Local/R/win-library/4.6"
.libPaths(if (a[1] == "base") c(base, ul) else c(a[1], base, ul))
suppressMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), "\n")
set.seed(29)
d <- data.frame(g = factor(rep(1:50, each = 3)), x = rnorm(150))
u <- rnorm(50, 0, 2)
d$y <- rpois(nrow(d), exp(-1 + 0.3 * d$x + u[d$g]))
w <- list()
f <- withCallingHandlers(
  frm(y ~ x + (1 | g), family = poisson(), data = d, quadrature = TRUE),
  warning = function(x) {
    w[[length(w) + 1L]] <<- x; invokeRestart("muffleWarning")
  })
cat(sprintf("logLik %.10f conv %d grad %.3e settle %s\n",
            as.numeric(logLik(f)), f$opt$convergence,
            max(abs(f$obj$gr(f$opt$par))),
            paste(format(f$opt$quad_settle %||% NA, digits = 3),
                  collapse = "/")))
for (x in w) cat("WARN", class(x)[1], ":",
                 gsub("[[:space:]]+", " ", conditionMessage(x)), "\n")
v <- frmtmb:::grad_verdict(f, drop(f$obj$gr(f$opt$par)), f$control)
cat(sprintf("verdict gmax %.3e proj %.3e headroom %s warn %s\n", v$gmax,
            v$proj, format(v$headroom, digits = 3), v$warn))
H <- f$obj$he(f$opt$par)
cat("AD hessian finite:", all(is.finite(H)), "\n")
print(H)
se <- tryCatch(sqrt(diag(frmtmb:::sdr_of(f)$cov.fixed)),
               error = function(e) conditionMessage(e))
cat("SE:", format(se, digits = 6), "\n")
invisible(suppressWarnings(frm(y ~ x + (1 | g), family = poisson(), data = d,
                               quadrature = TRUE, verbose = TRUE)))
