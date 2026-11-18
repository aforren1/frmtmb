# Review of lane sparsese: API probe on a small whitened REML fit.
.libPaths(c("C:/Users/adf44/source/r/wt-sparsese-lib",
            "C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("lib", dirname(find.package("frmtmb")), "\n")
set.seed(1)
d <- data.frame(x = stats::runif(200, 1, 2),
                g = factor(rep(1:20, length.out = 200)))
d$y <- sin(3 * d$x) + stats::rnorm(20, 0, 0.5)[d$g] + stats::rnorm(200, 0, 0.3)
for (k in 1:7) d[[paste0("x", k)]] <- d$x^k
fit <- suppressMessages(frm(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + (1 | g),
                            data = d, REML = TRUE))
obj <- fit$obj
cat("par_units:", format(fit$par_units), "\n")
cat("obj$par names:", names(obj$par), "\n")
cat("opt$par:", names(fit$opt$par), format(fit$opt$par), "\n")
cat("random len", length(obj$env$random), "\n")
cat("spHess formals:", names(formals(obj$env$spHess)), "\n")
cat("has he:", !is.null(obj$he), "\n")
print(args(frmtmb:::predict.frmtmb_fit))
print(names(formals(stats::predict)))
m <- methods(class = class(fit)[1])
print(m)
cat("whiten blocks:", length(fit$beta_whiten$blocks),
    "idx", fit$beta_whiten$blocks[[1]]$idx, "\n")
lp <- obj$env$last.par.best
cat("last.par.best names", unique(names(lp)), "\n")
H <- obj$env$spHess(lp, random = FALSE)
cat("H dim", dim(H), class(H), "\n")
