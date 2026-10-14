# Lane optstart: the stages of the wide design's fit (dev/optrap-p2-wide1.R,
# n 101000, 50 columns) with verbose = 2, so that each optimizer run,
# restart and Newton step is seen with its objective and code.
#   Rscript dev/optstart-wide-trace2.R <lib or ""> [seed]
a <- commandArgs(TRUE)
.libPaths(c(if (nzchar(a[1])) a[1], "C:/Users/adf44/source/r/rellib-r8",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
s <- if (length(a) > 1) as.integer(a[2]) else 1L
cat("lib", dirname(find.package("frmtmb")), "seed", s, "\n")
set.seed(s)
n <- 101000
p <- 49
Z <- matrix(rnorm(n * p), n, p)
colnames(Z) <- paste0("z", seq_len(p))
d <- data.frame(Z, g = factor(rep(seq_len(1000), length.out = n)))
d$y <- 1 + 1e4 * d$z1 + rnorm(1000, 0, 0.5)[d$g] + rnorm(n)
fo <- as.formula(paste("y ~", paste(colnames(Z), collapse = " + "),
                       "+ (1 | g)"))
f <- withCallingHandlers(frm(fo, data = d,
                             control = frmtmb_control(verbose = 2)),
  warning = function(w) {
    cat("WARNING:", substr(conditionMessage(w), 1, 100), "\n")
    invokeRestart("muffleWarning")
  })
cat("code", f$opt$convergence, f$opt$message, "| iterations",
    f$opt$iterations, "| objective", format(f$opt$objective, digits = 15),
    "\n")
g <- drop(f$obj$gr(f$opt$par))
u <- f$par_units %||% rep(1, length(g))
o <- order(-abs(g * u))[1:6]
cat("largest |grad * unit|:\n")
print(data.frame(par = names(f$opt$par)[o], value = f$opt$par[o],
                 grad = g[o], unit = u[o]))
cat("sigma", exp(f$opt$par[names(f$opt$par) == "betad"]), "group sd",
    exp(f$opt$par[names(f$opt$par) == "theta"]), "\n")
