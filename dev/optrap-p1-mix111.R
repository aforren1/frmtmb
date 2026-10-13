# Lane optrap, punch round 1 (m3): test-v18's group-level mixture with a
# random intercept, which a unit on the components' mu ran 111 units
# lower. Prints the start, the units and the optimum per build.
#   Rscript dev/optrap-p1-mix111.R <lib>
a <- commandArgs(TRUE)
.libPaths(c(a[1], "C:/Users/adf44/source/r/rellib-r7",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
set.seed(19)
ng <- 40
m <- 6
cls <- rbinom(ng, 1, 0.4)
g <- rep(seq_len(ng), each = m)
y <- rnorm(ng * m, ifelse(cls == 1, 2, -1)[g], 0.8)
dd <- data.frame(y = y, g = factor(g))
w <- character()
fre <- withCallingHandlers(frm(bf(y ~ 1 + (1 | g)) +
  mixture(gaussian(), gaussian(), groups = ~g), data = dd),
  warning = function(x) {
    w <<- c(w, conditionMessage(x))
    invokeRestart("muffleWarning")
  })
cat("lib", dirname(find.package("frmtmb")), "\n")
u <- tryCatch(frmtmb:::predictor_units(fre$frame), error = function(e) NULL)
cat("units:", deparse(u), " par_units:", deparse(fre$par_units), "\n")
cat("objective", fre$opt$objective, "code", fre$opt$convergence,
    fre$opt$message, "\n")
print(fre$opt$par, digits = 5)
cat("start theta:", deparse(fre$frame$par_template$theta), "\n")
cat("warnings:", substr(w, 1, 90), sep = "\n  ")
b <- fre$obj$env$last.par.best
b <- b[names(b) == "b"]
cat("modes: range", range(b), " by class:",
    tapply(b[seq_len(ng)], cls, mean), "\n")
