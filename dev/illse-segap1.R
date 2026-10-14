# Lane illse, item b: one silent fit of dev/illse-segap.R (ls ~ mo(inc),
# income category 2 never observed). The estimates, the simplex, the
# unit-diagonal spectrum of the AD Hessian with its eigenvectors, and
# the objective along the smallest direction.
#   Rscript dev/illse-segap1.R base|dev <seed> [formula 1|2|4]
args <- commandArgs(TRUE)
arm <- args[1]
seed <- as.integer(args[2])
which_fo <- if (length(args) > 2) as.integer(args[3]) else 1L
R <- "C:/Users/adf44/source/r"
.libPaths(c(file.path(R, "rellib-r8"),
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
if (arm == "dev") {
  suppressMessages(pkgload::load_all(file.path(R, "frmtmb-wt-illse"),
                                     quiet = TRUE))
} else suppressMessages(library(frmtmb))
se_mo_gap <- function(seed) {
  set.seed(seed)
  inc <- sample(c(0L, 1L, 3L), 100, TRUE)
  ls <- c(30, 60, 70, 75)[inc + 1L] + rnorm(100, sd = 7)
  data.frame(inc, ls, age = rnorm(100, mean = 40, sd = 10))
}
src <- readLines("C:/Users/adf44/source/r/frmtmb-wt-illse/dev/relrev070-interplay.R")
i <- grep("^mk <- function", src)
j <- grep("^designs <- list", src)
eval(parse(text = src[i:(j - 1L)]))
fo <- switch(which_fo, ls ~ mo(inc), ls ~ mo(inc) + age, ls ~ mo(inc) * age,
             y ~ mo(x2))
d <- if (which_fo == 4L) mk(seed) else se_mo_gap(seed)
w <- character()
f <- withCallingHandlers(frm(fo, data = d), warning = function(x) {
  w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
print(table(if (which_fo == 4L) d$x2 else d$inc))
cat("warnings:", length(w), "\n"); if (length(w)) cat(substr(w, 1, 200), sep = "\n")
p <- f$opt$par
print(p)
ns <- asNamespace("frmtmb")
cat("lost:", names(ns$sdr_of(f)$se_lost), "\n")
print(sqrt(diag(vcov(f, full = TRUE))))
H <- f$obj$he(p)
D <- sqrt(diag(H))
e <- eigen(H / outer(D, D), symmetric = TRUE)
cat("unit-diagonal spectrum / largest:\n"); print(e$values / max(e$values))
cat("smallest direction (natural units):\n")
v <- e$vectors[, length(p)] / D
print(setNames(v / max(abs(v)), names(p)))
f0 <- f$obj$fn(p)
for (s in c(1e-3, 1e-2, 0.1, 1, 3)) {
  cat(sprintf("step %g along it: objective change %.3e, %.3e\n", s,
              f$obj$fn(p + s * v / sqrt(sum(v^2))) - f0,
              f$obj$fn(p - s * v / sqrt(sum(v^2))) - f0))
}
