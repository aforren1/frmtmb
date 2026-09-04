## Elimination of the three hypotheses named in the frmtmb defect report:
##  (i)   static-pointer clobbering between two interleaved ODE tapes
##  (ii)  tmbstan evaluating through a different/retaped TMB object
##  (iii) R garbage collection invalidating the raw addresses held by the
##        C statics (F/X/Y in odesolve.c)
## plus the positive finding: deSolve "Returning early" hands back a matrix
## with fewer rows than `times`, which is what RTMB::ADjoint reports as
## "Wrong output length".
##
## Usage: Rscript 04-hypotheses.R <lib> <repro-dir>

args <- commandArgs(TRUE)
if (length(args) >= 1 && nzchar(args[1])) .libPaths(c(args[1], .libPaths()))
here <- if (length(args) >= 2) args[2] else "."
source(file.path(here, "lv-model.R"))

obs2 <- lv_data(2L)
obj2 <- lv_objective(obs2)
mode2 <- nlminb(obj2$par, obj2$fn, obj2$gr)$par

cat("== (i) two interleaved ODE tapes, direct path ==\n")
cat("fn:", obj2$fn(mode2), "\n")
g1 <- obj2$gr(mode2)
## Force many alternating evaluations; a clobbered static pointer would make
## repeated evaluations at the same point disagree.
set.seed(3)
agree <- TRUE
for (i in 1:20) {
    p <- mode2 + rnorm(length(mode2), sd = 0.05)
    a <- obj2$gr(p); b <- obj2$gr(p)
    if (!isTRUE(all.equal(as.numeric(a), as.numeric(b)))) agree <- FALSE
}
cat("repeated two-tape gradients reproducible:", agree, "\n")
cat("max|gr| at mode:", max(abs(g1)), "(two tapes, no error)\n")

cat("\n== (iii) gc pressure on the direct path ==\n")
obs1 <- lv_data(1L)
obj1 <- lv_objective(obs1)
mode1 <- nlminb(obj1$par, obj1$fn, obj1$gr)$par
ref <- obj1$gr(mode1)
gctorture2(step = 1)
tor <- tryCatch(obj1$gr(mode1), error = function(e) conditionMessage(e))
gctorture2(step = 0)
if (is.character(tor)) {
    cat("gctorture ERROR:", tor, "\n")
} else {
    cat("max abs diff vs untortured gradient:",
        max(abs(as.numeric(tor) - as.numeric(ref))), "\n")
}
gc()
cat("gradient after gc():", max(abs(as.numeric(obj1$gr(mode1)) - as.numeric(ref))), "\n")

cat("\n== the 'Wrong output length' mechanism ==\n")
## Trap the deSolve call and report the shape it actually returned.
shape <- NULL
trace(deSolve::ode, exit = quote({
    shape <<- c(shape, if (is.matrix(returnValue())) nrow(returnValue()) else NA)
}), print = FALSE, where = asNamespace("deSolve"))
bad <- mode1 + 100
res <- tryCatch(obj1$fn(bad), error = function(e) conditionMessage(e))
untrace(deSolve::ode, where = asNamespace("deSolve"))
cat("length(times) =", length(lv_times), "\n")
cat("nrow(sol) returned by deSolve:", paste(shape, collapse = ", "), "\n")
cat("objective result:", if (is.character(res)) paste("ERROR:", res) else res, "\n")
