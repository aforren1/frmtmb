# Reviewer check: the both-converged median |z| of v067 and r8 over the
# SAME models, from the lane's plausibility.rds files (same brms
# reference). The lane's headline compares two different sets.
root <- "C:/Users/adf44/source/r/frmtmb-wt-remeasure/dev/remeasure-out"
rd <- function(t) readRDS(file.path(root, t, "port", "plausibility.rds"))
ok2 <- function(P) {
  b <- P$table
  b$id[b$frmtmb == "OK" & b$brms == "OK" & !is.na(b$rhat) &
         b$rhat <= 1.05 & !is.na(b$fconv) & b$fconv == 0]
}
z <- function(P, ids) unlist(lapply(P$detail[ids], function(d) abs(d$fx$z)))
A <- rd("v067"); B <- rd("r8")
ia <- ok2(A); ib <- ok2(B)
cat("v067 set", length(ia), "models; r8 set", length(ib), "models\n")
cat("only v067:", setdiff(ia, ib), "\nonly r8:", setdiff(ib, ia), "\n")
f <- function(x) sprintf("n %d median %.4f", length(x), median(x))
cat("as published  v067:", f(z(A, ia)), " r8:", f(z(B, ib)), "\n")
cm <- intersect(ia, ib)
cat("common set    v067:", f(z(A, cm)), " r8:", f(z(B, cm)), "\n")
u <- union(ia, ib)
u <- u[!vapply(A$detail[u], is.null, NA) & !vapply(B$detail[u], is.null, NA)]
cat("union, every fit regardless of code  v067:", f(z(A, u)), " r8:",
    f(z(B, u)), "\n")
cat("identical z on the common set:",
    isTRUE(all.equal(z(A, cm), z(B, cm))), "\n")
