# Insert the bound-aware-covariance defect into dev/test-backlog.md,
# immediately above the "## Reference" heading, so the file is edited by a
# script rather than by a shell quoting dance.

p <- "dev/test-backlog.md"
old <- readLines(p, encoding = "UTF-8", warn = FALSE)
if (any(grepl("bound-aware", old, fixed = TRUE))) {
  stop("the entry is already there")
}
at <- grep("^## Reference", old)
stopifnot(length(at) == 1L)
new <- c(
"- **The covariance machinery is not bound-aware.** At a constrained",
"  optimum the UNCONSTRAINED Hessian can be indefinite, because the fit",
"  is only a minimum along the feasible directions. `sdreport()` reads",
"  the full Hessian, so such a fit reports `pdHess = FALSE` and NaN",
"  standard errors while being exactly right. Constructed",
"  (`dev/gradcheck-01-construct.R`, seed 101,",
"  `dev/gradcheck-rev-10-flip.R` for the cap sweep): gaussian",
"  `y ~ x` with a true slope of 2 under",
"  `set_prior(\"\", class = \"b\", ub = 0.1)` stops with `x` on the bound;",
"  the full Hessian's eigenvalues are 693.40, 73.10 and -29.63, and the",
"  same Hessian restricted to the two parameters no bound holds is",
"  positive definite, which is why the convergence check can measure a",
"  headroom there. Whether it bites depends on how far the bound is from",
"  the unconstrained optimum, not on bound-awareness: over caps 0.1 to",
"  1.99 on that design the bound holds `x` at every cap while `pdHess`",
"  is FALSE at 0.1 and 0.5 and TRUE at 1 and above. Fix: restrict the",
"  reported covariance to the free subspace, and say in the report that",
"  a bound-held parameter has no standard error rather than returning",
"  NaN for every parameter. Test: on the `ub = 0.1` fit, the free-set",
"  Hessian is positive definite while `pdHess` is FALSE, and the free",
"  parameters get finite standard errors.",
"")
writeLines(c(old[seq_len(at - 1L)], new, old[at:length(old)]), p,
           useBytes = TRUE)
cat("inserted", length(new), "lines above line", at, "\n")
