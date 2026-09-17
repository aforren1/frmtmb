## brms's own answers for the surfaces lane brmsnames changes, read off
## brms 2.23.0's bundled example fits so no Stan program is compiled.
##   Rscript dev/brmsnames-brms-ref.R
.libPaths(c("C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(brms))
cat("brms", as.character(packageVersion("brms")), "\n\n")
f <- brms:::brmsfit_example1
f <- brms:::restructure(f)
print(f$formula)
cat("\n== variables ==\n"); print(variables(f))
m <- function(nm) {
  g <- getS3method(nm, "brmsfit")
  cat(sprintf("%-20s (%s)\n", nm, paste(names(formals(g)), collapse = ", ")))
}
cat("\n== formals ==\n")
for (nm in c("fixef", "ranef", "coef", "VarCorr", "posterior_summary",
             "bayes_R2", "hypothesis", "summary", "as.matrix", "as.array",
             "as.data.frame", "as_draws_array", "as_draws_df",
             "as_draws_matrix", "as_draws_list", "as_draws_rvars",
             "as_draws", "variables", "loo_R2", "conditional_effects",
             "pairs", "pp_check")) m(nm)

sh <- function(x) {
  if (is.list(x) && !is.data.frame(x)) {
    return(paste0("list(", paste(names(x), collapse = ","), ")"))
  }
  d <- dim(x)
  paste(class(x)[1], if (is.null(d)) length(x) else paste(d, collapse = "x"))
}
cat("\n== fixef ==\n")
print(fixef(f)); str(fixef(f, summary = FALSE))
print(dimnames(fixef(f, summary = FALSE, robust = TRUE))[[2]])
print(fixef(f, probs = c(0.1, 0.9)))
print(fixef(f, robust = TRUE)[1, ])
cat("\n== ranef ==\n")
r <- ranef(f); cat(sh(r), "\n"); str(r)
rs <- ranef(f, summary = FALSE); cat(sh(rs), "\n"); str(rs)
cat("\n== coef ==\n")
cc <- coef(f); str(cc)
cs <- coef(f, summary = FALSE); str(cs)
cat("\n== VarCorr ==\n")
v <- VarCorr(f); str(v)
print(v)
vs <- VarCorr(f, summary = FALSE); str(vs)
cat("\n== posterior_summary ==\n")
ps <- posterior_summary(f); print(head(ps)); print(dim(ps))
print(posterior_summary(f, variable = "b_Intercept"))
print(posterior_summary(f, pars = "^b_"))
cat("\n== hypothesis ==\n")
h <- hypothesis(f, "Trt1 > 0")
print(class(h)); print(names(h)); str(h$hypothesis); print(h)
str(h$samples); print(h$alpha); print(class(h$prior_samples))
h2 <- hypothesis(f, "Trt1 = 0")
str(h2$hypothesis)
cat("\n== bayes_R2 ==\n")
print(bayes_R2(f)); print(bayes_R2(f, robust = TRUE))
str(bayes_R2(f, summary = FALSE))
cat("\n== as.data.frame / as.matrix / as.array names ==\n")
print(head(colnames(as.data.frame(f))))
print(dim(as.matrix(f))); print(dimnames(as.array(f))[2:3])
print(dim(as.matrix(f, variable = "b_Intercept")))
print(head(summary(f)$fixed))
print(summary(f)$random)
cat("DONE\n")
