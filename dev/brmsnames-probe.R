## What each call in lane brmsnames' scope does on one build. The same
## script is the BEFORE and the AFTER; nothing here is typed from memory.
##
##   Rscript dev/brmsnames-probe.R base > dev/brmsnames-log/probe-base.txt
##   Rscript dev/brmsnames-probe.R lane > dev/brmsnames-log/probe-lane.txt
##
## Fit side: the audit's epilepsy fit (brms's own data,
## count ~ zBase * Trt + (1 | patient), poisson), which is where items 8
## to 10 of dev/brms-suite-audit.md were found. Draws side:
## dev/stan-cache/brmsnames-draws-<arm>.rds from dev/brmsnames-draws.R.
arm <- commandArgs(trailingOnly = TRUE)[1L]
source("dev/brmsnames-libs.R")
brmsnames_libs(arm)
q <- function(e) suppressWarnings(suppressMessages(e))
q(library(frmtmb)); q(library(frmtmb.sample)); q(library(posterior))
cat("arm", arm, " frmtmb", format(packageVersion("frmtmb")), "from",
    dirname(system.file(package = "frmtmb")), "\n")
cat("frmtmb.sample", format(packageVersion("frmtmb.sample")), "from",
    dirname(system.file(package = "frmtmb.sample")), "\n\n")

shape <- function(v) {
  if (inherits(v, "try-error")) {
    m <- attr(v, "condition")$message
    return(paste0("ERROR: ", substr(gsub("\n", " ", m), 1, 150)))
  }
  d <- dim(v)
  paste0(class(v)[[1L]], " ",
         if (is.null(d)) paste0("len=", length(v)) else
           paste(d, collapse = "x"),
         if (!is.null(names(v)) && is.list(v))
           paste0(" names=", paste(names(v), collapse = ","))
         else "")
}
ev <- function(e) try(q(e), silent = TRUE)
show <- function(label, v) cat(sprintf("%-44s %s\n", label, shape(v)))

epilepsy <- brms::epilepsy
fit <- q(frm(count ~ zBase * Trt + (1 | patient), data = epilepsy,
             family = poisson()))

cat("== A8 posterior_summary on a fit ==\n")
show("posterior_summary(fit)", ev(posterior_summary(fit)))

cat("\n== A9 VarCorr on a fit ==\n")
v <- ev(VarCorr(fit))
show("VarCorr(fit)", v)
if (!inherits(v, "try-error")) {
  cat("names(VarCorr(fit)):", names(v), "\n")
  cat("class:", class(v), "\n")
  cat("element 1 names:", names(v[[1L]]), "\n")
  str(unclass(v[[1L]]))
}

cat("\n== A10 variables(fit) ==\n")
print(ev(variables(fit)))

cat("\n== A hypothesis(fit) object ==\n")
h <- ev(hypothesis(fit, "zBase > Trt1"))
show("hypothesis(fit, 'zBase > Trt1')", h)
if (!inherits(h, "try-error")) {
  cat("class:", class(h), "\nnames:", names(h), "\n")
  if (is.list(h$hypothesis)) str(h$hypothesis)
}
h2 <- ev(hypothesis(fit, "b_zBase > b_Trt1", class = NULL))
show("hypothesis(fit, 'b_zBase > b_Trt1', NULL)", h2)

cat("\n== B positional and summary = FALSE on a fit ==\n")
show("fixef(fit)", ev(fixef(fit)))
show("fixef(fit, FALSE)", ev(fixef(fit, FALSE)))
cat("  identical to fixef(fit):",
    identical(ev(fixef(fit, FALSE)), ev(fixef(fit))), "\n")
show("fixef(fit, summary = FALSE)", ev(fixef(fit, summary = FALSE)))
show("ranef(fit, FALSE)", ev(ranef(fit, FALSE)))
cat("  identical to ranef(fit):",
    identical(ev(ranef(fit, FALSE)), ev(ranef(fit))), "\n")
show("ranef(fit, summary = FALSE)", ev(ranef(fit, summary = FALSE)))
show("coef(fit, FALSE)", ev(coef(fit, FALSE)))
show("coef(fit, summary = FALSE)", ev(coef(fit, summary = FALSE)))
show("VarCorr(fit, NULL, FALSE)", ev(VarCorr(fit, NULL, FALSE)))
show("VarCorr(fit, summary = FALSE)", ev(VarCorr(fit, summary = FALSE)))

cat("\n== draws ==\n")
obj <- readRDS(sprintf("dev/stan-cache/brmsnames-draws-%s.rds", arm))
ds <- obj$ds
cat("draws", nrow(ds$draws), "x", ncol(ds$draws), "\n")
cat("variables(ds):", ev(variables(ds)), "\n")
cat("names(rhat(ds)):", names(ev(rhat(ds))), "\n")
cat("rownames(summary(ds)):", rownames(ev(summary(ds))), "\n")
cat("rownames(posterior_summary(ds)):",
    rownames(ev(posterior_summary(ds))), "\n")
cat("colnames(as.data.frame(ds)):", colnames(ev(as.data.frame(ds))), "\n")
cat("variables(as_draws_df(ds)):", variables(ev(as_draws_df(ds))), "\n")

cat("\n== B draws positional (the eleven and the six) ==\n")
pos <- list(
  c("as.matrix(ds, \"b_x\")", "as.matrix(ds)"),
  c("as.array(ds, \"b_x\")", "as.array(ds)"),
  c("as.data.frame(ds, NULL, TRUE, \"b_x\")", "as.data.frame(ds)"),
  c("as_draws_array(ds, \"b_x\")", "as_draws_array(ds)"),
  c("as_draws_df(ds, \"b_x\")", "as_draws_df(ds)"),
  c("as_draws_matrix(ds, \"b_x\")", "as_draws_matrix(ds)"),
  c("as_draws_list(ds, \"b_x\")", "as_draws_list(ds)"),
  c("as_draws_rvars(ds, \"b_x\")", "as_draws_rvars(ds)"),
  c("as_draws(ds, \"b_x\")", "as_draws(ds)"),
  c("fixef(ds, FALSE)", "fixef(ds)"),
  c("ranef(ds, FALSE)", "ranef(ds)"),
  c("coef(ds, FALSE)", "coef(ds)"),
  c("VarCorr(ds, NULL, FALSE)", "VarCorr(ds)"),
  c("summary(ds, FALSE, 0.5)", "summary(ds)"),
  c("bayes_R2(ds, NULL, TRUE, TRUE)", "bayes_R2(ds)"),
  c("posterior_summary(ds, \"^b_\")", "posterior_summary(ds)"),
  c("hypothesis(ds, \"x > 0\", \"b\")", "hypothesis(ds, \"x > 0\")"),
  c("pairs(ds, \"^b_\")", "pairs(ds)")
)
for (p in pos) {
  a <- ev(eval(str2lang(p[1L])))
  b <- ev(eval(str2lang(p[2L])))
  same <- !inherits(a, "try-error") && !inherits(b, "try-error") &&
    identical(a, b)
  cat(sprintf("%-40s %-48s %s\n", p[1L], substr(shape(a), 1, 48),
              if (same) "IDENTICAL to no-arg" else ""))
}
cat("DONE\n")
