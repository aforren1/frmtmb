## Reviewer: behavioural probes on the draws bmrev-diag.R saved.
## No refit. Covers:
##  1. summary(ds) against rhat(ds)/neff_ratio(ds)
##  2. the pars= semantics of rhat/neff_ratio against brms's
##  3. the two that "answered silently": refusal text, and over-refusal
##  4. contracts the truncate-at-dots criterion cannot see
LIB <- "C:/Users/adf44/source/r/bmrev-lib"
.libPaths(c(LIB,
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
q <- function(e) suppressWarnings(suppressMessages(e))
q(library(frmtmb)); q(library(frmtmb.sample))
q(library(posterior)); q(library(brms))
ds <- readRDS("dev/stan-cache/bmrev-draws.rds")
cat("draws ", nrow(ds$draws), " x ", ncol(ds$draws), "\n\n")

shape <- function(v) {
  if (inherits(v, "try-error"))
    return(paste0("ERROR: ", sub("\n.*$", "",
                                 conditionMessage(attr(v, "condition")))))
  d <- dim(v)
  paste0(paste(class(v), collapse = "/"), " ",
         if (is.null(d)) paste0("len=", length(v)) else
           paste(d, collapse = "x"))
}
go <- function(lbl, e) {
  v <- try(q(e), silent = TRUE)
  cat(sprintf("%-52s %s\n", lbl, shape(v)))
  invisible(v)
}

## ---- 1. summary(ds) vs the two moved diagnostics --------------------
cat("== 1. summary(ds) against rhat()/neff_ratio() ==\n")
s <- summary(ds)
cat("colnames(summary(ds)): ", paste(colnames(s), collapse = " "), "\n")
cat("rownames(summary(ds)): ", paste(rownames(s), collapse = " "), "\n")
r <- rhat(ds); n <- neff_ratio(ds)
if ("Rhat" %in% colnames(s)) {
  k <- intersect(rownames(s), names(r))
  cat("rows shared with rhat(ds): ", length(k), " of ", nrow(s), "\n")
  d <- data.frame(var = k,
                  summary_Rhat = unname(s[k, "Rhat"]),
                  rhat_ds = unname(r[k]))
  d$diff <- d$summary_Rhat - d$rhat_ds
  print(d, digits = 8, row.names = FALSE)
  cat("max |summary Rhat - rhat(ds)|: ",
      format(max(abs(d$diff)), digits = 8), "\n")
  cat("identical: ", identical(unname(s[k, "Rhat"]), unname(r[k])), "\n")
}
for (cn in c("n_eff", "Bulk_ESS", "Tail_ESS", "ESS")) {
  if (cn %in% colnames(s)) {
    k <- intersect(rownames(s), names(n))
    cat("\nsummary column ", cn, " vs neff_ratio(ds)*ndraws:\n", sep = "")
    d <- data.frame(var = k, summary = unname(s[k, cn]),
                    neff_ratio_x_N = unname(n[k]) * nrow(ds$draws))
    d$diff <- d$summary - d$neff_ratio_x_N
    print(d, digits = 8, row.names = FALSE)
  }
}
cat("\nprint(ds) header columns:\n")
print(utils::head(utils::capture.output(print(ds)), 30))

## ---- 2. pars= semantics against brms's rhat ------------------------
cat("\n== 2. pars= on rhat()/neff_ratio() ==\n")
cat("brms::rhat.brmsfit passes variable = pars to as_draws_array with\n")
cat("regex = FALSE, so a regex `pars` is an ERROR in brms.\n")
arr <- posterior::as_draws_array(ds)
cat("posterior::subset_draws(arr, variable = \"^x$\") -> ")
cat(shape(try(posterior::subset_draws(arr, variable = "^x$"),
               silent = TRUE)), "\n")
cat("posterior::subset_draws(arr, variable = \"x\")   -> ")
cat(shape(try(posterior::subset_draws(arr, variable = "x"),
               silent = TRUE)), "\n")
cat("brms:::as_draws_array.brmsfit formals: ",
    paste(names(formals(brms:::as_draws_array.brmsfit)),
          collapse = " "), "\n")
go("rhat(ds, \"^x$\")        [brms: ERROR]", rhat(ds, "^x$"))
go("rhat(ds, \"x\")          [brms: 1 variable]", rhat(ds, "x"))
go("rhat(ds, NULL)         [brms: ALL variables]", rhat(ds, NULL))
go("rhat(ds, NA)           [brms: ERROR]", rhat(ds, NA))
go("neff_ratio(ds, NULL)   [brms: ALL variables]", neff_ratio(ds, NULL))
cat("rhat(ds, \"x\") names: ",
    paste(names(q(rhat(ds, "x"))), collapse = " "), "\n")

## ---- 3. the two that answered silently ------------------------------
cat("\n== 3. the refusal, character for character ==\n")
ap <- variables(ds)
bm <- function(p) tryCatch({
  brms:::extract_pars(p, all_pars = ap, fixed = FALSE); "(no error)" },
  error = function(e) conditionMessage(e))
om <- function(e) tryCatch({ q(e); "(no error)" },
                           error = function(c) conditionMessage(c))
for (p in list(TRUE, 0.9, NULL, 1L)) {
  lbl <- paste(deparse(p), collapse = "")
  cat(sprintf("pars=%-6s brms: %s\n", lbl, bm(p)))
}
cat("\nas.mcmc(ds, TRUE)          ours: ", om(as.mcmc(ds, TRUE)), "\n")
cat("posterior_interval(ds,0.9) ours: ",
    om(posterior_interval(ds, 0.9)), "\n")
cat("as.mcmc(ds, NULL)          ours: ", om(as.mcmc(ds, NULL)), "\n")
cat("identical to brms's string: ",
    identical(om(as.mcmc(ds, TRUE)), bm(TRUE)), "\n")
cat("\n-- over-refusal check: a character pars must be ANSWERED --\n")
go("as.mcmc(ds, \"^x$\")", as.mcmc(ds, "^x$"))
go("posterior_interval(ds, \"^x$\")", posterior_interval(ds, "^x$"))
pi1 <- q(posterior_interval(ds, "^x$"))
cat("   rownames: ", paste(rownames(pi1), collapse = " "), "\n")
cat("   colnames: ", paste(colnames(pi1), collapse = " "), "\n")
mc <- q(as.mcmc(ds, "^x$"))
cat("   as.mcmc varnames: ",
    paste(colnames(mc[[1L]]), collapse = " "), "\n")
go("as.mcmc(ds, NA)", as.mcmc(ds, NA))
go("as.mcmc(ds, c(\"x\",\"sigma_Intercept\"), TRUE)",
   as.mcmc(ds, c("x", "sigma_Intercept"), TRUE))
go("posterior_interval(ds, NA, \"x\", 0.9)",
   posterior_interval(ds, NA, "x", 0.9))
go("mcmc_plot(ds, \"^x$\")", mcmc_plot(ds, "^x$"))

## ---- 3b. posterior_interval's default variable set ------------------
cat("\n== 3b. posterior_interval() default variables ==\n")
piD <- q(posterior_interval(ds))
cat("ours default rows:  ", paste(rownames(piD), collapse = " "), "\n")
cat("variables(ds):      ", paste(variables(ds), collapse = " "), "\n")
cat("dropped:            ",
    paste(setdiff(variables(ds), rownames(piD)), collapse = " "), "\n")
cat("brms would return every one of them (as.matrix(object, pars=NA)).\n")

## ---- 4. contracts the truncate-at-dots criterion cannot see ---------
cat("\n== 4. positional brms calls that land in our `...` ==\n")
cat("These 'agree as far as both go' only because OUR side reaches\n")
cat("`...` first. A positional brms call is silently ignored.\n\n")
m0 <- q(as.matrix(ds))
go("as.matrix(ds)                [brms: all vars]", as.matrix(ds))
mm <- go("as.matrix(ds, \"b_x\")         [brms: 1 column]",
         as.matrix(ds, "b_x"))
cat("   ours ncol unchanged by pars: ",
    identical(dim(mm), dim(m0)), "\n")
go("as.array(ds, \"x\")            [brms: 1 variable]",
   as.array(ds, "x"))
go("as_draws_array(ds, \"x\")      [brms: 1 variable]",
   as_draws_array(ds, "x"))
go("as_draws_df(ds, \"x\")         [brms: 1 variable]",
   as_draws_df(ds, "x"))
f1 <- go("fixef(ds)                    [brms: summary matrix]", fixef(ds))
f2 <- go("fixef(ds, FALSE)             [brms: RAW DRAWS matrix]",
         fixef(ds, FALSE))
cat("   fixef(ds, FALSE) identical to fixef(ds): ",
    identical(f1, f2), "\n")
r1 <- go("ranef(ds)                    [brms: summary array]", ranef(ds))
r2 <- go("ranef(ds, FALSE)             [brms: RAW DRAWS array]",
         ranef(ds, FALSE))
cat("   ranef(ds, FALSE) identical to ranef(ds): ",
    identical(r1, r2), "\n")
c1 <- go("coef(ds)", coef(ds))
c2 <- go("coef(ds, FALSE)              [brms: RAW DRAWS]", coef(ds, FALSE))
cat("   coef(ds, FALSE) identical to coef(ds): ",
    identical(c1, c2), "\n")
v1 <- go("VarCorr(ds)", VarCorr(ds))
v2 <- go("VarCorr(ds, NULL, FALSE)     [brms: RAW DRAWS]",
         VarCorr(ds, NULL, FALSE))
cat("   identical: ", identical(v1, v2), "\n")
p1 <- go("posterior_samples(ds, \"x\")", posterior_samples(ds, "x"))
p2 <- go("posterior_samples(ds, \"x\", TRUE) [brms: fixed=TRUE]",
         posterior_samples(ds, "x", TRUE))
cat("   identical: ", identical(p1, p2), "\n")
go("nsamples(ds, 10)             [brms: subset]", nsamples(ds, 10))
go("plot(ds, \"^x$\")              [brms: pars]", plot(ds, "^x$"))
go("summary(ds, NULL, 0.5)       [brms: priors, prob]",
   summary(ds, NULL, 0.5))
go("posterior_summary(ds, \"x\")   [brms: pars]",
   posterior_summary(ds, "x"))
go("loo(ds, ds)                  [brms: compare two models]",
   loo(ds, ds))
go("waic(ds, ds)                 [brms: compare two models]",
   waic(ds, ds))
go("pairs(ds, \"^x$\")             [brms: pars]", pairs(ds, "^x$"))
go("bayes_R2(ds, NULL, TRUE, TRUE) [brms: robust at 4]",
   bayes_R2(ds, NULL, TRUE, TRUE))
cat("DONE\n")
