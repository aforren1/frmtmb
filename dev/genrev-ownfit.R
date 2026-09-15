LIB <- commandArgs(trailingOnly = TRUE)[1]
.libPaths(c(LIB, "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
q <- function(...) suppressMessages(suppressWarnings(...))
q(library(posterior)); q(library(loo)); q(library(bayesplot))
q(library(frmtmb))
say <- function(...) cat(sprintf(...))
own <- function(f) { e <- environment(f); if (is.null(e)) return("base")
  n <- environmentName(topenv(e)); if (nzchar(n)) n else "<anon>" }
set.seed(2026)
dd <- data.frame(x = rnorm(200), g = factor(rep(1:10, 20)))
dd$y <- exp(rnorm(200, 8 + 0.4 * dd$x, 0.4))
f <- q(frm(bf(y ~ x + (1 | g)) + lognormal(), data = dd))
say("LIB %s\n", LIB)
gens <- c("as_draws","as_draws_df","as_draws_list","as_draws_rvars",
          "as_draws_array","as_draws_matrix","nvariables","niterations",
          "nchains","ndraws","variables","loo_compare","posterior_summary",
          "loo","waic","LOO","WAIC","bayes_R2","prior_summary",
          "hypothesis","conditional_effects","expose_functions","pp_check",
          "ngrps","fixef","ranef","VarCorr")
for (g in gens) {
  gf <- tryCatch(get(g), error = function(e) NULL)
  o <- if (is.function(gf)) own(gf) else "<absent>"
  m <- tryCatch(utils::getS3method(g, "frmtmb_fit", optional = TRUE),
                error = function(e) NULL)
  mo <- if (is.null(m)) "<none>" else own(m)
  v <- tryCatch({
        x <- do.call(g, list(f))
        paste0("OK:", class(x)[1])
      }, error = function(e) paste0("ERR:", substr(conditionMessage(e), 1, 55)))
  say("%-20s gen=%-11s meth=%-11s %s\n", g, o, mo, v)
}
cat("GENREVDONE\n")
