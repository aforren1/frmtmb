# Reviewer of lane defects: mutate one fix in the loaded namespace (no
# install) and run one test file against it. A mutant the tests do not
# catch is a gap in the tests.
#   Rscript dev/defects-rev-mutant.R <mutant> <pkg> <test file>
a <- commandArgs(trailingOnly = TRUE)
mut <- a[1]; p <- a[2]; tf <- a[3]
.libPaths(c("C:/Users/adf44/source/r/wt-defects-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true",
           FRMTMB_BRMSPORT_PKG = p,
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressMessages({library(testthat); library(frmtmb)})
mutants <- list(
  M1 = list("ord_linear_per_threshold", "as.numeric(eta) + CS",
            "as.numeric(eta) - CS"),
  M2 = list("prior_coef_rows", "out[order(out)]", "out"),
  M3 = list("predict_new_level_draw", "shared[[key]] <- if (is.null(pick)) ",
            "shared[[key]] <- if (TRUE) "),
  M4 = list("check_newdata_frame", 'attr(x, "contrasts") <- NULL', "NULL"),
  M5 = list("residuals_newdata", 'f[, "Est.Error"]', '0 * f[, "Est.Error"]'),
  M6 = list("fitted_point", 'if (scale == "linear" &&', "if (FALSE &&"),
  M7 = list("predict_new_level_spec",
            "pick[[key]] <- sample.int(bl[[key]]$bk[[\"n_levels\"]], 1L)",
            "pick[[key]] <- 1L"),
  M8 = list("summary_gp_frame", "tb$est_t[rg] <- tb$est_t[rg] - log(dmax)",
            "NULL")
)
m <- mutants[[mut]]
ns <- asNamespace("frmtmb")
f <- get(m[[1]], envir = ns)
src <- deparse(f, width.cutoff = 500L)
hit <- grepl(m[[2]], src, fixed = TRUE)
stopifnot("mutation site not found" = sum(hit) >= 1L)
src <- gsub(m[[2]], m[[3]], src, fixed = TRUE)
g <- eval(parse(text = src), envir = ns)
environment(g) <- ns
utils::assignInNamespace(m[[1]], g, ns = "frmtmb")
stopifnot(!identical(deparse(get(m[[1]], envir = ns)), deparse(f)))
suppressMessages(library(p, character.only = TRUE))
res <- as.data.frame(test_file(tf, package = p, env = test_env(p),
                               reporter = "silent", stop_on_failure = FALSE))
cat("MUTANT", mut, m[[1]], basename(tf), "pass=", sum(res$passed),
    "fail=", sum(res$failed), "err=", sum(res$error), "\n")
bad <- res[res$failed > 0 | res$error, "test"]
if (length(bad)) cat("  caught by:", paste(bad, collapse = " | "), "\n")
