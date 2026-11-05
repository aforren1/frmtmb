# Item 9: mutate each fix in this process and run the lane's new test file with
# the mutant: a test that passes under a mutant does not pin the fix. One
# mutant per process. Usage: Rscript dev/review-refuse-mutate.R <id>
# Also base: run a new test file on rellib-r16 for the behavioral failure.
id <- commandArgs(trailingOnly = TRUE)[1L]
source("dev/review-refuse-lib.R")
library(frmtmb.sample)
ns <- asNamespace("frmtmb")
mutate <- function(name, fn, nsp = "frmtmb") {
  assignInNamespace(name, fn, nsp)
}
tf <- switch(id,
  jac_off = "tests/testthat/test-mi-families.R",
  jac_off_l1 = "tests/testthat/test-mi-mm-brms-l1.R",
  chart_id = "tests/testthat/test-mi-families.R",
  chart_id_l1 = "tests/testthat/test-mi-mm-brms-l1.R",
  no_laplace_guard = "tests/testthat/test-mi-families.R",
  sds_class_reach = "tests/testthat/test-prior-sds.R",
  sds_term = "tests/testthat/test-prior-sds.R",
  sd_reach_smooth = "tests/testthat/test-prior-sds.R",
  gr_cov_ident = "tests/testthat/test-mm-id-cov.R",
  gr_cov_ident_l1 = "tests/testthat/test-mi-mm-brms-l1.R",
  mix_lcdf = "tests/testthat/test-refuse-pairs.R",
  mix_lcdf_l1 = "tests/testthat/test-mi-mm-brms-l1.R",
  se_resid = "tests/testthat/test-refuse-pairs.R",
  se_resid_l1 = "tests/testthat/test-mi-mm-brms-l1.R",
  plan_sds = "extensions/frmtmb.sample/tests/testthat/test-sds-prior.R",
  default_sds = "extensions/frmtmb.sample/tests/testthat/test-sds-prior.R",
  stop("unknown"))
pkg <- if (grepl("frmtmb.sample", tf)) "frmtmb.sample" else "frmtmb"
switch(id,
  jac_off = , jac_off_l1 = {
    orig <- ns$mi_chart_logjac
    mutate("mi_chart_logjac", function(u, chart) orig(u, chart) * 0)
  },
  chart_id = , chart_id_l1 = {
    # every family on the identity chart; a value can leave the support
    mutate("mi_chart_maps", function(chart) FALSE)
    mutate("mi_chart_fun", function(y, chart) y)
    mutate("mi_chart_inv", function(u, chart) u)
    mutate("mi_placeholder", function(obs, chart) mean(obs))
  },
  no_laplace_guard = mutate("check_mi_laplace", function(frame) invisible(NULL)),
  sds_class_reach = {
    # sds is treated as sd: class sds stops reaching a smooth, sd reaches it
    orig <- ns$sd_spec_reach
    fn <- orig
    body(fn)[[5]] <- quote(NULL)
    mutate("sd_spec_reach", fn)
  },
  sds_term = {
    # coef of sds compared loosely
    mutate("block_sm_term", function(bk) bk[["group_name"]])
  },
  sd_reach_smooth = {
    orig <- ns$sd_spec_reach
    mutate("sd_spec_reach", function(spec, frame, bk, s) {
      if (identical(bk[["covstruct"]], "smooth") && identical(s$class, "sd")) {
        s$class <- "sds"
      }
      orig(spec, frame, bk, s)
    })
  },
  gr_cov_ident = , gr_cov_ident_l1 = {
    orig <- ns$gr_cov_matrix
    mutate("gr_cov_matrix", function(expr, data2, data, env, lv, what) {
      A <- orig(expr, data2, data, env, lv, what)
      if (startsWith(what, "mm(")) diag(nrow(A)) else A
    })
  },
  mix_lcdf = , mix_lcdf_l1 = {
    # the bug the lane first wrote: lcdf a log-mixture, as lccdf
    orig <- ns$fam_lcdf
    mutate("fam_lcdf", function(fam, q, dpars, aterms, trunc = NULL) {
      if (!is.null(fam[["mix"]])) return(fam$lccdf(q, dpars, aterms))
      orig(fam, q, dpars, aterms, trunc)
    })
  },
  se_resid = , se_resid_l1 = {
    orig <- ns$resid_sd
    mutate("resid_sd", function(sigma, aterms) {
      if (!is.null(aterms[["se"]])) return(sigma)
      orig(sigma, aterms)
    })
  },
  plan_sds = {
    nss <- asNamespace("frmtmb.sample")
    orig <- nss$prior_plan
    fn <- orig
    txt <- deparse(fn)
    txt <- sub('tab\\$class == "sds" & tab\\$coef == sm_term\\[j\\]', 'tab$class == "sd" & tab$coef == sm_term[j]', txt)
    txt <- sub('tab\\$class == "sds" & !nzchar\\(tab\\$coef\\) & own', 'tab$class == "sd" & !nzchar(tab$coef) & own', txt)
    stopifnot(!identical(txt, deparse(fn)))
    fn2 <- eval(parse(text = txt)); environment(fn2) <- nss
    mutate("prior_plan", fn2, "frmtmb.sample")
  },
  default_sds = {
    nss <- asNamespace("frmtmb.sample")
    orig <- nss$default_sd_prefixes
    mutate("default_sd_prefixes", function(fit, cls = "sd") {
      if (identical(cls, "sds")) return(list())
      orig(fit, cls)
    }, "frmtmb.sample")
  })
env <- testthat::test_env(pkg)
library(pkg, character.only = TRUE)
res <- as.data.frame(testthat::test_file(tf, package = pkg, env = env, reporter = "silent", stop_on_failure = FALSE))
pr("MUT %-20s %-40s pass=%d fail=%d error=%d skip=%d; failing: %s", id, basename(tf), sum(res$passed), sum(res$failed), sum(res$error), sum(res$skipped),
   paste(substr(res$test[res$failed > 0 | res$error], 1, 60), collapse = " | "))
