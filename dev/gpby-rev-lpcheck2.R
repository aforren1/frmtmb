# Reviewer: check C on a gp(by = ) in sigma and in a multivariate model.
# The test translator (helper-brms.R) does not name a GP outside mu of a
# univariate model; it is patched HERE, in this process only, to read
# brms's dpar- and response-suffixed GP data (Kgp_sigma_1, Kgp_y_1).
.libPaths(c("C:/Users/adf44/source/r/wt-gpby-lib",
            "C:/Users/adf44/source/r/rellib-r5",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-gpby"
Sys.setenv(NOT_CRAN = "true", FRMTMB_BRMS_FIT_TESTS = "true",
           R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win",
           FRMTMB_STAN_CACHE = file.path(wt, "dev/gpby-rev-stan-cache"))
suppressPackageStartupMessages({library(frmtmb); library(testthat)})
cat("lib:", find.package("frmtmb"), "\n")
env <- testthat::test_env("frmtmb")
for (h in list.files(file.path(wt, "tests/testthat"), "^helper-.*[.]R$",
                     full.names = TRUE)) {
  txt <- readLines(h)
  if (basename(h) == "helper-brms.R") {
    i0 <- grep("i <- ip$idx[[1]]", txt, fixed = TRUE)[1]
    txt <- append(txt, paste0("      gsfx <- if (identical(ip$dpar, \"mu\"))",
                              " \"\" else paste0(ip$dpar, \"_\")"), i0)
    for (p in c("Kgp_", "Dgp_", "Igp_", "dmax_")) {
      txt <- gsub(paste0("paste0(\"", p, "\", i"),
                  paste0("paste0(\"", p, "\", gsfx, i"), txt, fixed = TRUE)
    }
    j0 <- grep("b[[\"covstruct\"]] %in% cs && identical(brms_block_dpar(b), dpar)",
               txt, fixed = TRUE)
    txt[j0] <- paste0(
      "    dd <- brms_block_dpar(b); ",
      "rr <- sub(\"[.][^.]*$\", \"\", b$components[[1]]$lp_key);",
      " b[[\"covstruct\"]] %in% cs && (identical(dd, dpar) || ",
      "(nzchar(rr) && identical(if (dd == \"mu\") rr else ",
      "paste0(dd, \"_\", rr), dpar)))")
    k0 <- grep("zgp_\\\\d+(_\\\\d+)?", txt, fixed = TRUE)
    txt <- gsub("zgp_\\\\d+(_\\\\d+)?", "zgp_([a-z0-9]+_)?[0-9]+(_[0-9]+)?",
                txt, fixed = TRUE)
    cat("patched lines:", i0, j0, k0, "\n")
  }
  eval(parse(text = txt), envir = env)
}
src <- parse(file.path(wt, "tests/testthat/test-gp-by.R"))
for (e in src) {
  if (is.call(e) && identical(e[[1]], as.name("<-"))) eval(e, env)
}
d <- env$gpby_data()
set.seed(31)
d$y2 <- 1 + cos(d$x) * (d$f == "a") + stats::rnorm(nrow(d), 0, 0.4)
d$ysig <- 0.5 + sin(d$x) + stats::rnorm(nrow(d), 0,
                                        exp(-1 + 0.3 * (d$f == "b") *
                                              sin(d$x)))
run <- function(lab, bform, fform, family, data) {
  cat("==", lab, "\n")
  fit <- frm(fform, data = data, family = family)
  r <- try(local({
    environment(env$brms_lp_check) <- env
    env$brms_lp_check(bform, family, data, fit, joint = TRUE)
  }))
  if (inherits(r, "try-error")) {
    cat(sprintf("RES %s | CHECK ERROR %s\n", lab,
                conditionMessage(attr(r, "condition"))))
  } else {
    cat(sprintf("RES %s | const %.3e | max_grad %.3e | ours %.6f\n", lab,
                r$measured_const, r$max_grad, r$ours))
  }
}
G <- gaussian()
run("lane form under the patch: gp(x, by = f, k = 8)",
    brms::bf(y ~ gp(x, by = f, k = 8)), bf(y ~ gp(x, by = f, k = 8)), G, d)
run("new: sigma ~ gp(x, by = f, k = 6)",
    brms::bf(ysig ~ x, sigma ~ gp(x, by = f, k = 6)),
    bf(ysig ~ x, sigma ~ gp(x, by = f, k = 6)), G, d)
run("new: mv gp(x, by = f, k = 8) | gp(x, k = 6)",
    brms::mvbf(brms::bf(y ~ gp(x, by = f, k = 8)),
               brms::bf(y2 ~ gp(x, k = 6))) + brms::set_rescor(FALSE),
    mvbf(bf(y ~ gp(x, by = f, k = 8)), bf(y2 ~ gp(x, k = 6))) +
      set_rescor(FALSE), G, d)
run("new: mv both by = f, k = 8 and k = 6",
    brms::mvbf(brms::bf(y ~ gp(x, by = f, k = 8)),
               brms::bf(y2 ~ gp(x, by = f, k = 6))) +
      brms::set_rescor(FALSE),
    mvbf(bf(y ~ gp(x, by = f, k = 8)), bf(y2 ~ gp(x, by = f, k = 6))) +
      set_rescor(FALSE), G, d)
cat("DONE\n")
