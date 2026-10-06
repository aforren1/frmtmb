# Reviewer: brms's Stan log density against frmtmb's joint density at
# frmtmb's estimates on gp(by = ) forms the lane did not list, plus
# three of the lane's forms as a reproduction. Flat priors both sides.
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
  sys.source(h, envir = env)
}
src <- parse(file.path(wt, "tests/testthat/test-gp-by.R"))
for (e in src) {
  if (is.call(e) && identical(e[[1]], as.name("<-"))) eval(e, env)
}
d <- env$gpby_data()
set.seed(31)
d$y2 <- 1 + cos(d$x) * (d$f == "a") + stats::rnorm(nrow(d), 0, 0.4)
d$g4 <- factor(rep(c("p", "q", "r", "s"), length.out = nrow(d)))
d$y4 <- d$y + 0.6 * sin(d$x) * (d$g4 == "q") - 0.4 * (d$g4 == "s") *
  cos(d$x)
d$fs <- d$g4
stats::contrasts(d$fs) <- stats::contr.sum(4)
d$ysig <- 0.5 + sin(d$x) + stats::rnorm(nrow(d), 0,
                                        exp(-1 + 0.3 * (d$f == "b") *
                                              sin(d$x)))
dsub <- d[d$f != "b", ]   # keeps level "b" with no rows
cat("dsub levels", levels(dsub$f), "rows per level",
    table(dsub$f), "\n")

run <- function(lab, bform, fform, family, data) {
  cat("==", lab, "\n")
  fit <- try(frm(fform, data = data, family = family))
  if (inherits(fit, "try-error")) {
    cat(sprintf("RES %s | FRM ERROR\n", lab)); return(invisible())
  }
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
run("lane: gp(x, by = f, k = 8)",
    brms::bf(y ~ gp(x, by = f, k = 8)), bf(y ~ gp(x, by = f, k = 8)), G, d)
run("lane: gp(x, by = f, k = 8, cmc = FALSE)",
    brms::bf(y ~ gp(x, by = f, k = 8, cmc = FALSE)),
    bf(y ~ gp(x, by = f, k = 8, cmc = FALSE)), G, d)
run("lane: gp(x, by = w, k = 8)",
    brms::bf(y ~ gp(x, by = w, k = 8)), bf(y ~ gp(x, by = w, k = 8)), G, d)
run("new: sigma ~ gp(x, by = f, k = 6)",
    brms::bf(ysig ~ x, sigma ~ gp(x, by = f, k = 6)),
    bf(ysig ~ x, sigma ~ gp(x, by = f, k = 6)), G, d)
run("new: gp(x, by = f, k = 8) + (1 | f)",
    brms::bf(y ~ gp(x, by = f, k = 8) + (1 | f)),
    bf(y ~ gp(x, by = f, k = 8) + (1 | f)), G, d)
run("new: gp(x, by = g4, k = 6, cmc = FALSE), 4 levels",
    brms::bf(y4 ~ gp(x, by = g4, k = 6, cmc = FALSE)),
    bf(y4 ~ gp(x, by = g4, k = 6, cmc = FALSE)), G, d)
run("new: gp(x, by = fs, k = 6, cmc = FALSE), contr.sum",
    brms::bf(y4 ~ gp(x, by = fs, k = 6, cmc = FALSE)),
    bf(y4 ~ gp(x, by = fs, k = 6, cmc = FALSE)), G, d)
run("new: empty level after subsetting, gp(x, by = f, k = 8)",
    brms::bf(y ~ gp(x, by = f, k = 8)), bf(y ~ gp(x, by = f, k = 8)), G,
    dsub)
run("new: empty level, cmc = FALSE",
    brms::bf(y ~ gp(x, by = f, k = 8, cmc = FALSE)),
    bf(y ~ gp(x, by = f, k = 8, cmc = FALSE)), G, dsub)
run("new: gp(x, by = f, k = 8, scale = FALSE)",
    brms::bf(y ~ gp(x, by = f, k = 8, scale = FALSE)),
    bf(y ~ gp(x, by = f, k = 8, scale = FALSE)), G, d)
run("new: gp(x, by = f, k = 8, c = 2)",
    brms::bf(y ~ gp(x, by = f, k = 8, c = 2)),
    bf(y ~ gp(x, by = f, k = 8, c = 2)), G, d)
run("new: gp(x, by = f, k = 8, c = 1.1)",
    brms::bf(y ~ gp(x, by = f, k = 8, c = 1.1)),
    bf(y ~ gp(x, by = f, k = 8, c = 1.1)), G, d)
run("new: mv gp(x, by = f, k = 8) | gp(x, k = 6)",
    brms::mvbf(brms::bf(y ~ gp(x, by = f, k = 8)),
               brms::bf(y2 ~ gp(x, k = 6))) + brms::set_rescor(FALSE),
    mvbf(bf(y ~ gp(x, by = f, k = 8)), bf(y2 ~ gp(x, k = 6))) +
      set_rescor(FALSE), G, d)
run("new: mv both by = f",
    brms::mvbf(brms::bf(y ~ gp(x, by = f, k = 8)),
               brms::bf(y2 ~ gp(x, by = f, k = 6))) +
      brms::set_rescor(FALSE),
    mvbf(bf(y ~ gp(x, by = f, k = 8)), bf(y2 ~ gp(x, by = f, k = 6))) +
      set_rescor(FALSE), G, d)
cat("DONE\n")
