# Reproduce the seven defects of lane wt-adefects on the BASE build
# (rellib-r3: frmtmb 0.60.0, frmtmb.sample 0.8.0, frmtmb.latent 0.4.0),
# with brms 2.23.0 on the same call wherever brms can answer without
# compiling Stan.
#
#   Rscript dev/adefects-repro-base.R > dev/adefects-log/repro-base.txt 2>&1
#
# Data seed 20260917 throughout, which is the seed
# dev/brmsport-defects.R uses, so the numbers are comparable to the ones
# in dev/brmsport-findings.md section 5.
# The reference build stays on the path BEHIND whichever build is under
# test, so that frmtmb.latent and the other extensions come from it in
# both arms. Without it the user library answered, with frmtmb.latent
# 0.3.0, and D5 read as a live defect that 0.60.0 had already fixed.
LIB <- Sys.getenv("ADEFECTS_LIB", "C:/Users/adf44/source/r/rellib-r3")
.libPaths(c(LIB, "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({
  library(frmtmb)
})
cat("lib", LIB, "\n")
cat("frmtmb", format(packageVersion("frmtmb")),
    "frmtmb.sample", format(packageVersion("frmtmb.sample")),
    "frmtmb.latent", format(packageVersion("frmtmb.latent")),
    "brms", format(packageVersion("brms")), "\n")

show <- function(label, expr) {
  warns <- character()
  out <- tryCatch({
    v <- withCallingHandlers(expr, warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }, message = function(m) invokeRestart("muffleMessage"))
    paste(utils::capture.output(print(v)), collapse = " | ")
  }, error = function(e) {
    paste0("ERROR[", paste(class(e), collapse = ","), "]: ",
           conditionMessage(e))
  })
  cat(sprintf("  %-44s %s\n", label, substr(gsub(" +", " ", out), 1, 260)))
  if (length(warns)) {
    cat(sprintf("  %-44s [warning] %s\n", "",
                substr(unique(warns)[1], 1, 180)))
  }
  invisible(NULL)
}

# ---------------------------------------------------------------- D1/D2
# The design of dev/brmsport-defects.R sections S2 and S3, so the
# logLik values here are comparable to the ones recorded there.
set.seed(20260917)
d <- data.frame(g1 = rep(c(1, 1, 2, 2), each = 6),
                g2 = rep(c(1, 2, 1, 2), each = 6),
                t = c(1:6, 1:6, 1:6, 7:12))
d$x <- rnorm(24)
d$g <- interaction(d$g1, d$g2)
d$y <- d$x + as.numeric(stats::filter(rnorm(24), 0.6, "recursive"))
d$gf1 <- factor(d$g1)
d$gf2 <- factor(d$g2)

cat("\n== D1 ar()/ma() accept an EXPRESSION as the time index\n")
show("frmtmb ar(t, g, cov = TRUE) control",
     round(logLik(frm(y ~ x + ar(t, g, cov = TRUE), d)), 6))
show("frmtmb ar(x + t, g, cov = TRUE)",
     round(logLik(frm(y ~ x + ar(x + t, g, cov = TRUE), d)), 6))
show("frmtmb ar(t - 10 * x, g, cov = TRUE)",
     round(logLik(frm(y ~ x + ar(t - 10 * x, g, cov = TRUE), d)), 6))
show("frmtmb ma(t, g, cov = TRUE) control",
     round(logLik(frm(y ~ x + ma(t, g, cov = TRUE), d)), 6))
show("frmtmb ma(x + t, g, cov = TRUE)",
     round(logLik(frm(y ~ x + ma(x + t, g, cov = TRUE), d)), 6))
show("frmtmb arma(x + t, g, cov = TRUE)",
     round(logLik(frm(y ~ x + arma(x + t, g, cov = TRUE), d)), 6))
show("frmtmb cosy(x + t, g)",
     round(logLik(frm(y ~ x + cosy(x + t, g), d)), 6))
show("frmtmb unstr(x + t, g)",
     round(logLik(frm(y ~ x + unstr(x + t, g), d)), 6))
show("brms ar(x + t, g)",
     names(brms::make_standata(y ~ x + ar(x + t, g, cov = TRUE), d)))
show("brms ma(x + t, g)",
     names(brms::make_standata(y ~ x + ma(x + t, g, cov = TRUE), d)))
show("brms arma(x + t, g)",
     names(brms::make_standata(y ~ x + arma(x + t, g, cov = TRUE), d)))
show("brms cosy(x + t, g)",
     names(brms::make_standata(y ~ x + cosy(x + t, g), d)))
show("brms unstr(x + t, g)",
     names(brms::make_standata(y ~ x + unstr(x + t, g), d)))
cat("  what brms ACCEPTS as time (control: bare names and NA)\n")
show("brms ar(t, g)", length(brms::make_standata(
  y ~ x + ar(t, g, cov = TRUE), d)))
show("brms ar(gr = g)", length(brms::make_standata(
  y ~ x + ar(gr = g, cov = TRUE), d)))
show("brms ar(NA, g)", length(brms::make_standata(
  y ~ x + ar(NA, g, cov = TRUE), d)))
show("brms ar(`odd name`, g)", tryCatch(length(brms::make_standata(
  y ~ x + ar(`t`, g, cov = TRUE), d)), error = conditionMessage))

cat("\n== D2 gr = grammar\n")
show("frmtmb ar(t, gr = g1/g2, cov = TRUE)",
     round(logLik(frm(y ~ x + ar(t, gr = g1/g2, cov = TRUE), d)), 6))
show("frmtmb ar(t, gr = g1:g2, cov = TRUE)",
     round(logLik(frm(y ~ x + ar(t, gr = g1:g2, cov = TRUE), d)), 6))
show("frmtmb ar(t, gr = gf1:gf2, cov = TRUE) factors",
     round(logLik(frm(y ~ x + ar(t, gr = gf1:gf2, cov = TRUE), d)), 6))
show("frmtmb ar(t, gr = interaction(g1,g2), cov = TRUE)",
     round(logLik(frm(y ~ x + ar(t, gr = interaction(g1, g2),
                                 cov = TRUE), d)), 6))
show("frmtmb ar(t, gr = g1 + g2, cov = TRUE)",
     round(logLik(frm(y ~ x + ar(t, gr = g1 + g2, cov = TRUE), d)), 6))
show("frmtmb ar(t, gr = g1 * g2, cov = TRUE)",
     round(logLik(frm(y ~ x + ar(t, gr = g1 * g2, cov = TRUE), d)), 6))
show("frmtmb ar(t, gr = factor(g), cov = TRUE)",
     round(logLik(frm(y ~ x + ar(t, gr = factor(g), cov = TRUE), d)), 6))
show("brms ar(t, gr = g1/g2)",
     names(brms::make_standata(y ~ x + ar(t, gr = g1/g2, cov = TRUE), d)))
show("brms ar(t, gr = g1:g2)",
     length(brms::make_standata(y ~ x + ar(t, gr = g1:g2, cov = TRUE), d)))
show("brms ar(t, gr = gf1:gf2) factors",
     length(brms::make_standata(y ~ x + ar(t, gr = gf1:gf2, cov = TRUE), d)))
show("brms ar(t, gr = interaction(g1, g2))",
     names(brms::make_standata(y ~ x + ar(t, gr = interaction(g1, g2),
                                          cov = TRUE), d)))
show("brms ar(t, gr = g1 + g2)",
     names(brms::make_standata(y ~ x + ar(t, gr = g1 + g2, cov = TRUE), d)))
show("brms ar(t, gr = g1 * g2)",
     names(brms::make_standata(y ~ x + ar(t, gr = g1 * g2, cov = TRUE), d)))
cat("  R's own reading of the two gr spellings\n")
show("class(d$g1 / d$g2)", c(class(d$g1 / d$g2), length(d$g1 / d$g2)))
show("eval(g1:g2) on numeric codes",
     utils::head(with(d, g1:g2), 12))
show("eval(gf1:gf2) on factors",
     utils::head(as.character(with(d, gf1:gf2)), 12))

cat("\n== D3 hypothesis() with no relation\n")
set.seed(20260917)
dh <- data.frame(y = rnorm(40), Age = rnorm(40), Trt = rnorm(40))
fh <- frm(y ~ Age + Trt, dh)
show("frmtmb hypothesis(fh, 'Age')", hypothesis(fh, "Age")$hypothesis)
show("frmtmb hypothesis(fh, 'Age = 0')",
     hypothesis(fh, "Age = 0")$hypothesis)
show("frmtmb hypothesis(fh, 'Age + Trt')",
     hypothesis(fh, "Age + Trt")$hypothesis)
show("frmtmb hypothesis(fh, 'Age x 0') loud",
     hypothesis(fh, "Age x 0"))
show("frmtmb hypothesis(fh, c('Age', 'Trt = 0'))",
     nrow(hypothesis(fh, c("Age", "Trt = 0"))$hypothesis))
be <- brms:::rename_pars(get("brmsfit_example1", envir = asNamespace("brms")))
show("brms hypothesis(ex1, 'Age')", brms::hypothesis(be, "Age"))
show("brms hypothesis(ex1, 'Age = 0')",
     brms::hypothesis(be, "Age = 0")$hypothesis)

cat("\n== D4 fit$data partial-matches fit$data2\n")
show("names(fh)", names(fh))
show("fh[['data', exact = TRUE]]", fh[["data", exact = TRUE]])
show("identical(fh$data, fh$data2)", identical(fh$data, fh$data2))
set.seed(20260917)
A <- diag(6)
A[A == 0] <- 0.3
dimnames(A) <- list(1:6, 1:6)
dd <- data.frame(g = factor(rep(1:6, each = 5)), x = rnorm(30))
dd$y <- dd$x + rnorm(6)[dd$g] + rnorm(30)
fa <- frm(y ~ x + (1 | gr(g, cov = A)), dd, data2 = list(A = A))
show("fit with data2: names(fa$data)", names(fa$data))
show("fit with data2: dim(fa$data$A)", dim(fa$data$A))
show("brms dim(ex1$data)", dim(be$data))
show("brms names(ex1$data)", names(be$data))
show("brms class(ex1$data)", class(be$data))
show("brms attr(ex1$data, 'terms') present",
     !is.null(attr(be$data, "terms")))
show("brms object.size(ex1$data) bytes",
     as.numeric(utils::object.size(be$data)))
show("frmtmb object.size(model.frame(fh)) bytes",
     as.numeric(utils::object.size(model.frame(fh))))
show("frmtmb object.size(fh) bytes", as.numeric(utils::object.size(fh)))

cat("\n== D5 frmtmb.latent::hmm_starts(1)\n")
show("hmm_starts(1)", frmtmb.latent::hmm_starts(1))
show("hmm_starts('a')", frmtmb.latent::hmm_starts("a"))
show("hmm_starts(NULL)", frmtmb.latent::hmm_starts(NULL))
show("hmm_starts(list())", frmtmb.latent::hmm_starts(list()))

cat("\n== D6 predict(newdata) without the grouping column\n")
set.seed(20260917)
dg <- data.frame(g = factor(rep(1:8, each = 5)), x = rnorm(40))
dg$y <- dg$x + rnorm(8)[dg$g] + rnorm(40)
fg <- frm(y ~ x + (1 | g), dg)
nd <- data.frame(x = c(0, 1))
nd_new <- data.frame(x = c(0, 1), g = factor(c("99", "98")))
show("frmtmb predict(allow_new_levels = TRUE)",
     predict(fg, newdata = nd, allow_new_levels = TRUE))
show("frmtmb predict(newdata with a NEW level)",
     predict(fg, newdata = nd_new, allow_new_levels = TRUE))
show("frmtmb predict(re_formula = NA)",
     predict(fg, newdata = nd, re_formula = NA))
show("frmtmb se.fit, filled column",
     predict(fg, newdata = nd, allow_new_levels = TRUE,
             se.fit = TRUE)$se.fit)
show("frmtmb se.fit, re_formula = NA",
     predict(fg, newdata = nd, re_formula = NA, se.fit = TRUE)$se.fit)
show("frmtmb fitted(allow_new_levels = TRUE)",
     fitted(fg, newdata = nd, allow_new_levels = TRUE))
show("frmtmb predict(no allow_new_levels)",
     predict(fg, newdata = nd))
# brms on the same question, on its own stored example fit. Fixture 1
# has `(1 + Trt | visit)`, so dropping `visit` from newdata is exactly
# the case above. Nothing compiles: the fit is brms's own.
b1 <- brms:::rename_pars(get("brmsfit_example1", envir = asNamespace("brms")))
b_full <- b1$data[1:3, ]
b_nogrp <- b_full[, setdiff(names(b_full), "visit"), drop = FALSE]
show("brms epred, full newdata", dim(brms::posterior_epred(
  b1, newdata = b_full, allow_new_levels = TRUE)))
show("brms epred, no visit column, allow_new_levels",
     dim(brms::posterior_epred(b1, newdata = b_nogrp,
                               allow_new_levels = TRUE)))
show("brms epred, no visit column, NO allow_new_levels",
     dim(brms::posterior_epred(b1, newdata = b_nogrp,
                               allow_new_levels = FALSE)))
show("brms validate_newdata fills visit with",
     brms:::validate_newdata(b_nogrp, b1, allow_new_levels = TRUE)$visit)

cat("\n== D7 log_lik() on a frequentist fit\n")
show("frmtmb log_lik(fh)", log_lik(fh))
show("frmtmb loo(fh)", loo(fh))
show("frmtmb waic(fh)", waic(fh))
show("exists log_lik generic",
     c(frmtmb = exists("log_lik", asNamespace("frmtmb")),
       sample = exists("log_lik", asNamespace("frmtmb.sample"))))

cat("\nDONE\n")
