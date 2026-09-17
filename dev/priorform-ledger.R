# brms 2.23.0's own assertions in this lane's scope, run LITERALLY
# against frmtmb, one at a time, and joined to the verdict recorded for
# each. The script stops if a verdict disagrees with the run: an
# assertion recorded as a pass that fails, or one recorded as a
# divergence that passes, is a stale record rather than a count.
#
#   Rscript dev/priorform-ledger.R lane|ref
#
# Sources: dev/brms-suite/brms/tests/testthat/tests.brmsformula.R (all
# 16), tests.priors.R (all 37), and the two tests.brm.R assertions for
# items 12 and 13. The code is brms's, with only `brm(` replaced by a
# frame-only frm() call and data loaded from brms's own .rda files.
# brms's suite runs under testthat edition 2, so this does too.
args <- commandArgs(trailingOnly = TRUE)
which_lib <- if (length(args)) args[1] else "lane"
lib <- switch(which_lib,
  ref = "C:/Users/adf44/source/r/rellib-r3",
  lane = "C:/Users/adf44/source/r/priorform-lib")
.libPaths(c(lib, "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(testthat))
suppressMessages(library(frmtmb))
# the option local_edition() sets, without a frame to defer a reset to
options(testthat.edition = 2)
load("dev/brms-suite/brms/data/epilepsy.rda")
load("dev/brms-suite/brms/data/inhaler.rda")
# brms's tests call brm() for these two; the refusal under test fires
# while the frame is built, so a frame-only fit reaches it
brm <- function(formula, data, ...) {
  frmtmb::frm(formula, data = data, dry_run = "frame", ...)
}

E <- list()
add <- function(id, file, block, verdict, code) {
  E[[length(E) + 1L]] <<- list(id = id, file = file, block = block,
                               verdict = verdict, code = code)
}
PASS <- "pass"

# ---- tests.brmsformula.R ------------------------------------------------
f <- "tests.brmsformula.R"
b <- "validates formulas of non-linear parameters"
add("F1", f, b, PASS, quote(
  expect_error(bf(y ~ a, ~ 1, a ~ 1), "Additional formulas must be named")))
add("F2", f, b, PASS, quote(
  expect_error(bf(y ~ a^x, a.b ~ 1), "not contain dots or underscores")))
add("F3", f, b, PASS, quote(
  expect_error(bf(y ~ a^(x+b), a_b ~ 1), "not contain dots or underscores")))
b <- "validates formulas of auxiliary parameters"
add("F4", f, b, PASS, quote(
  expect_error(bf(y ~ a, ~ 1, sigma ~ 1), "Additional formulas must be named")))
b <- "detects use if '~~'"
add("F5", f, b, PASS, quote(expect_error(bf(y~~x), "~~")))
b <- "does not change a 'brmsformula' object"
add("F6", f, b, PASS, quote({
  form <- bf(y ~ a, sigma ~ 1)
  expect_identical(form, bf(form))
}))
add("F7", f, b, PASS, quote({
  form <- bf(y ~ a, sigma ~ 1, a ~ x, nl = TRUE)
  expect_identical(form, bf(form))
}))
b <- "detects auxiliary parameter equations"
eq <- paste0(
  "cannot transfer: frmtmb has no equating of one dpar to ",
  "another (bf(y ~ x, sigma1 = \"sigma2\")); the call is ",
  "refused as an uninterpretable argument, a feature outside ",
  "this lane")
add("F8", f, b, eq, quote(expect_error(bf(y~x, sigma1 = "sigmaa2"),
  "Can only equate parameters of the same class")))
add("F9", f, b, eq, quote(expect_error(bf(y~x, mu3 = "mu2"),
  "Equating parameters of class 'mu' is not allowed")))
add("F10", f, b, eq, quote(expect_error(bf(y~x, sigma1 = "sigma1"),
  "Equating 'sigma1' with itself is not meaningful")))
add("F11", f, b, eq, quote(expect_error(bf(y~x, shape1 ~ x, shape2 = "shape1"),
  "Cannot use predicted parameters on the right-hand side")))
add("F12", f, b, eq, quote(expect_error(
  bf(y~x, shape1 = "shape3", shape2 = "shape1"),
  "Cannot use fixed parameters on the right-hand side")))
b <- "update_adterms works correctly"
ua <- paste0(
  "cannot transfer: frmtmb exports no update_adterms(), a ",
  "formula-editing helper outside this lane")
add("F13", f, b, ua, quote(expect_equal(
  update_adterms(y | trials(size) ~ x, ~ trials(10)), y | trials(10) ~ x)))
add("F14", f, b, ua, quote(expect_equal(
  update_adterms(y | trials(size) ~ x, ~ weights(w)),
  y | trials(size) + weights(w) ~ x)))
add("F15", f, b, ua, quote(expect_equal(
  update_adterms(y | trials(size) ~ x, ~ weights(w), action = "replace"),
  y | weights(w) ~ x)))
add("F16", f, b, ua, quote(expect_equal(
  update_adterms(y ~ x, ~ trials(10)), y | trials(10) ~ x)))

# ---- tests.priors.R -----------------------------------------------------
f <- "tests.priors.R"
b <- "default_prior finds all classes for which priors can be specified"
add("P1", f, b, paste0(
  "divergence: frmtmb lists 5 class theta rows (its internal ",
  "covariance parameters, a real class) and no per-coefficient ",
  "sd rows (its class sd addresses a block and refuses coef)"), quote(
  expect_equal(sort(default_prior(
    count ~ zBase * Trt + (1|patient) + (1+Trt|visit),
    data = epilepsy, family = "poisson")$class),
    sort(c(rep("b", 4), c("cor", "cor"), "Intercept", rep("sd", 6))))))
add("P2", f, b, paste0(
  "cannot transfer: frmtmb's sratio() takes no threshold = ",
  "\"equidistant\" and has no cse() alias (family surface, lane",
  " wt-famlink)"), quote(
  expect_equal(sort(default_prior(rating ~ treat + period + cse(carry),
    data = inhaler, family = sratio(threshold = "equidistant"))$class),
    sort(c(rep("b", 4), "delta", rep("Intercept", 1))))))
b <- "set_prior allows arguments to be vectors"
add("P3", f, b, paste0(
  "divergence: the object is a frmtmb_priorlist, which holds ",
  "parsed densities; $ reads brms's columns from it"), quote({
  bprior <- set_prior("normal(0, 2)", class = c("b", "sd"))
  expect_is(bprior, "brmsprior")
}))
add("P4", f, b, PASS, quote({
  bprior <- set_prior("normal(0, 2)", class = c("b", "sd"))
  expect_equal(bprior$prior, rep("normal(0, 2)", 2))
}))
add("P5", f, b, PASS, quote({
  bprior <- set_prior("normal(0, 2)", class = c("b", "sd"))
  expect_equal(bprior$class, c("b", "sd"))
}))
b <- "print for class brmsprior works correctly"
add("P6", f, b, PASS, quote(expect_output(print(set_prior("normal(0,1)")),
  fixed = TRUE, "b ~ normal(0,1)")))
add("P7", f, b, PASS, quote(expect_output(
  print(set_prior("normal(0,1)", coef = "x")),
  "b_x ~ normal(0,1)", fixed = TRUE)))
add("P8", f, b, PASS, quote(expect_output(
  print(set_prior("cauchy(0,1)", class = "sd", group = "x")),
  "sd_x ~ cauchy(0,1)", fixed = TRUE)))
add("P9", f, b, paste0(
  "cannot transfer: set_prior(check = FALSE) passes Stan code ",
  "through, and frmtmb builds no Stan program"), quote(
  expect_output(
    print(set_prior("target += normal_lpdf(x | 0,1))", check = FALSE)),
    "target += normal_lpdf(x | 0,1))", fixed = TRUE)))
b <- "default_prior returns correct nlpar names for random effects pars"
add("P10", f, b, PASS, quote({
  set.seed(1)
  dat <- data.frame(y = rnorm(10), x = rnorm(10), g = rep(1:2, 5))
  bform <- bf(y ~ a - b^x, a + b ~ (1+x|g), nl = TRUE)
  gp <- default_prior(bform, data = dat)
  expect_equal(sort(unique(gp$nlpar)), c("", "a", "b"))
}))
b <- "default_prior returns correct fixed effect names for GAMMs"
add("P11", f, b, paste0(
  "divergence: a smooth's unpenalized column is named s(x).fx1 ",
  "where brms names it sx_1 (coefficient naming, lane wt-",
  "brmsnames)"), quote({
  set.seed(1)
  dat <- data.frame(y = rnorm(10), x = rnorm(10), z = rnorm(10),
                    g = rep(1:2, 5))
  prior <- default_prior(y ~ z + s(x) + (1|g), data = dat)
  expect_equal(prior[prior$class == "b", ]$coef, c("", "sx_1", "z"))
}))
add("P12", f, b, paste0(
  "divergence: as P11, and the nonlinear intercept is listed as",
  " (Intercept) where brms writes Intercept; set_prior() takes ",
  "both"), quote({
  set.seed(1)
  dat <- data.frame(y = rnorm(10), x = rnorm(10), z = rnorm(10),
                    g = rep(1:2, 5))
  prior <- default_prior(bf(y ~ lp, lp ~ z + s(x) + (1|g), nl = TRUE),
                         data = dat)
  expect_equal(prior[prior$class == "b", ]$coef,
               c("", "Intercept", "sx_1", "z"))
}))
b <- "default_prior returns correct prior names for auxiliary parameters"
add("P13", f, b, paste0(
  "divergence: frmtmb lists the class-wide sd row once with no ",
  "dpar and no per-coefficient sd row, so the phi rows are b, b",
  " z, Intercept and sd g"), quote({
  set.seed(1)
  dat <- data.frame(y = rnorm(10), x = rnorm(10), z = rnorm(10),
                    g = rep(1:2, 5))
  bform <- bf(y ~ 1, phi ~ z + (1|g), family = Beta())
  dat$y <- stats::runif(10)
  prior <- default_prior(bform, data = dat)
  prior <- prior[prior$dpar == "phi", ]
  pdata <- data.frame(class = c("b", "b", "Intercept", rep("sd", 3)),
                      coef = c("", "z", "", "", "", "Intercept"),
                      group = c(rep("", 4), "g", "g"),
                      stringsAsFactors = FALSE)
  pdata <- pdata[with(pdata, order(class, group, coef)), ]
  expect_equivalent(prior[, c("class", "coef", "group")], pdata)
}))
b <- "default_prior returns correct priors for multivariate models"
mvdat <- quote({
  set.seed(1)
  dat <- data.frame(y1 = rnorm(10), y2 = c(1, rep(1:3, 3)),
                    x = rnorm(10), g = rep(1:2, 5))
})
add("P14", f, b, PASS, bquote({
  .(mvdat)
  bform <- bf(mvbind(y1, y2) ~ x + (x|ID1|g)) + set_rescor(TRUE)
  prior <- default_prior(bform, dat, family = gaussian())
  expect_equal(prior[prior$resp == "y1" & prior$class == "b", "coef"],
               c("", "x"))
}))
add("P15", f, b, paste0(
  "divergence: default_prior() describes frm() unless route = ",
  "\"sample\", and frm() is flat in every slot. With frmtmb.sample ",
  "attached, route = \"sample\" ALSO gives (flat) here where brms ",
  "gives lkj(1): a filed frmtmb.sample defect, not fixed in this lane ",
  "(dev/priorform-findings.md section 6)"), bquote({
  .(mvdat)
  bform <- bf(mvbind(y1, y2) ~ x + (x|ID1|g)) + set_rescor(TRUE)
  prior <- default_prior(bform, dat, family = gaussian())
  expect_equal(prior[prior$class == "rescor", "prior"], "lkj(1)")
}))
mvfam <- paste0(
  "cannot transfer: frmtmb takes no list of families as ",
  "`family` (a multivariate model attaches one family per ",
  "bf()), family surface outside this lane")
for (k in 1:4) {
  expr <- switch(k,
    quote(expect_true(any(with(prior, class == "sigma" & resp == "y1")))),
    quote(expect_true(any(with(prior, class == "ar" & resp == "y1")))),
    quote(expect_true(any(with(prior, class == "phi" & resp == "y2")))),
    quote(expect_true(!any(with(prior, class == "ar" & resp == "y2")))))
  add(paste0("P", 15 + k), f, b, mvfam, bquote({
    .(mvdat)
    family <- list(gaussian, Beta())
    bform <- bf(y1 ~ x + (x|ID1|g) + ar()) + bf(y2 ~ 1)
    prior <- default_prior(bform, dat, family = family)
    .(expr)
  }))
}
b <- "default_prior returns correct priors for categorical models"
add("P20", f, b, paste0(
  "cannot transfer: categorical() refuses this numeric response",
  " ('fewer than two categories'), response-type validation ",
  "owned by lane wt-famlink"), quote({
  set.seed(1)
  dat <- data.frame(y2 = c(1, rep(1:3, 3)), x = rnorm(10),
                    g = rep(1:2, 5))
  prior <- default_prior(y2 ~ x + (x|ID1|g), data = dat,
                         family = categorical())
  expect_equal(prior[prior$dpar == "mu2" & prior$class == "b", "coef"],
               c("", "x"))
}))
b <- "set_prior alias functions produce equivalent results"
add("P21", f, b, PASS, quote(expect_equal(
  set_prior("normal(0, 1)", class = "sd"),
  prior(normal(0, 1), class = sd))))
add("P22", f, b, PASS, quote(expect_equal(
  set_prior("normal(0, 1)", class = "sd", nlpar = "a"),
  prior(normal(0, 1), class = "sd", nlpar = a))))
add("P23", f, b, PASS, quote(expect_equal(
  set_prior("normal(0, 1)", class = "sd", nlpar = "a"),
  prior_(~normal(0, 1), class = ~sd, nlpar = quote(a)))))
add("P24", f, b, PASS, quote(expect_equal(
  set_prior("normal(0, 1)", class = "sd"),
  prior_string("normal(0, 1)", class = "sd"))))
b <- "external interface of validate_prior works correctly"
vp <- quote({
  prior1 <- prior(normal(0,10), class = b) + prior(cauchy(0,2), class = sd)
  prior1 <- validate_prior(prior1, count ~ zAge + zBase * Trt + (1|patient),
                           data = epilepsy, family = poisson())
})
add("P25", f, b, PASS, bquote({
  .(vp)
  expect_true(all(c("b", "Intercept", "sd") %in% prior1$class))
}))
add("P26", f, b, paste0(
  "divergence: 10 rows against brms's 9, for the reasons in P1:",
  " two class theta rows added, one per-coefficient sd row ",
  "absent"), bquote({
  .(vp)
  expect_equal(nrow(prior1), 9)
}))
b <- "overall intercept priors are adjusted for the intercept"
add("P27", f, b, paste0(
  "divergence: the fit route is flat. With frmtmb.sample attached, ",
  "route = \"sample\" gives student_t(3, 2, 2.5) where brms gives ",
  "student_t(3, -8, 2.5): the sampling default does not subtract the ",
  "offset, a filed frmtmb.sample defect not fixed in this lane ",
  "(dev/priorform-findings.md section 6)"), quote({
  dat <- data.frame(y = rep(c(1, 3), each = 5), off = 10)
  prior1 <- default_prior(y ~ 1 + offset(off), dat)
  int_prior <- prior1$prior[prior1$class == "Intercept"]
  expect_equal(int_prior, "student_t(3, -8, 2.5)")
}))
b <- "as.brmsprior works correctly"
abp <- quote({
  dat <- data.frame(prior = "normal(0,1)", x = "test", coef = c("a", "b"))
  bprior <- as.brmsprior(dat)
})
add("P28", f, b, PASS, bquote({
  .(abp); expect_equal(bprior$prior, rep("normal(0,1)", 2)) }))
add("P29", f, b, PASS, bquote({
  .(abp); expect_equal(bprior$class, rep("b", 2)) }))
add("P30", f, b, PASS, bquote({
  .(abp); expect_equal(bprior$coef, c("a", "b")) }))
add("P31", f, b, PASS, bquote({
  .(abp); expect_equal(bprior$x, NULL) }))
add("P32", f, b, PASS, bquote({
  .(abp); expect_equal(bprior$lb, rep(NA_character_, 2)) }))
b <- "prior tags are correctly applied"
tg <- paste0(
  "cannot transfer: set_prior() takes no tag, and a brms table ",
  "carrying one is refused by name, because a tag names a prior",
  " inside a Stan program")
tagged <- quote({
  prior1 <- prior(normal(0, 1), class = sd, tag = "prior_tag1")
  prior2 <- prior(normal(0, 5), class = b, tag = "prior_tag2")
  prior3 <- prior(normal(0, 0.5), coef = "Trt1", tag = "prior_tag3")
  prior4 <- prior(normal(0, 10), class = "Intercept", tag = "prior_tag4")
  prior5 <- prior(lkj_corr_cholesky(3), class = "L", group = "visit",
                  tag = "prior_tag5")
  v <- validate_prior(
    c(prior1, prior2, prior3, prior4, prior5),
    formula = count ~ zBase * Trt + (1 | patient) + (1 + Trt | visit),
    data = epilepsy, family = poisson())
})
tag_expect <- list(
  quote(expect_equal(v[which(v$class == "sd"),]$tag[[1]], "prior_tag1")),
  quote(expect_equal(v[which(v$class == "b" & v$coef != "Trt1"),]$tag[[1]],
                     "prior_tag2")),
  quote(expect_equal(v[which(v$class == "b" & v$coef == "Trt1"),]$tag,
                     "prior_tag3")),
  quote(expect_equal(v[which(v$class == "Intercept"),]$tag, "prior_tag4")),
  quote(expect_equal(v[which(v$class == "L"),]$tag[[2]], "prior_tag5")))
for (k in 1:5) {
  add(paste0("P", 32 + k), f, b, tg, bquote({
    .(tagged)
    .(tag_expect[[k]])
  }))
}

# ---- tests.brm.R ----------------------------------------------------------
f <- "tests.brm.R"
b <- "brm produces expected errors"
add("B1", f, b, PASS, quote({
  set.seed(1)
  dat <- data.frame(y = rnorm(10), x = rnorm(10), g = rep(1:5, 2))
  expect_error(brm(y ~ x + (1|g) + (x|g), dat),
               "Duplicated group-level effects are not allowed")
}))
add("B2", f, b, PASS, quote({
  set.seed(1)
  dat <- data.frame(y = rnorm(10), x = rnorm(10), g = rep(1:5, 2))
  expect_error(brm(y~x*cs(g), dat), fixed = TRUE,
               "The term 'x:cs(g)' is invalid")
}))

# an assertion passes when its code runs to the end without a failed
# expectation or an error; testthat signals a failure as a condition
run_one <- function(e) {
  tryCatch({
    eval(e$code, new.env(parent = globalenv()))
    "pass"
  }, expectation_failure = function(c) "fail",
     error = function(c) "fail")
}

rows <- lapply(E, function(e) {
  got <- run_one(e)
  kind <- sub(":.*", "", e$verdict)
  data.frame(id = e$id, file = e$file, block = e$block, run = got,
             verdict = kind, reason = sub("^[^:]*: ?", "", e$verdict),
             stringsAsFactors = FALSE)
})
tab <- do.call(rbind, rows)
tab$reason[tab$verdict == "pass"] <- ""
utils::write.table(tab, paste0("dev/priorform-ledger-", which_lib, ".tsv"),
                   sep = "\t", row.names = FALSE, quote = FALSE)

cat("frmtmb from", find.package("frmtmb"), "\n")
stale <- tab[(tab$verdict == "pass") != (tab$run == "pass"), ]
if (which_lib == "lane" && nrow(stale)) {
  print(stale[, c("id", "run", "verdict", "reason")], row.names = FALSE)
  stop(nrow(stale), " ledger verdicts disagree with the run")
}
cat("\n<!-- priorform-ledger-", which_lib, ":begin -->\n", sep = "")
cat(sprintf("%d brms assertions run against frmtmb (%s): %d pass, %d fail\n",
            nrow(tab), which_lib, sum(tab$run == "pass"),
            sum(tab$run == "fail")))
if (which_lib == "lane") {
  cat(sprintf(paste0("verdicts: %d pass, %d deliberate divergence, ",
                     "%d cannot transfer\n\n"),
              sum(tab$verdict == "pass"), sum(tab$verdict == "divergence"),
              sum(tab$verdict == "cannot transfer")))
  cat("| id | brms file | block | run | verdict | reason |\n")
  cat("|---|---|---|---|---|---|\n")
  for (i in seq_len(nrow(tab))) {
    cat(sprintf("| %s | %s | %s | %s | %s | %s |\n", tab$id[i], tab$file[i],
                tab$block[i], tab$run[i], tab$verdict[i], tab$reason[i]))
  }
} else {
  cat("ids passing on this build:",
      paste(tab$id[tab$run == "pass"], collapse = " "), "\n")
}
cat("<!-- priorform-ledger-", which_lib, ":end -->\n", sep = "")
