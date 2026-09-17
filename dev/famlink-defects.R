# Lane wt-famlink: every behavioral check this lane changes, run against
# ONE frmtmb library chosen on the command line, so the same script shows
# the base failing and the lane passing.
#
#   Rscript dev/famlink-defects.R base > dev/famlink-defects-base-log.txt
#   Rscript dev/famlink-defects.R lane > dev/famlink-defects-lane-log.txt
#
# Each check states what brms 2.23.0 does (dev/famlink-brms-behavior-log.txt
# is the record of that) and reports PASS when frmtmb does the same.
# Seeds are fixed per construction and printed beside it.
arm <- commandArgs(trailingOnly = TRUE)[1]
lib <- switch(arm,
              base = "C:/Users/adf44/source/r/rellib-r3",
              lane = "C:/Users/adf44/source/r/famlink-lib",
              stop("arm must be base or lane"))
.libPaths(c(lib, "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("arm:", arm, " frmtmb", format(packageVersion("frmtmb")), "from",
    dirname(find.package("frmtmb")), "\n\n")

results <- data.frame(item = character(0), check = character(0),
                      pass = logical(0), got = character(0))

# Runs `expr` capturing its value, error, warnings and messages.
run <- function(expr) {
  msgs <- character(0)
  warns <- character(0)
  val <- NULL
  err <- NULL
  withCallingHandlers(
    tryCatch(val <- force(expr), error = function(e) err <<- e),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    },
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  list(value = val, error = err, messages = msgs, warnings = warns)
}

record <- function(item, check, pass, got) {
  got <- paste(gsub("\n", " ", got), collapse = " | ")
  results[nrow(results) + 1L, ] <<- list(item, check, isTRUE(pass),
                                         substr(got, 1, 160))
  cat(sprintf("[%s] item %s: %s\n    got: %s\n",
              if (isTRUE(pass)) "PASS" else "FAIL", item, check,
              substr(got, 1, 300)))
}

describe <- function(r) {
  if (!is.null(r$error)) return(paste("ERROR:", conditionMessage(r$error)))
  v <- r$value
  out <- if (is.character(v) && length(v) <= 3) paste0("\"", v, "\"", collapse = ", ")
         else paste("value of class", paste(class(v), collapse = "/"))
  if (length(r$warnings)) out <- paste(out, "; WARNING:", r$warnings)
  if (length(r$messages)) out <- paste(out, "; MESSAGE:", r$messages)
  out
}

expect_value <- function(item, label, expr, want) {
  r <- run(expr)
  record(item, paste(label, "is", deparse(want)),
         is.null(r$error) && identical(r$value, want), describe(r))
}

expect_refusal <- function(item, label, expr, pattern) {
  r <- run(expr)
  record(item, paste(label, "is refused matching", deparse(pattern)),
         !is.null(r$error) &&
           grepl(pattern, conditionMessage(r$error), fixed = TRUE),
         describe(r))
}

# ---- item 1: brms's fields on the family object -------------------------
expect_value(1, "student()$link", student()$link, "identity")
expect_value(1, "beta_binomial()$link_phi", beta_binomial()$link_phi, "log")
expect_value(1, "zero_inflated_binomial()$link_zi",
             zero_inflated_binomial()$link_zi, "logit")
expect_value(1, "student()$link_nu", student()$link_nu, "logm1")
expect_value(1, "cumulative('probit')$link", cumulative("probit")$link,
             "probit")
# the fields are computed on read, not stored, so this checks a read
expect_value(1, "is.function(student()$linkinv)",
             is.function(student()$linkinv), TRUE)
expect_value(1, "is.function(stats::gaussian()-converted $linkinv)",
             is.function(frmtmb:::as_frmtmb_family(stats::poisson())$linkinv),
             TRUE)
r <- run(student()$lin)
record(1, "student()$lin (a partial name) does not silently return a field",
       !is.null(r$error) || is.null(r$value), describe(r))

# ---- item 2: an unquoted link -------------------------------------------
expect_value(2, "student(identity)$link", student(identity)$link, "identity")
expect_value(2, "negbinomial(sqrt)$link", negbinomial(sqrt)$link, "sqrt")
expect_value(2, "exponential(log)$link", exponential(log)$link, "log")
expect_value(2, "geometric(identity)$link", geometric(identity)$link,
             "identity")
expect_value(2, "zero_inflated_poisson(log)$link",
             zero_inflated_poisson(log)$link, "log")
lk <- "probit"
expect_value(2, "bernoulli(lk) with lk <- 'probit'", bernoulli(lk)$link,
             "probit")

# ---- item 3: the mean link is validated per family -----------------------
expect_refusal(3, "bernoulli('sqrt')", bernoulli("sqrt"),
               "'sqrt' is not a supported link for family 'bernoulli'")
expect_refusal(3, "exponential('cloglog')", exponential("cloglog"),
               "'cloglog' is not a supported link for family 'exponential'")
expect_refusal(3, "Beta('1/mu^2')", Beta("1/mu^2"),
               "'1/mu^2' is not a supported link for family 'beta'")
expect_refusal(3, "zero_inflated_negbinomial('logit')",
               zero_inflated_negbinomial("logit"),
               "is not a supported link for family 'zero_inflated_negbinomial'")
expect_refusal(3, "beta_binomial('log')", beta_binomial("log"),
               "'log' is not a supported link for family 'beta_binomial'")
expect_refusal(3, "brmsfamily('poisson', link = 'inverse')",
               brmsfamily("poisson", link = "inverse"),
               "'inverse' is not a supported link for family 'poisson'")
expect_refusal(3, "stats::poisson('inverse') through frm()",
               frm(y ~ 1, data.frame(y = rpois(20, 3)),
                   family = stats::poisson("inverse")),
               "'inverse' is not a supported link for family 'poisson'")
# the silent wrong answer from dev/brmssuite-defects.R D2, same seed
set.seed(1)
d <- data.frame(y = rbinom(60, 1, 0.4), x = rnorm(60))
r <- run(frm(y ~ x, d, family = bernoulli("sqrt")))
record(3, "frm(y ~ x, bernoulli('sqrt')) is refused, seed 1, n = 60",
       !is.null(r$error),
       if (is.null(r$error)) paste("FITTED: mu coef",
                                   paste(format(unlist(fixef(r$value)), digits = 8),
                                         collapse = " "), "; warnings:",
                                   paste(r$warnings, collapse = " / "))
       else paste("ERROR:", conditionMessage(r$error)))

# ---- item 5: mixture support --------------------------------------------
expect_refusal(5, "mixture(poisson(), categorical())",
               mixture(poisson(), categorical()),
               "Some of the families are not allowed in mixture models")
expect_refusal(5, "mixture(lognormal, exgaussian, poisson())",
               mixture(lognormal, exgaussian, poisson()),
               "Cannot mix families with real and integer support")
expect_refusal(5, "mixture(poisson, 'cumulative')",
               mixture(poisson, "cumulative"),
               "Cannot mix ordinal and non-ordinal families")

# ---- item 11: an unordered factor response to an ordinal family ---------
set.seed(11)
dd <- data.frame(y = factor(sample(c("high", "low", "mid"), 60, TRUE)),
                 x = rnorm(60))
for (fn in c("cumulative", "sratio", "cratio", "acat")) {
  expect_refusal(11, paste0("frm(unordered factor ~ x, ", fn, "()), seed 11"),
                 frm(y ~ x, dd, family = get(fn)()),
                 paste0("Family '", fn, "' requires either positive integers ",
                        "or ordered factors"))
}

# ---- item 15: the bernoulli suggestion -----------------------------------
set.seed(15)
db <- data.frame(y = rbinom(40, 1, 0.5), x = rnorm(40),
                 y2 = sample(1:2, 40, TRUE), y3 = sample(1:3, 40, TRUE),
                 f2 = factor(sample(c("a", "b"), 40, TRUE)))
bern <- "Only 2 levels detected so that family 'bernoulli' might be a more efficient choice"
msg_check <- function(label, expr, want) {
  r <- run(expr)
  hit <- any(grepl(bern, r$messages, fixed = TRUE))
  record(15, paste(label, if (want) "messages" else "does not message",
                   "the bernoulli suggestion"),
         is.null(r$error) && hit == want, describe(r))
}
msg_check("binomial y | trials(1), seed 15",
          frm(y | trials(1) ~ x, db, family = binomial()), TRUE)
msg_check("beta_binomial y | trials(1), seed 15",
          frm(y | trials(1) ~ x, db, family = beta_binomial()), TRUE)
msg_check("cumulative on 2 categories, seed 15",
          frm(y2 ~ x, db, family = cumulative()), TRUE)
msg_check("categorical on 2 levels, seed 15",
          frm(f2 ~ x, db, family = categorical()), TRUE)
msg_check("cumulative on 3 categories, seed 15",
          frm(y3 ~ x, db, family = cumulative()), FALSE)

# ---- item 16: brmsfamily() and c(family, link) ---------------------------
expect_value(16, "brmsfamily('gaussian', inverse)$link",
             brmsfamily("gaussian", inverse)$link, "inverse")
expect_value(16, "brmsfamily('geometric', 'identity')$family",
             brmsfamily("geometric", "identity")$family, "geometric")
expect_value(16, "brmsfamily('zi_poisson')$link_zi",
             brmsfamily("zi_poisson")$link_zi, "logit")
expect_value(16, "frmtmb:::as_frmtmb_family(c('weibull', 'log'))$link",
             frmtmb:::as_frmtmb_family(c("weibull", "log"))$link, "log")
expect_refusal(16, "frmtmb:::as_frmtmb_family(c('weibull', 'sqrt'))",
               frmtmb:::as_frmtmb_family(c("weibull", "sqrt")),
               "'sqrt' is not a supported link for family 'weibull'")
set.seed(16)
dw <- data.frame(y = rweibull(50, 2, 1), x = rnorm(50))
r <- run(frm(y ~ x, dw, family = c("weibull", "log")))
record(16, "frm(family = c('weibull', 'log')) fits, seed 16",
       is.null(r$error) && inherits(r$value, "frmtmb_fit"), describe(r))
expect_refusal(16, "frm(family = c('categorical', 'probit'))",
               frm(y ~ x, dw, family = c("categorical", "probit")),
               "'probit' is not a supported link for family 'categorical'")

cat("\n==== summary (", arm, ") ====\n", sep = "")
by_item <- split(results$pass, as.numeric(results$item))
cat(sprintf("item %-3s pass %d of %d\n", names(by_item),
            vapply(by_item, sum, 0L), lengths(by_item)), sep = "")
cat("TOTAL PASS", sum(results$pass), "of", nrow(results), "\n")
write.table(results, paste0("dev/famlink-defects-", arm, ".tsv"), sep = "\t",
            row.names = FALSE, quote = FALSE)
