# frmtmb and brms in the same session.
#
# This is the one script in dev/brms-vignettes/ that attaches brms on
# purpose. Every other script reaches brms with the `brms::` prefix so
# that a post-processing name must resolve to frmtmb or fail. Here the
# question is the opposite one: a user who already has brms attached
# installs frmtmb, and needs to know what changes.
#
# Three things get measured.
#   1. Which names collide, and which package wins in which attach
#      order.
#   2. Whether the brms generics still dispatch to frmtmb's methods when
#      brms is the one that supplies the generic.
#   3. Whether a plain frm() call still works with brms on the path.
#
# Read the comments, not just the code: every `<EDGE>` marker on a bv()
# call is a place where the collision costs the user something.

source(file.path(Sys.getenv(
  "BV_DIR",
  unset = "C:/Users/adf44/source/r/frmtmb-wt-brmsvig/dev/brms-vignettes"
), "_harness.R"))
bv_load()
bv_init("brms-coexistence")

options(mc.cores = 1)
set.seed(1234)

## ---- data -----------------------------------------------------------

bv("data", "two tiny data sets", {
  d <- data.frame(x = rnorm(80), g = gl(8, 10))
  d$y <- 1 + 2 * d$x + rnorm(80)
  d$z <- rbinom(80, 1, plogis(0.5 * d$x))
  assign("d", d, envir = globalenv())
  utils::str(d)
}, NA_character_, "")

## ---- attaching brms after frmtmb -------------------------------------

bv("data", "library(brms) with frmtmb already attached", {
  msg <- utils::capture.output(library(brms), type = "message")
  cat(paste(msg, collapse = "\n"), "\n")
  assign("attach_msg", msg, envir = globalenv())
  "attached"
}, NA_character_, "")

# The collision is large and entirely one-directional in this order:
# brms sits above frmtmb on the search path, so every shared name
# resolves to brms.
bv("post", "ML: the masking table (brms attached second)", {
  fe <- getNamespaceExports("frmtmb")
  be <- getNamespaceExports("brms")
  both <- sort(intersect(fe, be))
  # find() lists every attached position a name lives at, in search
  # order, so its first entry is the one an unqualified call reaches.
  wins <- vapply(both, function(n) find(n)[1L], "")
  cat("frmtmb exports:", length(fe), "  brms exports:", length(be), "\n")
  cat("names exported by BOTH:", length(both), "\n")
  cat("resolving to brms:", sum(wins == "package:brms"), "\n")
  cat("resolving to frmtmb:", sum(wins == "package:frmtmb"), "\n")
  cat("resolving elsewhere:", sum(!wins %in%
        c("package:brms", "package:frmtmb")), "\n")
  cat("\nthe shared names:\n")
  print(both)
  length(both)
}, "SPELLING",
"92 of frmtmb's 134 exports are also brms exports, and attaching brms second masks all 92; the attach message lists them but does not say what to do")

# The exact message a user sees, quoted, because it is the only warning
# they get.
bv("post", "ML: the exact masking message", {
  cat(paste(attach_msg[grepl("masked", attach_msg)], collapse = "\n"), "\n")
  cat("\nfirst line of the object list:\n")
  i <- grep("masked from 'package:frmtmb'", attach_msg)
  cat(attach_msg[i + 1L], "\n")
  "quoted"
}, NA_character_, "")

## ---- do the brms generics still dispatch on a frmtmb fit? ------------

# The model calls themselves. frm() is a frmtmb-only name so it is never
# masked, but bf() is masked, which is the first thing a copied line
# hits.
fit_g <- bv("model", "ML: plain frm() with brms attached", {
  frm(y ~ x + (1 | g), data = d, family = gaussian())
}, NA_character_, "")

fit_b <- bv("model", "ML: frm() binomial with brms attached", {
  frm(z ~ x + (1 | g), data = d, family = binomial())
}, NA_character_, "")

# bf() is the highest-traffic collision, and frm() catches it: the
# formula object is a brmsformula and the refusal says so and names the
# fix.
bv("model", "ML: frm(bf(...)) where bf() is brms's", {
  frm(bf(y ~ x + (1 | g)), data = d, family = gaussian())
}, "REFUSAL",
"the message names frmtmb::bf() and the reverse attach order, and it is the right pointer; without it the failure would be an unreadable structure error")

bv("model", "ML: frm(frmtmb::bf(...)) [the fix]", {
  frm(frmtmb::bf(y ~ x + (1 | g)), data = d, family = gaussian())
}, "SPELLING", "the bf() call needs the frmtmb:: prefix once brms is attached")

# The generics. Each of these is a name brms owns, and each is supposed
# to reach frmtmb's method by S3 dispatch even though the generic is
# brms's. Both a gaussian and a binomial fit are exercised, because a
# dispatch failure could be family dependent.
generic_row <- function(nm, fit, what, edge = NA_character_, why = "") {
  bv("post", paste0("ML: ", nm, "(", what, ") with brms attached"), {
    grDevices::pdf(NULL); on.exit(grDevices::dev.off())
    fn <- get(nm)
    val <- if (identical(nm, "hypothesis")) fn(fit, "x") else fn(fit)
    stopifnot(!inherits(val, "brmsfit"))
    cat("dispatched to a frmtmb method; class:",
        paste(class(val), collapse = ", "), "\n")
    class(val)[1]
  }, edge, why)
}

for (nm in c("conditional_effects", "hypothesis", "pp_check", "ngrps",
             "fixef", "ranef", "VarCorr")) {
  generic_row(nm, fit_g, "gaussian")
  generic_row(nm, fit_b, "binomial")
}

# loo() is the one generic that does not reach a frmtmb method, and the
# reason is not the collision: frmtmb registers loo only on
# frmtmb_draws. So brms's generic is reached and finds nothing.
generic_row("loo", fit_g, "gaussian", "MISSING",
            "loo has no frmtmb_fit method with or without brms attached; the error is a bare dispatch failure that names neither frm_sample() nor AIC()")
generic_row("loo", fit_b, "binomial", "MISSING",
            "same as the gaussian fit")

## ---- the collisions that are NOT generics ---------------------------
#
# A masked generic still finds frmtmb's method. A masked CONSTRUCTOR
# does not: it just builds the wrong kind of object, and the failure
# surfaces later, or not at all.

bv("post", "ML: set_prior() resolves to brms", {
  p <- set_prior("normal(0, 1)", class = "b")
  cat("class:", paste(class(p), collapse = ", "), "\n")
  stopifnot(inherits(p, "brmsprior"))
  "brmsprior"
}, "SPELLING",
"the two set_prior() signatures agree closely enough that the call succeeds and returns a brmsprior; frm(priors =) then has to refuse it, so the collision is only found one step later. Prefix with frmtmb::")

bv("post", "ML: frmtmb::set_prior() [the fix]", {
  p <- frmtmb::set_prior("normal(0, 1)", class = "b")
  cat("class:", paste(class(p), collapse = ", "), "\n")
  p
}, "SPELLING", "the frmtmb:: prefix is needed for every masked constructor")

bv("post", "ML: get_prior() resolves to brms", {
  p <- get_prior(y ~ x, data = d)
  stopifnot(inherits(p, "brmsprior"))
  cat("class:", paste(class(p), collapse = ", "), "\n")
  "brmsprior"
}, "SPELLING", "brms's get_prior() answers, so a user reads brms's default priors and believes they describe the frmtmb fit; nothing warns")

bv("post", "ML: frmtmb::get_prior() [the fix]", {
  p <- frmtmb::get_prior(y ~ x, data = d, family = gaussian())
  print(utils::head(p))
  "data.frame"
}, "SPELLING", "the frmtmb:: prefix again; the two return different classes so the confusion is silent")

bv("post", "ML: prior() is brms-only", {
  p <- prior(normal(0, 1), class = "b")
  stopifnot(inherits(p, "brmsprior"))
  cat("prior() is not a frmtmb export, so it is not masked: brms owns it\n")
  "brmsprior"
}, "MISSING",
"frmtmb has no prior(); the brms idiom prior(normal(0, 1), class = b) resolves to brms and produces an object frm() cannot use, with no collision warning at all because there is nothing to collide with")

bv("post", "ML: mixture() resolves to brms", {
  m <- suppressMessages(mixture(gaussian(), gaussian()))
  stopifnot(inherits(m, "brmsfamily"))
  cat("class:", paste(class(m), collapse = ", "), "\n")
  "brmsfamily"
}, "SPELLING", "a brms mixfamily comes back where a frmtmb_family was wanted; frmtmb::mixture() is the fix")

bv("post", "ML: custom_family() resolves to brms", {
  cf <- custom_family("f", dpars = "mu", links = "identity", type = "real")
  stopifnot(inherits(cf, "customfamily"))
  cat("class:", paste(class(cf), collapse = ", "), "\n")
  "customfamily"
}, "SPELLING",
"the two custom_family() signatures are different enough that brms's accepts arguments frmtmb's would reject, so the wrong call succeeds; frmtmb::custom_family() is the fix")

bv("post", "ML: lf() resolves to brms", {
  l <- lf(sigma ~ x)
  cat("class:", paste(class(l), collapse = ", "), "\n")
  stopifnot(!inherits(l, "frmtmb_lf"))
  "brms lf"
}, "SPELLING", "a brms lf object; frmtmb::lf() returns frmtmb_lf")

bv("post", "ML: nlf() resolves to brms", {
  l <- nlf(s ~ a * x)
  cat("class:", paste(class(l), collapse = ", "), "\n")
  stopifnot(!inherits(l, "frmtmb_nlf"))
  "brms nlf"
}, "SPELLING", "same as lf()")

# The families are the quietest collision of the lot: brms's
# constructors return a brmsfamily, and frm() has to recognize it.
bv("post", "ML: family constructors resolve to brms", {
  f <- student()
  cat("student() class:", paste(class(f), collapse = ", "), "\n")
  fit <- frm(y ~ x, data = d, family = f)
  cat("frm() accepted a brmsfamily:", family(fit)$family, "\n")
  "accepted"
}, "BEHAVIOR",
"29 family constructors are masked, so family = student() hands frm() a brmsfamily rather than a frmtmb_family; frm() converts it silently, which is the good outcome but means a link or dpar difference between the two packages would never be reported")

## ---- the reverse attach order ---------------------------------------
#
# The same questions with frmtmb attached second, which needs a fresh R
# process because the search path order is fixed at attach time. The
# child prints one line per check and this row relays them.

bv("post", "ML: reversed attach order (brms first, then frmtmb)", {
  child <- tempfile(fileext = ".R")
  writeLines(c(
    'suppressMessages(library(brms))',
    'suppressMessages(pkgload::load_all(',
    '  Sys.getenv("BV_PKG", unset = "C:/Users/adf44/source/r/frmtmb-wt-brmsvig"),',
    '  quiet = TRUE, export_all = FALSE))',
    'both <- intersect(getNamespaceExports("frmtmb"), getNamespaceExports("brms"))',
    'w <- vapply(both, function(n) find(n)[1L], "")',
    'cat("shared names:", length(both), " now resolving to frmtmb:",',
    '    sum(w == "package:frmtmb"), "\n")',
    'set.seed(1234)',
    'd <- data.frame(x = rnorm(80), g = gl(8, 10))',
    'd$y <- 1 + 2 * d$x + rnorm(80)',
    'r <- function(l, e) cat(sprintf("  %-34s %s\n", l,',
    '  tryCatch(paste(class(e), collapse = ","),',
    '           error = function(err) paste("ERR:", conditionMessage(err)))))',
    'r("bf(y ~ x)", bf(y ~ x))',
    'r("lf(sigma ~ x)", lf(sigma ~ x))',
    'r("nlf(s ~ a * x)", nlf(s ~ a * x))',
    'r("mixture(gaussian(), gaussian())", suppressMessages(mixture(gaussian(), gaussian())))',
    'r("get_prior(y ~ x, data = d)", get_prior(y ~ x, data = d))',
    'r("set_prior(\\"normal(0, 1)\\", class = \\"b\\")", set_prior("normal(0, 1)", class = "b"))',
    'r("custom_family(...)", custom_family("f", dpars = "mu",',
    '  links = list(mu = "identity"),',
    '  lpdf = function(y, dpars, aterms) RTMB::dnorm(y, dpars$mu, 1, log = TRUE)))',
    'r("prior(normal(0, 1), class = b)", prior(normal(0, 1), class = "b"))',
    'r("student()", student())',
    'f <- frm(bf(y ~ x + (1 | g)), data = d, family = gaussian())',
    'r("frm(bf(...))", f)',
    'grDevices::pdf(NULL)',
    'for (n in c("conditional_effects", "pp_check", "ngrps", "fixef",',
    '            "ranef", "VarCorr")) r(paste0(n, "(fit)"), get(n)(f))',
    'r("hypothesis(fit, \\"x\\")", hypothesis(f, "x"))',
    'grDevices::dev.off()'
  ), child)
  out <- system2(file.path(R.home("bin"), "Rscript"), shQuote(child),
                 stdout = TRUE, stderr = TRUE)
  cat(paste(out, collapse = "\n"), "\n")
  "ran"
}, "BEHAVIOR",
"with frmtmb attached second every shared name resolves to frmtmb, so bf(), lf(), nlf(), mixture(), get_prior(), set_prior() and custom_family() all give the frmtmb object and the copied brms line silently changes meaning in the other direction; prior() still resolves to brms because frmtmb does not export it")

## ============ PATH 2: sampling (frm_sample) ============
#
# The FIT route, because the model needs no priors. The point here is
# not the posterior: it is that every generic below reaches the
# frmtmb_draws method rather than brms's brmsfit method, which is the
# case most likely to go wrong because a draws object is the kind of
# thing brms's methods expect to own.

s_g <- bv("model", "SAMPLE: frm_sample(fit_g) with brms attached", {
  frm_sample(fit_g, chains = 1, iter = 400, warmup = 200,
             seed = 1, cores = 1, refresh = 0)
}, "BEHAVIOR",
"frm_sample() prints a note that the random-effect block stays centered because its sd has a flat prior, and names set_prior(class = 'sd'); brms would have supplied that prior itself")

s_b <- bv("model", "SAMPLE: frm_sample(fit_b) with brms attached", {
  frm_sample(fit_b, chains = 1, iter = 400, warmup = 200,
             seed = 1, cores = 1, refresh = 0)
}, NA_character_, "")

draws_row <- function(nm, obj, what, edge = NA_character_, why = "") {
  bv("post", paste0("SAMPLE: ", nm, "(", what, " draws) with brms attached"), {
    grDevices::pdf(NULL); on.exit(grDevices::dev.off())
    fn <- get(nm)
    val <- if (identical(nm, "hypothesis")) fn(obj, "x > 0") else fn(obj)
    stopifnot(!inherits(val, "brmsfit"))
    cat(nm, "->", paste(class(val), collapse = ", "), "\n")
    class(val)[1]
  }, edge, why)
}

for (nm in c("summary", "fixef", "ranef", "VarCorr", "posterior_epred",
             "posterior_predict", "posterior_linpred", "log_lik", "loo",
             "waic", "bayes_R2", "pp_check", "mcmc_plot",
             "conditional_effects", "ngrps", "nchains", "ndraws",
             "posterior_summary", "posterior_interval", "as.mcmc",
             "hypothesis")) {
  draws_row(nm, s_g, "gaussian")
}

# The binomial fit gets the shorter list: what is being checked is that
# dispatch does not depend on the family.
for (nm in c("summary", "fixef", "pp_check", "loo", "conditional_effects",
             "posterior_predict", "hypothesis")) {
  draws_row(nm, s_b, "binomial")
}

# posterior's own entry points, which brms also registers methods on.
bv("post", "SAMPLE: posterior::as_draws_df() with brms attached", {
  x <- posterior::as_draws_df(s_g)
  cat("dim:", paste(dim(x), collapse = " x "), "\n")
  class(x)[1]
}, NA_character_, "")

# The deprecated brms spellings are the one place where the collision
# changes the ANSWER rather than the object: brms's generic reaches
# frmtmb's method, and frmtmb's method refuses on purpose.
for (nm in c("LOO", "WAIC", "parnames", "posterior_samples", "nsamples")) {
  local({
    n <- nm
    bv("post", paste0("SAMPLE: ", n, "(draws) with brms attached"), {
      suppressWarnings(get(n)(s_g))
    }, "REFUSAL",
    "the deprecated brms spelling is refused with the modern frmtmb name; the message points right, and it is reached through brms's own generic")
  })
}

bv_done()
