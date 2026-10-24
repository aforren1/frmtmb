# Release fixes from the release review (dev/reviews/2026-10-23-release.md,
# verdict FIX), sourced by dev/rel074-release-fix.R, which
# dev/rel074-resolve.R sources, so a merge rerun applies them. Each item
# is by the review's finding number; the coordinator asked for each.
# The NEWS items are in dev/rel074-news-split.R and the leads
# (dev/rel074-news-*-lead.md), since the NEWS merge runs after this
# file.

# Finding 1 (major): the README did not know frmtmb.spatial and listed
# car() and spde() as core's. The user decided on 2026-10-23 that there
# is no meta-package, so the README and the brms-migration vignette
# document one way to install every extension
edit_file("README.md", r"---[Seven companion packages live in this repository, each installed the
same way and each documented on its own part of the site.]---",
r"---[Eight companion packages live in this repository, each installed the
same way and each documented on its own part of the site.]---")
edit_file("README.md", r"---[| [frmtmb.survival](https://aforren1.github.io/frmtmb/frmtmb.survival/) | survival families: brms's `cox()` with an M-spline baseline, and the flexible parametric `royston_parmar()` |]---",
r"---[| [frmtmb.survival](https://aforren1.github.io/frmtmb/frmtmb.survival/) | survival families: brms's `cox()` with an M-spline baseline, and the flexible parametric `royston_parmar()` |
| [frmtmb.spatial](https://aforren1.github.io/frmtmb/frmtmb.spatial/) | spatial random fields: brms's `car()` and the SPDE Matern field `spde()` |]---")
install_all <- r"---[# every extension, after frmtmb itself; frmtmb.eam goes before
# frmtmb.learn, which imports it
for (p in c("eam", "latent", "ode", "survival", "spatial", "coupling",
            "sample", "learn")) {
  remotes::install_github("aforren1/frmtmb",
                          subdir = paste0("extensions/frmtmb.", p))
}
```]---"
edit_file("README.md", r"---[remotes::install_github("aforren1/frmtmb", subdir = "extensions/frmtmb.survival")
remotes::install_github("aforren1/frmtmb", subdir = "extensions/frmtmb.learn")
remotes::install_github("aforren1/frmtmb", subdir = "extensions/frmtmb.coupling")
```]---",
paste0(r"---[remotes::install_github("aforren1/frmtmb", subdir = "extensions/frmtmb.survival")
remotes::install_github("aforren1/frmtmb", subdir = "extensions/frmtmb.spatial")
remotes::install_github("aforren1/frmtmb", subdir = "extensions/frmtmb.learn")
remotes::install_github("aforren1/frmtmb", subdir = "extensions/frmtmb.coupling")
```

frmtmb has no package that installs everything. To install every
extension at once:

```r
]---", install_all))
edit_file("README.md", r"---[- `car(M, gr = g, type =)` fits spatial Gaussian Markov random
  fields with brms's spelling and all four of its types (`escar`,
  `esicar`, `icar`, and `bym2` with brms's scaling convention).
  `spde(fem, gr = node)` fits a Matern field over a finite-element
  mesh, with `fmesher`/INLA matrices as fixed data. The package
  assembles precisions on the tape from fixed sparse matrices.
  Every normalizing constant is analytic.]---",
r"---[- brms's `car(M, gr = g, type =)` spatial Gaussian Markov random
  field, with all four of its types (`escar`, `esicar`, `icar`, and
  `bym2` with brms's scaling convention), and `spde(fem, gr = node)`,
  a Matern field over a finite-element mesh, are in frmtmb.spatial.
  They register through core's covariance-structure seam,
  `frmtmb_register_covstruct()`.]---")
edit_file("vignettes/brms-migration.Rmd", r"---[frequentist (Wald and profile intervals, likelihood-ratio tests, AIC).

## What ports directly]---",
paste0(r"---[frequentist (Wald and profile intervals, likelihood-ratio tests, AIC).

## Installing frmtmb and its extensions

brms holds every family and term in one package. frmtmb holds the
general ones in core and the ones for one domain in extension
packages: for example `cox()` is in frmtmb.survival, `car()` in
frmtmb.spatial and `wiener()` in frmtmb.eam. Without the package, such
a spelling stops with a message that names it. frmtmb has no package
that installs everything. To install frmtmb and every extension:

```r
remotes::install_github("aforren1/frmtmb")
]---", install_all, r"---[

Then load the extensions a model needs, for example
`library(frmtmb.survival)`.

## What ports directly]---"))

# Finding 2: the README's first example stopped with an error (since
# before this round); the review's working spelling, run at the
# consolidation (dev/rel074-readme-example.R)
edit_file("README.md", r"---[hypothesis(fit, "sd_Subject__Days^2 / (sd_Subject__Days^2 + sigma^2)",
           method = "boot")]---",
r"---[hypothesis(fit, "sd_Subject__Days^2 / (sd_Subject__Days^2 + sigma^2) > 0",
           class = NULL, method = "boot")]---")

# Finding 3: a cox() fit an older frmtmb built records frmtmb as its
# package, so require_model_packages() let it through and update()
# failed on a function core no longer has. Refused by name now, with
# what to do (dev/rel074-oldcox.R reads a real 0.73.0 fit)
edit_file("R/extension-home.R", r"---[  for (fam in fams) {
    if (!inherits(fam, "frmtmb_family")) next
    for (p in family_packages(fam)) {]---",
r"---[  for (fam in fams) {
    if (!inherits(fam, "frmtmb_family")) next
    home <- old_core_family_home(fam)
    if (!is.null(home)) {
      frm_stop("The ", fam[["family"]], " family of this model was built ",
               "by an older frmtmb, when the family was in frmtmb. It is ",
               "in the extension package ", home, " now, and the ",
               "functions this fit calls are no longer in frmtmb. Refit ",
               "the model with library(", home, ") loaded",
               class = "frmtmb_extension_missing", call. = FALSE)
    }
    for (p in family_packages(fam)) {]---")
edit_file("R/extension-home.R", r"---[#' Check every family of a fit, a frame, a spec or a single family with
#' `require_extension()`.]---",
r"---[#' The extension that now holds a family an older frmtmb built, or
#' `NULL`. Such a family has a name in `extension_family_homes`, records
#' frmtmb as its package, and has its density in frmtmb's namespace: a
#' `cox()` fit saved by frmtmb 0.73.0 or earlier. A family the extension
#' built records the extension, and a user's family of the same name has
#' its density elsewhere, so neither is matched; a `wiener()` fit an
#' older frmtmb.eam built records frmtmb.eam and still works.
#'
#' @noRd
old_core_family_home <- function(fam) {
  nm <- fam[["family"]]
  if (!is.character(nm) || length(nm) != 1L ||
        !nm %in% names(extension_family_homes) ||
        !identical(frm_family_package(fam), "frmtmb")) {
    return(NULL)
  }
  f <- fam[["lpdf"]]
  if (!is.function(f) ||
        !identical(topenv(environment(f)), asNamespace("frmtmb"))) {
    return(NULL)
  }
  unname(extension_family_homes[[nm]])
}

#' Check every family of a fit, a frame, a spec or a single family with
#' `require_extension()`.]---")
edit_file("tests/testthat/test-extension-home.R",
r"---[test_that("a family from an absent package is refused at the fit", {]---",
r"---[test_that("a family an older frmtmb built, now an extension's, is refused", {
  # A core family under a moved family's name stands for a cox() fit
  # saved by frmtmb 0.73.0: frmtmb is its package and its density is in
  # frmtmb's namespace. frmtmb's dev/rel074-oldcox.R reads a real one
  fam <- owned_by(as_frmtmb_family(gaussian()), "frmtmb")
  old <- fam
  old[["family"]] <- "cox"
  e <- tryCatch(require_model_packages(old), error = identity)
  expect_s3_class(e, "frmtmb_extension_missing")
  expect_match(conditionMessage(e),
               "The cox family of this model was built by an older frmtmb",
               fixed = TRUE)
  expect_match(conditionMessage(e),
               "Refit the model with library(frmtmb.survival) loaded",
               fixed = TRUE)
  # the guard absent: the core family itself; the same name recorded as
  # the extension's (a frmtmb.eam wiener() fit records frmtmb.eam); a
  # user's family of that name, whose density is not frmtmb's
  expect_null(old_core_family_home(fam))
  expect_true(require_model_packages(fam))
  ext <- old
  ext[["family"]] <- "wiener"
  expect_identical(old_core_family_home(ext), "frmtmb.eam")
  expect_null(old_core_family_home(owned_by(ext, "frmtmb.eam")))
  usr <- old
  lp <- usr[["lpdf"]]
  environment(lp) <- globalenv()
  usr[["lpdf"]] <- lp
  expect_null(old_core_family_home(usr))
})

test_that("a family from an absent package is refused at the fit", {]---")

# Finding 5: the one-session runner did not watch the covariance
# structure registry this round added (the review's probe test leaked
# "covstruct registry: 22 -> 23 entries" with the line, nothing without)
edit_file("dev/release/run-session.R",
r"---[    `aterm registry` = reg("frmtmb_aterm_registry")$reg,]---",
r"---[    `aterm registry` = reg("frmtmb_aterm_registry")$reg,
    `covstruct registry` = reg("frmtmb_covstruct_registry")$reg,]---")

# Finding 6: a CI job keeps the two helper copies identical
wr(rd("dev/rel074-absent-helper-identical.yaml"),
   ".github/workflows/absent-helper-identical.yaml")

# Finding 9: messages
edit_file("R/parse.R", r"---[                 deparse1(tm), ". The supported residual correlation ",
                 "terms are ar(), ma(), arma(), cosy(), unstr() and car()",]---",
r"---[                 deparse1(tm), ". The supported residual correlation ",
                 "terms are ar(), ma(), arma(), cosy() and unstr(); the ",
                 "car() field is in the extension package frmtmb.spatial",]---")
edit_file("R/priors.R", r"---[             "class = \"simo\" (a mo() simplex) or class = \"sbhaz\" (a ",
             "cox() baseline), not class = \"", class, "\"", call. = FALSE)]---",
r"---[             "class = \"simo\" (a mo() simplex) or class = \"sbhaz\" (the ",
             "baseline of a cox() model, frmtmb.survival), not class = \"",
             class, "\"", call. = FALSE)]---")
edit_file("R/priors.R", r"---[      frm_stop("class = \"sbhaz\" is the whole baseline simplex of a cox() ",
               "response, so coef and group name nothing it can narrow to",
               call. = FALSE)]---",
r"---[      frm_stop("class = \"sbhaz\" is the whole baseline simplex of a cox() ",
               "response (frmtmb.survival), so coef and group name nothing ",
               "it can narrow to", call. = FALSE)]---")
edit_file("R/priors.R", r"---[             "class = \"simo\" (a mo() simplex), class = \"sbhaz\" (a ",
             "cox() baseline) or a simplex a family declares, not class = ",]---",
r"---[             "class = \"simo\" (a mo() simplex), class = \"sbhaz\" (the ",
             "baseline of a cox() model, frmtmb.survival) or a simplex a ",
             "family declares, not class = ",]---")

# Finding 10: site and CI comments, and a Suggests version
edit_file("_pkgdown.yml", r"---[  - ["reference/cox_baseline.html", "frmtmb.survival/reference/cox_baseline.html"]]---",
r"---[  - ["reference/cox_baseline.html", "frmtmb.survival/reference/cox_baseline.html"]
  - ["reference/mixture_mvn.html", "frmtmb.latent/reference/mixture_mvn.html"]]---")
# no lane changed this file, so it starts from its HEAD blob
restore(".github/workflows/check-frmtmb-learn.yaml")
edit_file(".github/workflows/check-frmtmb-learn.yaml", r"---[# ONE SUGGEST IS DELIBERATELY NOT INSTALLED, unlike the spline job,
# which installs every one of its.]---",
r"---[# ONE SUGGEST IS DELIBERATELY NOT INSTALLED, unlike the survival and
# spatial jobs, which install every one of theirs.]---")
edit_file(".github/workflows/check-frmtmb-learn.yaml",
r"---[      # purpose, as in the spline job. Nothing in this suite is gated on]---",
r"---[      # purpose, as in the survival job. Nothing in this suite is gated on]---")
edit_file(".github/workflows/check-frmtmb-spatial.yaml", r"---[# No extra-repositories here: every dependency of this package (Matrix,
# RTMB, brms, rstan) is on CRAN, and the package compiles nothing. brms]---",
r"---[# No extra-repositories here: every dependency of this package (Matrix,
# RTMB, brms, fmesher, rstan) is on CRAN, and the package compiles
# nothing. fmesher serves the spde() examples and tests. brms]---")
edit_file("extensions/frmtmb.sample/DESCRIPTION", "    frmtmb.latent,\n",
          "    frmtmb.latent (>= 0.7.0),\n")
cat("RELEASE REVIEW FIXES APPLIED\n")
