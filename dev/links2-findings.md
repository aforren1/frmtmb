# links2 lane: making the link registry reachable and findable

Worktree `frmtmb-wt-links2`, branch `wt-links2`, base `c18253e`.
Core only. brms 2.23.0 read from the user library for every comparison.

## Survey: the four defects, verified

### Defect 1, no documented home. CONFIRMED

    grep -rln "probit|softplus|cauchit|squareplus|softit" man/   ->  0 of 82 files

`man/frmtmb-families.Rd` documents the argument as `@param link Link
for mu.` (R/families.R:4646) and stops. Nothing enumerates the roster.
The only user-facing mentions are one NEWS bullet and the
brms-migration vignette. `docs/reference/` inherits the same hole
because it is built from `man/`.

### Defect 2, no link on any dpar but the mean. CONFIRMED

Every constructor is `function(link = "...")` and every non-mu link is
a string literal in the family body:

    R/families.R:821   fam_gaussian     links = list(mu = link, sigma = "log")
    R/families.R:1096  fam_negbinomial  links = list(mu = lk, shape = "log")
    R/families.R:1058  fam_student      links = list(mu = link, sigma = "log", nu = "logm1")

`grep -rn "link_sigma|link_shape|link_nu" R/ vignettes/` returns
nothing: the brms spelling is absent from the package entirely.

The good news, measured before touching anything: the machinery below
the grammar is ALREADY link-generic, so this is a signature and
validation change, not a plumbing one.

    R/objective.R:326,344   dparv[[..]] <- lp[["link"]]$linkinv(eta)
    R/fit.R:1950            val <- lp[["link"]]$linkfun(raw)      # init
    R/frame.R:2485          lp[["link"]]$linkfun(lp[["constant"]])
    R/predict.R (12 sites)  lp[["link"]]$linkinv(eta)
    R/priors.R:1911         coef_placement() special-cases log and
                            identity and falls through to the link
                            object's own linkinv/mu_eta for every other

That last one is the priors interaction the task warns about, and it
is already written generically. It still has to be TESTED, because
"generic by inspection" is not "generic in fact".

### Defect 3, ordinal links are a hand-written switch. CONFIRMED

Current lines (the task's 2352/2600 have moved):

    R/families.R:2401  fam_cumulative()  inline switch(link, logit=, probit=)
    R/families.R:2576  ord_num_cdf()     switch, for the simulators
    R/families.R:2649  ord_link_cdf()    switch, for the taped density

and `robust <- identical(link, "logit")` at :2408, :2681, :2729.

### Defect 4, summary() never names the link. CONFIRMED

    R/methods-fit.R:169  cat("Family:", x$family[["family"]], "\n")
    R/methods-fit.R:52   print.frmtmb_fit, family name only

`print.frmtmb_family` (R/families.R:4877) DOES print links, so the
information is one line away from the two places a reader looks.

## brms 2.23.0, per-dpar link roster (measured, not recalled)

`brms:::links_dpars(dpar)`. The first entry is brms's default.

    dpar       links accepted
    ---------- ---------------------------------------
    sigma      log, identity, softplus, squareplus
    shape      log, identity, softplus, squareplus
    phi        log, identity, softplus, squareplus
    kappa      log, identity, softplus, squareplus
    beta       log, identity, softplus, squareplus
    disc       log, identity, softplus, squareplus
    bs         log, identity, softplus, squareplus
    ndt        log, identity, softplus, squareplus
    nu         logm1, identity
    zi         logit, identity
    hu         logit, identity
    zoi        logit, identity
    coi        logit, identity
    bias       logit, identity
    quantile   logit, identity
    xi         log1p, identity
    alpha      identity, log, softplus, squareplus
    theta      identity
    mu         (from the family table, not this one)

`nu` is the trap: `links_dpars("nu")` is student's degrees of freedom
and answers `logm1`, but frmtmb's `compois` also calls its dispersion
`nu` and it is an ordinary positive parameter. Keying a validator on
the dpar NAME alone would put `logm1` on a compois dispersion. The
allowed set therefore has to come from the call site, which knows the
parameter's support, not from a name lookup.

## Defect 3 done: the ordinal families go through the registry

`ord_link()` (R/families.R:2585) resolves an ordinal `link` against the
distribution functions its thresholds can be read through, and
`ord_cdf()` (R/families.R:2612) is the one CDF both the taped density
and the plain simulators use. The three hand-written switches are
gone. Measured acceptance now, against brms 2.23.0:

    link            cumulative  sratio   cratio   acat     brms agrees
    -------------   ----------  -------  -------  -------  -----------
    logit           ok          ok       ok       ok       yes
    probit          ok          ok       ok       refused  yes
    probit_approx   ok          ok       ok       refused  yes
    cloglog         ok          ok       ok       refused  yes
    cauchit         ok          ok       ok       refused  yes
    softit          ok          refused  refused  refused  yes
    log             refused     refused  refused  refused  yes

softit for cumulative and not for the sequential pair is brms's own
asymmetry, kept rather than tidied.

### acat is refused off the logit, and that is not laziness

`brms:::inv_link_acat()` BRANCHES on the link. On the logit the
category probability is the log-linear form frmtmb implements; off it,
brms switches to a product of CDFs times a reversed product of
survivals. That is a different density, not a different link, so
routing the registry through `fam_acat()` would silently fit a model
the user did not ask for. Refused by name, with the reason in the code
at R/families.R:2860.

### The robust branch generalized, and it is what makes the new links usable

The logit-only log-space form is now taken at the link's own log odds.
`logit_eta` is by definition the log odds of that link's CDF, so
`F(x) = plogis(logit_eta(x))` exactly, and the logistic algebra
`log(F(a) - F(b)) = logspace_sub(-qb, -qa) + log_inv_logit(qa) +
log_inv_logit(qb)` then holds for every link that carries the field.
For the logit `q` is the identity, so the code reduces to what was
there and the logit path is unchanged.

That assumption was checked against a 500-bit Rmpfr reference:

    is plogis(logit_eta(x)) == linkinv(x), x over [-30, 30]?
      logit 0        probit 2.8e-14   probit_approx 3.6e-14
      cloglog 2.6e-15   softit 1.8e-15

and the resulting log density, interior category, against the same
reference. "plain" is what the old code would do off the logit:

    link            eta range     plain        robust
    -------------   -----------   ----------   ----------
    logit           [-2, 2]       8.9e-16      4.4e-16
    logit           [-8, 8]       2.5e-13      8.9e-16
    logit           [-25, 25]     1.2e-06      3.6e-15
    probit          [-8, 8]       1.6e-04      3.6e-15
    probit_approx   [-8, 8]       Inf          7.1e-15
    probit_approx   [-25, 25]     Inf          1.1e-13
    cloglog         [-8, 8]       Inf          1.8e-15
    cloglog         [-25, 25]     Inf          3.6e-15
    softit          [-25, 25]     4.0e-14      8.0e-15

`Inf` is a log density of `-Inf` with a dead gradient. Without the
generalization cloglog and probit_approx would be reachable and
useless past a linear predictor of 2, so this is not a refinement, it
is the feature. cauchit has no `logit_eta`, needs none, and keeps the
plain path: its tails are polynomial.

`ord_tau_init()` (R/families.R:2550) now maps the observed cumulative
proportions through the family's OWN link rather than always through
`qlogis`, so a probit fit no longer starts on logit-scale thresholds.

One latent bug fixed on the way: `fam_cratio()`'s robust branch read
its pair at `M` and relied on `1 - F(-x) = F(x)`. That is true of the
logistic, the normal, the Cauchy and probit_approx's odd cubic, and
FALSE of the cloglog. It now reads at `-M` and holds for an asymmetric
CDF (R/families.R:2770).

## Defect 2 done: a link on every distributional parameter

26 internal constructors and 22 exported ones gained `link_<dpar>`
arguments, one per non-mean parameter, validated against the set that
parameter's support admits (`dpar_link()`, R/links.R:471). Measured:
every exported constructor accepts every link brms allows for each of
its dpars, and refuses one outside the set, 0 failures over the whole
cross product.

    constructor                  arguments added
    ---------------------------  ------------------------------------
    student                      link_sigma, link_nu
    skew_normal                  link_sigma, link_alpha
    exgaussian                   link_sigma, link_beta
    shifted_lognormal            link_sigma, link_ndt
    asym_laplace                 link_sigma, link_quantile
    zero_inflated_asym_laplace   link_sigma, link_quantile, link_zi
    hurdle_lognormal             link_sigma, link_hu
    lognormal, huber             link_sigma
    negbinomial, weibull         link_shape
    hurdle_gamma                 link_shape, link_hu
    zero_inflated_negbinomial    link_shape, link_zi
    nbinom1, Beta, tweedie,      link_phi
      beta_binomial
    zero_inflated_beta           link_phi, link_zi
    compois                      link_nu
    von_mises                    link_kappa
    zero_inflated_poisson,       link_zi
      zero_inflated_binomial
    hurdle_poisson               link_hu

### Two helpers had to be fixed first, or this would have been silent

`log_dpar()` returned the linear predictor AS the log of the dpar. That
is right only on a log link, and its own comment said so ("Do not call
it on a dpar whose link the user chooses"). The moment `link_shape`
exists, a softplus predictor would have been read as a log and
`negbinomial()` would have returned a wrong density with no error. It
now goes through the link's `log_eta` and falls back to `log(value)`
(R/families.R:447). `gate_logs()` had the identical problem for `zi`,
`hu` and `quantile`, assuming logit (R/families.R:458). 13 call sites
updated.

### What was already generic, measured rather than assumed

Nothing below the grammar needed changing. The objective, predict,
init and the frame all apply `link$linkinv` / `link$linkfun` from the
resolved link object already.

### The priors pairing, tested

`set_prior(class = "sigma")` is a density on sigma itself, so the
objective, which works on the link scale, carries the Jacobian
`|d sigma / d eta|`. That Jacobian has to be the dpar's OWN link. It
is. The MAP was compared against an independent posterior computed
outside frmtmb as `loglik + log p(sigma) + log|mu_eta(eta)|` with beta
profiled out at OLS:

    link_sigma   frmtmb        independent   abs diff
    log          1.44384397    1.44384448    5.1e-07
    identity     1.44050864    1.44050827    3.8e-07
    softplus     1.44199680    1.44199537    1.4e-06
    squareplus   1.44267077    1.44267122    4.5e-07

All four differ, which is what a correctly applied Jacobian does to a
MAP; all four match their own link's independent calculation to
optimizer tolerance. Had the placement assumed a log link, the
softplus row would have matched the log row instead.

### The stats-owned families needed one new export

`gaussian()`, `poisson()`, `binomial()` and `Gamma()` come from
'stats'. They have no `link_sigma` and `stats::gaussian()` also
refuses softplus, squareplus, probit_approx and softit for the mean.
frmtmb does not export replacements, and must not: a `gaussian()` that
answered a `frmtmb_family` would break every `glm()` call in an
attached session. brms has the same hole and fills it with
`brmsfamily()`. `frm_family()` (R/families.R:4595) is the equivalent:

    frm_family("gaussian", link_sigma = "softplus")
    frm_family("Gamma", link = "inverse", link_shape = "identity")
    frm_family("poisson", link = "softplus")

## Defect 1 done: the links have a documented home

`?\`frmtmb-links\`` (roxygen at the head of R/links.R, rendered to
man/frmtmb-links.Rd) carries four sections:

* the roster, one row per link, with `g`, `g^-1`, the range the
  inverse maps onto, and whether the link carries a robust field;
* which links each family group takes for its mean, which is brms
  2.23.0's table as measured from `.family_*()`;
* which set each distributional parameter admits, from
  `brms:::links_dpars()`;
* what a robust field buys, and why cauchit, inverse and `1/mu^2` do
  not have one.

Cross-references added: the shared `@param link` of every family
constructor, `frmtmb_family(links =)`, `frm_family()`, and the
pkgdown reference index under Families. `get_link()`'s unknown-link
error now ends "See ?`frmtmb-links` for what each one maps and which
families take it" (R/links.R:441); it already listed the names, and
what it was missing was where to find out what they mean.

## Defect 4 done: summary() names the link

`family_link_str()` (R/families.R:4940) renders one family's links as
`mu = identity; sigma = log`. `print.summary.frmtmb_fit()` and
`print.frmtmb_fit()` show it on a ` Links:` line under `Family:`, and
a multivariate fit prefixes each response.

An ordinal fit reports `cdf = cloglog` instead. Its `mu` genuinely
carries an identity link, because the distribution function applies to
`tau_j - eta` and not to `eta`, so printing `mu = identity` would
answer a question the reader did not ask.

## Item 5 done: the bcm probit workaround is retired

`bcm_probit()`, `bcm_binomial_probit()` and `bcm_gaussian_probit()`
are deleted from inst/bcm/binomial-extras.R. Both families collapsed
into core constructors; neither had to stay a custom family.

    bcm_binomial_probit()  ->  binomial(link = "probit")
    bcm_gaussian_probit()  ->  gaussian(link = "probit") + se(sd)
    bcm_mpt_pairs(link = bcm_probit())  ->  bcm_mpt_pairs(link = "probit")

`stats::binomial()` and `stats::gaussian()` both accept the name
"probit" through `make.link()`, and frmtmb resolves it against its own
registry, so no new export was needed for these.

Equivalence measured before switching the call sites, on the same
data:

    binomial probit   d logLik 1.4e-14   d coef 7.1e-12   d se 8.4e-14
    gaussian probit   d logLik 0         d coef 0         d se 0

The gaussian collapse needed `vreal(sd)` to become `se(sd)`. That does
NOT touch seam 2: R/frame.R is unchanged and `se()` is still gated on
the family name. The model simply stopped being a custom family, so
the gate stopped applying to it. Seams 2, 3 and 4 are untouched, and
dev/bcm-findings.md now says so at the head of its Reproductions
section.

### Stan identity residuals, before and after

Three files gated, `FRMTMB_BRMS_FIT_TESTS=true`, against my own copy
of the Stan cache. Zero failures both times. Eleven identities; the
stated constant is 0 in all of them.

    #   file                       before            after
    1   signal-detection SDT_1     -2.344791e-13     -2.593481e-13
    2   signal-detection SDT_2      0                 0
    3   signal-detection SDT_3      0                 0
    4   binomial Rate_3            -6.217249e-15     -6.217249e-15
    5   binomial Rate_4            -1.243450e-14     -1.243450e-14
    6   binomial Rate_5            -1.332268e-14     -1.332268e-14
    7   binomial (4th)             -4.440892e-16     -4.440892e-16
    8   binomial (5th)              0                 0
    9   binomial Survey            -7.283063e-14     -7.283063e-14
    10  esp OptionalStopping        3.552714e-15      3.552714e-15
    11  esp Extraversion           -2.853717e-12     -2.841283e-12

Nine of eleven are bit-identical. The two that moved are the two the
change touches: SDT_1 is the probit binomial, now scored through
`dbinom_robust()` off the registry's `logit_eta` rather than through a
plain round trip, and Extraversion is the gaussian that moved from
`vreal()` to `se()`. Both moved within 1e-13 and both are still far
inside the harness tolerance.

## What I refused, and why

1. **`acat()` off the logit.** brms's `inv_link_acat()` branches: the
   logit gets the log-linear form frmtmb implements, and every other
   link gets a product of distribution functions times a reversed
   product of survivals. Substituting a CDF into the existing
   expression would fit a different model and report it as the same
   one. Implementing brms's second form is a new density with its own
   robustness question, not a routing change.

2. **Gating a family's MEAN link.** brms refuses `poisson(link =
   "logit")`; frmtmb accepts any registry link for `mu` and still
   does. An extension family is free to mean something else by its own
   `mu`, and the gate would be a behavior change nobody asked for that
   could refuse a working extension model. The per-family table is
   documented as what PORTS rather than as what is allowed. SEAM: if a
   later lane wants the gate, the table is already measured and in
   `?\`frmtmb-links\``.

3. **The ordinal `disc` parameter.** brms's cumulative, sratio, cratio
   and acat carry a `disc` (discrimination) dpar with a log link;
   frmtmb's carry `mu` alone. That is a missing PARAMETER, not a
   missing link, so `link_disc` would be an argument for something
   that does not exist. SEAM, and the bigger of the two ordinal gaps.

4. **Exporting `gaussian()`, `poisson()`, `binomial()` or `Gamma()`.**
   It is the obvious way to give them `link_sigma`, and it would break
   every `glm()` call in a session with frmtmb attached, because the
   constructor would answer a `frmtmb_family` rather than a
   `stats::family`. brms does not do it either. `frm_family()` instead.

5. **A `link_power` on `tweedie()`.** The tweedie power lives on the
   open interval from one to two, `power12` is the only link in the
   roster that maps onto it, and brms has no tweedie to port from, so
   the argument would have exactly one legal value.

6. **`log1p`, which no constructor can reach.** It is in the registry
   for brms's `xi` of `gen_extreme_value`, and frmtmb has no
   `gen_extreme_value` family. The link is documented and correct;
   there is nothing to attach it to. SEAM, for whoever adds that
   family.

7. **Anything owned by a sibling lane.** `R/frame.R`'s `se()` and
   `cens()` gates and the addition-term allow-list are untouched, and
   so is `R/predict.R`'s band code. Retiring `bcm_gaussian_probit()`
   moved one model from `vreal()` to `se()`, which USES the gate that
   already admits gaussian; it did not change the gate.

8. **Rebuilding `docs/`.** The pkgdown site is regenerated as its own
   release step (see the last two commits on main), and this lane
   changed no rendered page beyond adding two. `_pkgdown.yml` lists
   `frmtmb-links` and `frm_family` under Families, so the next site
   build picks them up. Until then `docs/reference/` has no
   frmtmb-links.html, which is worth knowing when reading the defect 1
   verification.

## The decisive ordinal cross-check

The taped log-density and the numeric category distribution are
written independently: the lpdf takes the generalized log-space form
off each link's log odds, `ord_cat_probs()` builds the whole
distribution from the plain CDF. Agreement of `logLik()` with
`sum(log(P[cbind(i, y_i)]))` therefore checks both. 600 rows, four
categories:

    family       link            logLik            |ll - sum log P[y]|
    -----------  --------------  ----------------  -------------------
    cumulative   logit           -713.2846451249   1.1e-13
    cumulative   probit          -715.0652318056   0
    cumulative   probit_approx   -715.0887012252   4.5e-13
    cumulative   cloglog         -725.5188785548   2.3e-13
    cumulative   cauchit         -719.2326343665   4.5e-13
    cumulative   softit          -731.9617405992   1.1e-13
    sratio       logit           -719.7289034755   2.3e-13
    sratio       probit          -719.7202303673   4.5e-13
    sratio       probit_approx   -719.7150526748   4.5e-13
    sratio       cloglog         -725.5188785548   4.5e-13
    sratio       cauchit         -720.7539944544   0
    cratio       logit           -719.7289034767   1.1e-13
    cratio       probit          -719.7202303684   1.0e-12
    cratio       probit_approx   -719.7150526748   0
    cratio       cloglog         -717.6209811328   2.3e-13
    cratio       cauchit         -720.7539944556   8.0e-13

One row of that table is worth reading twice. `cumulative cloglog` and
`sratio cloglog` land on the SAME log-likelihood to ten decimal
places, -725.5188785548, from two code paths with nothing in common:
an interior difference of two distribution functions on one side, a
product of hazards on the other. That is the textbook identity that
the proportional-hazards cumulative model and the stopping-ratio model
coincide on the complementary log-log link. Neither implementation was
written with the other in mind, and the logit pair differs by 6.4 nats
on the same data, so it is not an artifact of the fit. Both are in
`tests/testthat/test-ordinal.R`.

## The brms identity for the new ordinal links

Defect 3 asked for an identity against brms wherever brms accepts the
same model. brms accepts the whole roster, so five rows were added to
`tests/testthat/test-brms-likelihood.R:255-268`, inside the existing
"row 12: ordinal families" block that already pinned the threshold
parameterization:

    cumulative("probit")     cumulative("cloglog")   cumulative("cauchit")
    sratio("cloglog")        cratio("cloglog")

`brms_lp_check()` compares frmtmb's joint log density with
`rstan::log_prob()` of the Stan program brms GENERATES for the same
formula, family and data, at frmtmb's own estimate mapped into Stan's
parameter blocks, and asserts the stated constant and a vanishing
gradient. Five new Stan programs compiled into my own cache copy
(91 to 96 entries).

    gated test-brms-likelihood.R:  404 pass, 0 fail, 0 error, 0 skip

`cratio("cloglog")` is the row that earns its keep. It is exactly the
case the old robust branch got wrong, because that branch read its
distribution function at `tau - eta` and relied on `1 - F(-x) = F(x)`,
which the cloglog does not satisfy. brms has no such shortcut, so this
row is what would have caught the bug had it shipped.

The threshold parameterization holds unchanged for every new link: the
same block asserts it, and `brms_ord_thresholds()` needed no change,
because whether a family stores `(tau_1, log increments)` or the raw
thresholds is a property of the family and not of the link.

## Why both long runs were restarted

The first core-suite and `R CMD check` runs were discarded rather than
reported. `test-message-uniqueness.R` caught a real regression while
they were in flight: `frm_family()` reused `as_frmtmb_family()`'s
"Unsupported family name" string verbatim, and the suite asserts every
condition message template in `R/` is unique. Fixed at
R/families.R:4640, where it now reads "frm_family(): no family called
'...'", which also names the function the caller actually used.

That fix changed `R/` mid-run, so both jobs were then measuring a
source tree that no longer existed. Both were stopped and restarted
from the final tree. Only processes carrying this lane's own `l2-`
prefix were killed.

The restarted check also fixes an environmental NOTE the first run
produced: "Files 'README.md' or 'NEWS.md' cannot be checked without
'pandoc' being installed". pandoc 3.8.3 is present and runs, but
`R CMD check` calls the binary directly and never reads
`RSTUDIO_PANDOC`, so the second run puts it on `PATH` as well.

## R CMD check --as-cran

Built and checked from the final tree, R 4.6.1, pandoc 3.8.3 and
TinyTeX on PATH, manual built.

    * checking tests ... [24m] OK
    Status: 2 NOTEs          (no ERROR, no WARNING)

Both NOTEs named with their cause:

1. `checking HTML version of manual ... NOTE`
   "Skipping checking math rendering: package 'V8' unavailable".
   V8 is not installed on this machine. This is the expected NOTE.

2. `checking examples ... [60s] NOTE`
   Examples over five seconds: `residuals.frmtmb_fit` 6.68s and
   `profile.frmtmb_fit` 5.19s. Pre-existing and not this lane's:
   those topics are documented in R/predict.R and R/confint.R, and
   this lane's diff touches only R/families.R, R/links.R and
   R/methods-fit.R. None of the topics added here
   (`frmtmb-links`, `frm_family`) appears in the list, so the
   examples written this round all run under five seconds. The
   elapsed figures are also load-sensitive: the same NOTE listed four
   topics at 106s while three jobs shared the machine, and two at 60s
   once they did not.

A third NOTE from the first attempt is GONE and was environmental:
"Files 'README.md' or 'NEWS.md' cannot be checked without 'pandoc'
being installed". pandoc 3.8.3 is installed and runs, but R CMD check
invokes the binary directly and never reads `RSTUDIO_PANDOC`. With
pandoc on PATH the step reports `checking top-level files ... OK`.

## Verification, counts audited by name

### Core suite, one file per process

Clean run from the final tree, one `Rscript` per file, 131 processes.

    RESULT lines               131
    distinct files reported    131
    files that ran twice         0
    test files on disk         131
    files with no result line    0
    files with a failure         0

    totals   pass=7705  fail=0  error=0  warn=1  skip=133

The one warning is `test-prior-compat.R:505`, "Optimizer did not
report convergence: singular convergence (7)", inside the test "coef
and group narrow the classes that read them". It is not this lane's.
Measured rather than argued: the base commit c18253e, extracted to
scratch and run through the same harness, gives the identical
`pass=195 fail=0 warn=1` on the same test. That file exercises
R/priors.R, which this lane does not touch.

The 133 skips are the gated tiers (brms, Stan, bcm) that a plain run
does not enable.

### Gated tiers

    test-brms-likelihood.R          404 pass, 0 fail, 0 error, 0 skip
    test-bcm-signal-detection.R,
      test-bcm-binomial.R,
      test-bcm-esp.R                0 fail, 11 identities, residuals
                                    tabulated above

### Documentation

`roxygen2::roxygenise()` run twice on the final tree; the second run
wrote no file, so it is idempotent. `man/frmtmb-links.Rd` renders all
five sections and three tables; `man/frm_family.Rd` renders with its
usage block. Every formal of every function on the
`frmtmb-families` page is documented: 18 of 18, no undocumented
arguments and no documented-but-absent ones.

# Punch round, review 2026-09-08

Verdict PUNCH, three items, all of them defects in this lane's own new
artifacts and all on the findability axis the lane was sent to fix.
That is the right criticism to have received: the substantive work was
correct and the reference material I wrote to explain it was not.

## Item 1, BLOCKER. acat()'s refusal stated a false reason. FIXED

`fam_acat()` routed through the generic `ord_link()` message, so
`acat(link = "probit")` answered:

    acat(): link must be one of "logit". An ordinal link names the
    distribution function the thresholds are read through, so it has
    to map onto (0, 1); you gave character "probit".

Wrong twice. Probit DOES map onto (0, 1), and brms accepts it here.
Measured, brms 2.23.0:

    brms::acat("probit")         ACCEPTED
    brms::acat("probit_approx")  ACCEPTED
    brms::acat("cloglog")        ACCEPTED
    brms::acat("cauchit")        ACCEPTED
    brms::acat("softit")         ACCEPTED

The true reason existed at R/families.R:4861 and in NEWS, but not
where the user lands. `acat()` now has its own resolver,
`acat_link()` (R/families.R:2691), and says:

    acat() takes the 'logit' link only, and not because character
    "probit" is a bad link: brms accepts "probit", "probit_approx",
    "cloglog", "cauchit" and "softit" here, and they all map onto
    (0, 1). frmtmb refuses them because off the logit brms computes a
    category probability from a SECOND expression, a product of
    distribution functions times a reversed product of survivals. It
    agrees with acat's log-linear form when the distribution function
    is logistic, so it is the same model generalized, but it is a
    density frmtmb has not written, and substituting a distribution
    function into the log-linear form does not reach it. cumulative(),
    sratio() and cratio() do take these links. See ?`frmtmb-links`

I took the reviewer's correction to my framing. Calling brms's second
branch "a different density, not a different link" overstated it: that
branch REDUCES to the log-linear form when F is logistic, so it
generalizes the same model. The honest statement, and the one the
message and the docs now make, is that reaching it needs a second
expression written and taped, and this lane declined to write one.

`cumulative()`, `sratio()` and `cratio()` keep the general message,
which is true for them, and a test asserts both halves
(tests/testthat/test-ordinal.R): for all five links acat refuses, the
message must NOT contain "map onto (0, 1); you gave" and MUST contain
"brms accepts" and "has not written", while `cumulative("log")` must
still contain "map onto (0, 1)".

## Item 2. The frmtmb-links topic never recorded the acat gap. FIXED

The topic's table listed six links for `acat` and its prose said the
ordinal families are the enforced ones, so the new reference page
promised exactly what item 1 refuses. The table row is now split
(R/links.R:69-71):

    | cumulative        | logit, probit, probit_approx, cloglog, cauchit, softit |
    | sratio, cratio    | logit, probit, probit_approx, cloglog, cauchit |
    | acat              | logit ONLY. brms takes the same six as cumulative;
                          this is the one place frmtmb departs, and the
                          reason is below |

and a paragraph under the ordinal note gives the reason in full, in the
corrected framing above.

## Item 3. The brms table claimed a pairing that does not port. FIXED

R/links.R:65 put `beta_binomial` in a row ending `identity`, `log`,
under a heading saying the table says which pairings PORT. Measured:

    brms beta_binomial            refuses log
    brms binomial                 ACCEPTS log
    brms beta                     ACCEPTS log
    brms zero_inflated_binomial   ACCEPTS log
    brms zero_inflated_beta       ACCEPTS log
    brms bernoulli                ACCEPTS log

`brms:::.family_beta_binomial()$links` ends at `identity`. One row, one
family. `beta_binomial` is now its own row saying the same list WITHOUT
`log`.

## The two optional findings, both taken

Finding 4, the mixture `Links:` line at 106 characters. TAKEN, because
the line exists to tell a reader which scale a coefficient is on, and a
line the terminal breaks mid-name fails at that. `cat_family_links()`
(R/methods-fit.R:59) now wraps, and breaks BETWEEN entries rather than
on any space: `strwrap()` would leave "nu2" ending one line and
"= logm1" starting the next.

    before  Links: mu1 = identity; ... nu2 = logm1; theta1 = identity   106
    after   Links: mu1 = identity; sigma1 = log; mu2 = identity;
                   sigma2 = squareplus;                                  74
                   nu2 = logm1; theta1 = identity                        38

Measured widths now: distributional 53, ordinal 22, multivariate 50,
two-component mixture 74 then 38. Every print line is within 80.

Finding 5, a multivariate `summary()` printed no `Links:` line at all,
because `print.summary.frmtmb_fit` read the single `x$family` slot,
which is empty for a multivariate fit. TAKEN: `summary()` now stores
`links` rendered off `object$spec$responses` the way `print()` does
(R/methods-fit.R:146), and the line prints:

     Links: y: mu = identity; sigma = log; z: mu = log

The blank `Family:` on that same path is left alone. It is pre-existing
on the base, it is about the family NAME rather than the link, and
changing it would alter multivariate `summary()` output for tests this
lane does not own. SEAM, recorded.

## Count corrected

The lane said 22 exported constructors gained a `link_<dpar>`
argument. Measured, the number is 23; I had missed `huber()`, which
was wired by hand rather than by the script that did the other 22.
23 exported constructors, 34 (constructor, dpar) pairs, 26 internal
`fam_*`. NEWS and this file now say 23.

## Finding 3 of the review, the stale built site. NOT TAKEN, recorded

`docs/articles/bayesian-cognitive-modeling.md` still documents
`bcm_probit()` and the two retired families at lines 331, 377, 382 and
479. The vignette SOURCE is correct; `docs/` is generated.

Not rebuilt, for the same reason as before and one more. `docs/` is
regenerated as a release step, and main has just taken in `wt-xspec`
and `wt-loose`, so a rebuild from this worktree would bake one lane's
tree into a site that is meant to describe the merged one. Hand-editing
a generated file would be overwritten by the next real build. The right
move is one site build after the lanes land; this is flagged so it is
not forgotten.

## Merge with main 0ea59df, verified read-only

    my base c18253e is an ancestor of 0ea59df   yes
    files main changed since the base           51
    files this lane changed                     18
    overlap                                      1   _pkgdown.yml
    main's hunk in it                     line 18   navbar
    this lane's hunk in it                line 62   reference
    R/ files touched by both                     0
    tests/ files touched by both                 0

Clean. Nothing committed, per the lane's standing instruction.

## Punch-round reverification

### Gated tiers, all five files

    test-brms-likelihood.R          404 pass  0 fail  0 error  0 skip
    test-bcm-signal-detection.R      30 pass  0 fail  0 error  0 skip
    test-bcm-binomial.R              37 pass  0 fail  0 error  0 skip
    test-bcm-esp.R                   17 pass  0 fail  0 error  0 skip
    test-bcm-model-selection.R       29 pass  0 fail  0 error  0 skip

Every count matches the reviewer's independently measured figures.

The Stan identity residuals are BIT-IDENTICAL to the pre-punch run, all
eleven of them, which is the check that matters: the punch round
changed an error message, two documentation tables and a print
formatter, and nothing that enters a density. Had any of it touched
arithmetic, these would have moved.

    -2.593480986e-13   0   0   -6.217248938e-15   -1.243449788e-14
    -1.33226763e-14    -4.440892099e-16   0   -7.283063042e-14
    3.552713679e-15    -2.841282765e-12

Three further identities come from `test-bcm-model-selection.R`, which
is in the gated batch this round and was not in the last one:
-1.250555215e-12, -7.389644452e-13, 0. Fourteen identities in all,
every one holding.

### R CMD check --as-cran, punch round

    * checking top-level files ... OK
    * checking examples ... [41s] OK
    * checking examples with --run-donttest ... [71s] OK
    * checking tests ... [22m] OK
    Status: 1 NOTE          (no ERROR, no WARNING)

ONE NOTE, and it is the expected one:

    * checking HTML version of manual ... NOTE
      Skipping checking math rendering: package 'V8' unavailable

The examples NOTE of the previous run is GONE, which settles what it
was. Same tree, same examples: it listed four topics at 106s while
three jobs shared the machine, two topics at 60s with two jobs, and
does not fire at all at 41s with one. It was elapsed time under load,
never a property of the package, and none of the topics it named were
this lane's.

### Core suite, punch round, one file per process

    RESULT lines            131
    distinct files          131
    reported twice            0
    test files on disk      131
    files with no result      0
    files with a failure      0

    totals   pass=7721  fail=0  error=0  warn=1  skip=133

Reconciles to the previous round exactly: 7705 + 16 = 7721, and the 16
are the acat regression guard added to `test-ordinal.R` this round,
93 to 109. `test-families.R` 216 and `test-message-uniqueness.R` 6 are
unchanged.

The one warning is the same pre-existing `test-prior-compat.R`
convergence warning, still reproducing identically on the base.

Note on the reviewer's `pass=7691 error=5`: their single-file harness
runs outside `test_check()`, so `local_mocked_bindings()` has no
pkgload and five `test_that()` blocks abort in three files this lane
does not touch. This harness calls `pkgload::load_all()` per process,
so those five run, which is the whole of the difference and why this
count is 0 errors.
