#' @keywords internal
#'
#' @details
#' frmtmb fits regression models written in a 'brms'-style formula grammar
#' by maximum likelihood. Random effects are integrated out with the
#' Laplace approximation. Each model's objective is generated as an R
#' closure and differentiated on the 'RTMB' tape, so no compilation
#' happens at run time.
#'
#' Start with [frm()] to fit a model, [bf()] to write one, and
#' `vignette("frmtmb")` for a tour. `vignette("diagnostics")` covers
#' convergence, residuals, and the reliability of the approximation.
#'
#' @section Algorithm provenance:
#' frmtmb implements established methods; it introduces no new estimator.
#'
#' - The Laplace approximation to the marginal likelihood of a
#'   latent-variable model, and its automatic differentiation, come from
#'   'TMB' / 'RTMB' (Kristensen et al. 2016). frmtmb generates the
#'   objective; 'RTMB' differentiates and integrates it.
#' - The formula grammar follows brms (Buerkner 2017). frmtmb is a
#'   frequentist implementation of a documented subset of that grammar.
#'   Formulas port from brms without change.
#' - Random-effect syntax and the structured covariance vocabulary follow
#'   lme4 (Bates et al. 2015) and glmmTMB (Brooks et al. 2017).
#' - Smooth terms are built with mgcv (Wood 2017, 2011) through
#'   `smoothCon()` and `smooth2random()`.
#' - Approximate Gaussian processes use the Hilbert-space basis of
#'   Riutort-Mayol et al. (2023). Spatial Gaussian Markov random fields
#'   use the SPDE representation of Lindgren et al. (2011) and the BYM2
#'   parameterization of Riebler et al. (2016).
#' - One-step-ahead residuals come from Thygesen et al. (2017).
#' - Pooling across multiple imputations follows Rubin (1987) and the
#'   references given in [anova.frmtmb_multiple()].
#'
#' The contribution is architectural, not statistical. Existing
#' frequentist packages put a formula front end on one fixed likelihood;
#' frmtmb compiles the formula into the objective. This makes several
#' combinations ordinary code paths that are structural dead ends
#' elsewhere: random effects in any distributional parameter, nonlinear
#' predictors, per-response families in a multivariate model, monotonic
#' effects, in-model imputation, latent-class mixtures, and custom
#' families written as plain R log-densities. `vignette("compatibility")`
#' reports which combinations are supported, and which are refused.
#'
#' @section Life cycle:
#' frmtmb is maturing and is not yet on CRAN.
#'
#' The formula grammar is stable, because it follows brms. The fitted
#' object API is stable: [frm()], the accessor methods, and the family
#' constructors keep their current behavior. Fields of the fitted object
#' that no exported method reaches are internal and can change; use the
#' accessors.
#'
#' Multivariate coverage of the post-fit methods, the mixture families,
#' and [frm_sample()] still gain features between releases. Version
#' numbers stay below 1.0 until the CRAN release. Breaking changes are
#' listed in `NEWS.md`. The package is actively developed and maintained.
#'
#' @references
#' Bates, D., Maechler, M., Bolker, B., and Walker, S. (2015). Fitting
#' Linear Mixed-Effects Models Using lme4. *Journal of Statistical
#' Software*, 67(1), 1-48. \doi{10.18637/jss.v067.i01}
#'
#' Brooks, M. E., Kristensen, K., van Benthem, K. J., et al. (2017).
#' glmmTMB Balances Speed and Flexibility Among Packages for Zero-inflated
#' Generalized Linear Mixed Modeling. *The R Journal*, 9(2), 378-400.
#' \doi{10.32614/RJ-2017-066}
#'
#' Buerkner, P.-C. (2017). brms: An R Package for Bayesian Multilevel
#' Models Using Stan. *Journal of Statistical Software*, 80(1), 1-28.
#' \doi{10.18637/jss.v080.i01}
#'
#' Kristensen, K., Nielsen, A., Berg, C. W., Skaug, H., and Bell, B. M.
#' (2016). TMB: Automatic Differentiation and Laplace Approximation.
#' *Journal of Statistical Software*, 70(5), 1-21.
#' \doi{10.18637/jss.v070.i05}
#'
#' Lindgren, F., Rue, H., and Lindstrom, J. (2011). An explicit link
#' between Gaussian fields and Gaussian Markov random fields: the
#' stochastic partial differential equation approach. *Journal of the
#' Royal Statistical Society: Series B*, 73(4), 423-498.
#' \doi{10.1111/j.1467-9868.2011.00777.x}
#'
#' Riebler, A., Sorbye, S. H., Simpson, D., and Rue, H. (2016). An
#' intuitive Bayesian spatial model for disease mapping that accounts for
#' scaling. *Statistical Methods in Medical Research*, 25(4), 1145-1165.
#' \doi{10.1177/0962280216660421}
#'
#' Riutort-Mayol, G., Buerkner, P.-C., Andersen, M. R., Solin, A., and
#' Vehtari, A. (2023). Practical Hilbert space approximate Bayesian
#' Gaussian processes for probabilistic programming. *Statistics and
#' Computing*, 33(1), 17. \doi{10.1007/s11222-022-10167-2}
#'
#' Rubin, D. B. (1987). *Multiple Imputation for Nonresponse in Surveys*.
#' Wiley, New York. \doi{10.1002/9780470316696}
#'
#' Thygesen, U. H., Albertsen, C. M., Berg, C. W., Kristensen, K., and
#' Nielsen, A. (2017). Validation of ecological state space models using
#' the Laplace approximation. *Environmental and Ecological Statistics*,
#' 24(2), 317-339. \doi{10.1007/s10651-017-0372-4}
#'
#' Tierney, L., and Kadane, J. B. (1986). Accurate Approximations for
#' Posterior Moments and Marginal Densities. *Journal of the American
#' Statistical Association*, 81(393), 82-86.
#' \doi{10.1080/01621459.1986.10478240}
#'
#' Wood, S. N. (2011). Fast stable restricted maximum likelihood and
#' marginal likelihood estimation of semiparametric generalized linear
#' models. *Journal of the Royal Statistical Society: Series B*, 73(1),
#' 3-36. \doi{10.1111/j.1467-9868.2010.00749.x}
#'
#' Wood, S. N. (2017). *Generalized Additive Models: An Introduction with
#' R*, 2nd edition. Chapman and Hall/CRC.
#' \doi{10.1201/9781315370279}
#'
#' @srrstats {G1.0} Primary references for every implemented method are
#'   listed in the `@references` section above: Kristensen et al. (2016)
#'   for the Laplace approximation and AD backend, Tierney and Kadane
#'   (1986) for the approximation itself, Buerkner (2017) for the formula
#'   grammar, Bates et al. (2015) and Brooks et al. (2017) for the
#'   random-effect vocabulary, Wood (2011, 2017) for smooths,
#'   Riutort-Mayol et al. (2023) for approximate Gaussian processes,
#'   Lindgren et al. (2011) and Riebler et al. (2016) for the spatial
#'   fields, Thygesen et al. (2017) for one-step-ahead residuals, and
#'   Rubin (1987) for multiple-imputation pooling.
#' @srrstats {G1.1} The "Algorithm provenance" section states that
#'   frmtmb introduces no new estimator. It is a new implementation, in
#'   R, of established methods, and the contribution is architectural:
#'   the objective is generated from the formula rather than selected
#'   from a fixed set of compiled likelihoods. The section names the
#'   feature combinations this makes reachable that existing R
#'   implementations cannot express.
#' @srrstats {G1.2} The "Life cycle" section above states the current
#'   state (maturing, pre-CRAN), which parts of the interface are stable,
#'   which are still moving, and the maintenance intent. `CONTRIBUTING.md`
#'   and the README carry the same statement.
#' @srrstats {G1.4} All exported functions are documented with roxygen2.
#'   `NAMESPACE` holds 72 `export()` and 100 `S3method()` entries, and
#'   every one originates from a roxygen block; the package is built with
#'   `roxygen2::roxygenise()` and carries `Config/roxygen2/version` in
#'   `DESCRIPTION`.
#' @srrstats {G1.4a} All internal functions are documented with roxygen2
#'   blocks that state the contract (inputs, outputs, and the reason the
#'   helper exists) and end in `@noRd`. There are about 290 of them; the
#'   pass covering every internal function landed in v0.30.0.
#' @srrstats {G5.0} Tests use standard data sets with known properties
#'   wherever a reference implementation is being checked:
#'   `lme4::sleepstudy` and `lme4::cbpp`, `datasets::faithful`,
#'   `nlme::Soybean`, and `brms::inhaler` and `brms::kidney`. Models with
#'   no standard data set use seeded simulators shared through
#'   `tests/testthat/helper-reference.R`.
#' @srrstats {G5.10} The extended test tiers are switched on by
#'   environment variables, in the same testthat framework as the regular
#'   tests: `FRMTMB_FUZZ=true` runs the pairwise grammar fuzzer
#'   (`FRMTMB_FUZZ_N` caps the plan size) and `FRMTMB_BRMS_FIT_TESTS=true`
#'   runs the tier that compares against compiled brms Stan programs.
#'   `NOT_CRAN=true` runs the heavy reference-validation files gated by
#'   `skip_on_cran()`.
#' @srrstats {G5.12} `CONTRIBUTING.md` documents how to run the extended
#'   tests, which packages they need, and that they need no special
#'   hardware and download nothing.
#'
#' @importFrom stats anova coef confint cooks.distance deviance dfbeta
#'   dfbetas df.residual drop1 extractAIC family fitted formula
#'   influence logLik model.frame model.matrix na.action nobs predict
#'   profile residuals sigma simulate terms update vcov weights
#' @importFrom graphics plot
"_PACKAGE"
