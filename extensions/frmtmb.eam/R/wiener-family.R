#' The Wiener first-passage time family
#'
#' A drift-diffusion model for a two-choice decision: a noisy evidence
#' accumulator starts between two boundaries and the response time is
#' the first time it touches one of them. The family models the response
#' time; which boundary was touched is data, supplied through the
#' `dec()` addition term.
#'
#' The parameterization is brms's `wiener` family, name for name:
#'
#' \describe{
#'   \item{`mu`}{Drift rate, the mean rate of evidence accumulation
#'     (brms and the literature also call it `v`). Identity link, so it
#'     is signed: positive drift favors the upper boundary.}
#'   \item{`bs`}{Boundary separation, the distance between the two
#'     boundaries (`a`). Log link.}
#'   \item{`ndt`}{Non-decision time, the part of the response time spent
#'     encoding and moving rather than deciding (`t0` or `tau`). Bounded
#'     link; see Non-decision time below.}
#'   \item{`bias`}{Relative start point in (0, 1), the fraction of the
#'     boundary separation the accumulator starts at (`w`). Logit link.
#'     0.5 is unbiased.}
#' }
#'
#' `variability` adds Ratcliff's three across-trial variability
#' parameters to that set; see Across-trial variability below.
#'
#' @section The decision indicator:
#' The boundary a trial ended at is data, and it reaches the density as
#' an addition term:
#'
#' ```
#' frm(bf(rt | dec(response) ~ condition), family = wiener(), data = dat)
#' ```
#'
#' `dec()` is spelled as brms spells it and takes what brms takes: a
#' factor or character vector whose SECOND level is the upper boundary
#' (so `"lower"`/`"upper"` and `c(FALSE, TRUE)` both work as they read),
#' or a numeric 0/1 column. The package contributes the term to frmtmb's
#' addition-term registry when it loads.
#'
#' `vint()` also still works, and carries the indicator as a plain 0/1
#' integer column:
#'
#' ```
#' frm(bf(rt | vint(upper) ~ condition), family = wiener(), data = dat)
#' ```
#'
#' Supplying neither is refused with a message that says so, because the
#' failure is otherwise silent.
#'
#' @section The unit of the response:
#' SECONDS. Every default in this package assumes it: the starting
#' values, the bound the `ndt` link is scaled onto, and
#' `gddm_control(dt = )` and the window it takes from the data. Nothing
#' in any of the likelihoods refuses milliseconds, and a fit to
#' millisecond data converges and reports a boundary separation three
#' orders of magnitude out.
#'
#' So every family here warns when the fastest response in the data is
#' above 20, names milliseconds as the likely cause and says to divide
#' by 1000. It is a ceiling on the fastest response in the WHOLE data
#' set rather than on any one trial, so one slow trial does not reach
#' it.
#'
#' It is a warning rather than a refusal because it CAN fire on a
#' correct model: a slow task with a short session. Measured, the rate
#' is zero over 192 designs at the parameters a two-choice task
#' produces, and rises to 0.07 at 20 trials of a task whose median
#' response is 39 seconds, 0.56 at 60 seconds and 1.00 at 100 seconds.
#' A deliberation or matrix-reasoning design is where that lands. The
#' warning carries the class `frmtmb_eam_units_warning` so that such a
#' design can silence this one condition and keep the rest. `NEWS.md`
#' carries the full tables.
#'
#' @section Non-decision time, and what its coefficients mean:
#' The density is zero for a response time at or below `ndt`, so the
#' likelihood has a hard edge at the fastest response and an ordinary
#' log link would let the optimizer walk straight over it. The `ndt`
#' link makes the constraint structural instead, and it does so in one
#' of two ways.
#'
#' **Without `ndt_group()`, `ndt` is a TIME**, on a logit scaled onto
#' `(0, ub)` with `ub` the fastest response in the whole data set, or
#' `max_ndt` when you give one. This is the parameterization the family
#' has always had. `frm_linpred(dpar = "ndt", type = "response")` reports
#' seconds, a `prior(class = "ndt")` is a density on those seconds, and
#' a `bf(ndt = 0.2)` constant is 0.2 seconds.
#'
#' **With `ndt_group()`, `ndt` is a FRACTION of the row's own bound**,
#' on a plain logit, and the density multiplies it by that bound. The
#' bound is the fastest response of the row's group. Write
#'
#' ```
#' frm(bf(rt | dec(response) + ndt_group(subject) ~ coherence,
#'        ndt ~ 1 + (1 | subject), bias = 0.5),
#'     family = wiener(), data = dat)
#' ```
#'
#' and each subject's non-decision time is bounded by its own fastest
#' response. That is what a random effect on `ndt` needs; the next
#' section is what one global bound does to it. The price is that `ndt`
#' is on a different scale: `frm_linpred(dpar = "ndt", type = "response")`
#' reports the fraction, and [ndt_time()] reports the time for either
#' parameterization. A `prior(class = "ndt")` and a `bf(ndt = )`
#' constant are fractions under a grouping too.
#'
#' Give the grouping as a factor, a character vector, a logical or
#' integer codes. It is keyed on the group's LABEL, so subsetting,
#' `droplevels()`, `relevel()` and a prediction grid you build yourself
#' all pair a row with the same bound the fit used. `max_ndt` and
#' `ndt_group()` cannot be combined, because they set the same bound to
#' different things, and an `ndt_group()` no family reads is refused
#' rather than carried into the fit unused.
#'
#' A `max_ndt` above the smallest response time is refused, because for
#' this family alone it admits parameter values at which some observed
#' row has no likelihood. Inside a [frmtmb::mixture()] that is exactly
#' what the other component is for, so `allow_unreachable = TRUE` lifts
#' the refusal; see Mixtures. A mixture never finalizes its components,
#' so a component needs `max_ndt` given up front and cannot use
#' `ndt_group()`.
#'
#' Whichever parameterization a model is in, the bound is fixed when
#' the model frame is assembled and is a property of the FITTED data,
#' so a prediction on new rows is scaled by the bound the fit used
#' rather than by one re-derived from the new rows. A group the fit
#' never saw has no bound and is refused.
#'
#' Past a linear predictor of about 37 the logit saturates in double
#' precision. Nothing guards that, and nothing needs to: the density
#' falls off a cliff as the decision time goes to zero, so the log
#' likelihood is already unreachable long before the link runs out of
#' digits.
#'
#' @section Why one bound is the wrong constraint under a random effect:
#' A single bound is the GLOBAL fastest response, so a subject whose
#' true `ndt` is above it cannot be represented at any value of the
#' random effect. At 30 subjects by 400 trials with a between-subject
#' spread of 26 ms on a mean of 250 ms, 20 of the 30 subjects are in
#' that position and NONE is inconsistent with its own data. Without
#' `variability` such a fit does not converge; with
#' `variability = "sv"` it converges, reports nothing from `diagnose()`,
#' and returns a population `ndt` pinned at the bound and wrong by ten
#' percent with a standard error of 7.2e-06 on it.
#'
#' What settles it is what happens as data accumulates. On that design,
#' the per-subject root mean squared error of the fitted non-decision
#' times:
#'
#' \tabular{rll}{
#'   trials \tab with ndt_group(s) \tab one global bound \cr
#'   100 \tab 20.91 ms \tab 27.23 ms \cr
#'   200 \tab 14.08 ms \tab 20.51 ms \cr
#'   400 \tab 7.67 ms \tab 30.37 ms
#' }
#'
#' The per-group bound converges on the truth and the global bound does
#' not, because more data lowers the global minimum and tightens the
#' ceiling on every subject at once.
#'
#' @section What sd(ndt) does and does not tell you:
#' A grouped model estimates `ndt` as a fraction of each group's own
#' floor, so the fitted per-subject times vary with those floors even
#' when the variance component is exactly zero. On the design above an
#' estimator with NO random effect on `ndt` returns a between-subject
#' standard deviation of 0.02748 against a truth of 0.02629, where the
#' full model returns 0.02543, and at 100 trials per subject the two
#' are identical in every digit.
#'
#' The variance component is not empty: at 400 trials it buys 8.44
#' log-likelihood units and cuts the per-subject error from 11.90 ms to
#' 7.67 ms. But a between-subject standard deviation is the wrong
#' statistic to read that off. Compare the per-subject non-decision
#' times from [ndt_time()] against what you believe, and compare
#' log-likelihoods, rather than reading `VarCorr()`'s `ndt` entry as
#' evidence that the component was estimated.
#'
#' The bound's own quality is a function of the group's trial count.
#' The group's fastest response overshoots its true non-decision time
#' by about 47 ms at 400 trials and 78 ms at 50, and a group of one
#' trial gets that trial's own response time as its bound. Nothing
#' refuses a small group, because any threshold would fire on correct
#' models; the trial counts are recorded on the fitted family, at
#' `family(fit)$ndt_bound$sizes`.
#'
#' @section Recovery at 30 subjects by 400 trials:
#' Measured on 60 replicate data sets at the design named in the
#' realistic-scale table of `dev/extension-gaps-plan.md`: 30 subjects,
#' 400 trials each, two conditions, `mu ~ cond + (1 | s)`,
#' `bs ~ 1 + (1 | s)` and `ndt ~ 1 + (1 | s)` with `ndt_group(s)`. The
#' truths are a drift of 0.4 and 1.3, a boundary separation of 1.4, a
#' non-decision time of 250 ms, and between-subject standard deviations
#' of 0.35 on the drift, 0.20 on the log boundary and 0.12 on the log
#' non-decision time. Seeds 20260910 to 20260969.
#' `dev/eamhier-findings.md` holds the tables and
#' `dev/eamhier-scripts/` the harness.
#'
#' All 60 fits converged with a positive definite Hessian, no `NaN`
#' standard error and nothing from `diagnose()`, and every subject's
#' fitted non-decision time stayed below its own fastest response.
#'
#' \preformatted{
#'   quantity              truth    mean     mcse     Wald coverage
#'   mu intercept          0.4      0.4019   0.0094   54/60   90.0\%
#'   mu condition effect   0.9      0.9004   0.0035   57/60   95.0\%
#'   bs                    1.4      1.3900   0.0067   56/60   93.3\%
#'   sd(mu | s)            0.35     0.3474   0.0059   55/60   91.7\%
#'   sd(log bs | s)        0.20     0.1898   0.0037   54/60   90.0\%
#' }
#'
#' **What 60 replicates can resolve.** The rule applied is that a
#' coverage misses when its Wilson interval excludes 0.95, which at
#' n = 60 rejects on 53 or fewer. Against a nominal rate that is known
#' rather than estimated, the exact binomial power of that rule is 0.82
#' at a true 85 percent, 0.59 at 88 and 0.39 at 90, with a size of
#' 0.030. So the table can show a ten-point shortfall and cannot show
#' seven. No coefficient's interval excludes 95.
#'
#' **Both variance components come back about as low as maximum
#' likelihood is expected to put them.** The shrinkage at 30 groups,
#' including the term for reporting a standard deviation rather than a
#' variance, is `sqrt(1 - 1/30) * (1 - 1/116) = 0.9747`, that is 2.5
#' percent low. Against that expectation `sd(mu | s)` sits 1.06 Monte
#' Carlo standard errors high and `sd(log bs | s)` 1.39 low: one common
#' shrinkage covers both, and neither is a separate finding. Read a
#' variance component from 30 groups as a slight underestimate.
#'
#' The non-decision time is read per subject rather than as a spread,
#' for the reason the section above gives. Its root mean squared error
#' is 7.9 ms against a between-subject spread of 30.2 ms, it is smaller
#' than that spread on 60 of 60 replicates, and its correlation with
#' the drawn values is 0.967. The population value comes back about
#' 2 ms high.
#'
#' The same 60 draws with `variability = "sv"` recover the same
#' quantities to the same accuracy. What they add is in the next
#' section.
#'
#' @section Across-trial variability:
#' Ratcliff's full diffusion model draws three of the four parameters
#' afresh on every trial. `variability` names which of those to
#' estimate, and each one it names becomes an ordinary distributional
#' parameter that takes its own formula:
#'
#' ```
#' frm(bf(rt | dec(response) ~ coherence, bias = 0.5),
#'     family = wiener(variability = c("sv", "sz", "st")), data = dat)
#' ```
#'
#' \describe{
#'   \item{`sv`}{Standard deviation of a normal drift rate. Log link.}
#'   \item{`sz`}{Width of a uniform relative start point, centered on
#'     `bias`, on the same (0, 1) scale as `bias`. Logit link, so the
#'     width is below 1 and the start point stays inside the boundaries
#'     whenever `bias` is 0.5.}
#'   \item{`st`}{Width of a uniform non-decision time, centered on
#'     `ndt`, and on whichever scale `ndt` is on: a duration in the
#'     units of the response, on a logit scaled onto `(0, 2 * bound)`,
#'     or a FRACTION of twice the row's own bound under
#'     `ndt_group()`.}
#' }
#'
#' The likelihood is the analytic Wiener density averaged over those
#' distributions, and the three are done three different ways because
#' they are three different integrals. The drift integral is Gaussian
#' against an exponential-quadratic and is evaluated in CLOSED FORM: it
#' is exact, it takes no nodes, and there is nothing to tune. The other
#' two are uniform and are evaluated by fixed-node Gauss-Legendre
#' quadrature, whose node counts are the `nodes` argument.
#'
#' The node positions and counts are decided when the family object is
#' built and are constants from then on, because a node count that moved
#' with a parameter would be a branch on a parameter and an
#' automatic-differentiation tape cannot record one. A parameter only
#' rescales the interval the fixed nodes are mapped onto.
#'
#' `fitted()` and `simulate()` follow the variability rather than
#' ignoring it. Both condition on the boundary a row ended at, and
#' conditioning reweights which per-trial parameters that row could have
#' had: the trials that reached a boundary are not a fair sample of the
#' drift rates that could have produced them. So the fitted mean is a
#' ratio of two quadratures and the simulator accepts a drawn drift rate
#' and start point with the boundary probability they imply. The plain
#' closed-form mean is not a usable approximation here: at an unbiased
#' start point it returns the same number for both boundaries, and the
#' full model does not.
#'
#' Two limits of the parameterization are worth knowing. The uniform
#' start point stays inside the boundaries by construction only when
#' `bias` is 0.5, which is the usual case and the one the logit link on
#' `sz` is scaled for; at a strongly biased start a wide `sz` can push
#' the range past a boundary, where the density is near zero and the
#' likelihood is a barrier rather than a cliff. And nothing holds
#' `ndt - st / 2` above zero, so a fit is free to report a range that
#' includes a negative non-decision time. Neither can be made structural
#' from outside frmtmb: both are joint constraints on two distributional
#' parameters, and a link is a property of one.
#'
#' **How well `sv` is identified, measured.** On 60 replicates of one
#' subject with 12,000 trials, no random effects, and a true `sv` of
#' 0.4, the estimate ranges from 0.194 to 0.599 and its Wald interval
#' covers on 57 of the 59 fits that report a usable one. The estimated
#' `sv` and an estimated drift contrast are correlated, 0.59 by rank: a
#' draw that returns a low `sv` returns a low contrast with it, and
#' every replicate that missed on the contrast had an `sv` below the
#' truth. Read the two together rather than one at a time. In the
#' hierarchical design of "Recovery at 30 subjects by 400 trials" both
#' cover at 57 of 60.
#'
#' **One fit in 60 lost `sv` entirely and nothing reported it.** At one
#' seed the estimate ran to the log link's floor, 4.4e-04, with a Wald
#' interval of (-1163, 1147) on the log scale, convergence code 0, a
#' positive definite Hessian and `diagnose()` printing "No convergence
#' problems detected". The likelihood is genuinely flat there on that
#' draw: dropping the `sv` term moves the log-likelihood by 7e-07,
#' while holding `sv` at its true value costs 4.05 units. So the fit is
#' not wrong, but the interval it prints is not a statement about
#' anything. Read the interval `confint()` gives `sv` before you read
#' the estimate. The signature is inside that one fit and needs no
#' replicates: the standard error on `log sv` was 589.25 where the
#' largest of the other four was 0.0501, a ratio of 11,763. Nothing
#' else in the fit reports it.
#'
#' @section Mixtures:
#' A contaminant component covers the trials the diffusion process
#' cannot produce, which is the standard treatment for fast guesses. The
#' Wiener component then wants a non-decision time that some observed
#' rows fall below, so pass the bound and lift the refusal:
#'
#' ```
#' frm(bf(rt | dec(response) ~ 1, bias1 = 0.5),
#'     family = mixture(wiener(max_ndt = 0.4, allow_unreachable = TRUE),
#'                      lognormal()),
#'     data = dat)
#' ```
#'
#' A row below the non-decision time gets a log density of about
#' `-1 / delta` where `delta` is a billionth of the smallest response
#' time: it exponentiates to exactly zero, which is the right likelihood
#' for the component, and it differentiates to exactly zero, which a
#' true `-Inf` would not.
#'
#' @section Censoring and truncation:
#' `cens()` and `trunc()` work on the response time, on the plain model.
#' What a censored row means depends on its code, because the boundary
#' is known on some and not on others:
#'
#' \describe{
#'   \item{right}{The trial had not reached EITHER boundary by the
#'     recorded time: a deadline passed with no response. Its
#'     probability is the survival over both boundaries,
#'     `P(T > t)`. Its `dec()` value is required, since a declaration
#'     cannot depend on a censoring code, and is NOT read: any 0 or 1
#'     gives the same fit, bit for bit.}
#'   \item{left and interval}{The trial DID reach a boundary, and
#'     `dec()` says which; only its time is coarse. Its probability is
#'     that boundary's defective distribution function,
#'     `P(T <= t, boundary b)`, or its difference across the interval.
#'     `dec()` is read. A row whose `dec()` is missing is dropped by
#'     frmtmb's `na.action` with a message, before this family sees it;
#'     this family cannot refuse it by name.}
#' }
#'
#' A deadline design is
#'
#' ```
#' frm(bf(rt | dec(response) + cens(censored) ~ condition, bias = 0.5),
#'     family = wiener(), data = dat)
#' ```
#'
#' with `censored` set to `"right"` on the trials with no response and
#' `rt` set to the deadline there.
#'
#' `trunc()` divides by the window's mass over both boundaries. Left or
#' interval censoring cannot be combined with `trunc()`: frmtmb forms
#' the censored rows and the truncation window from the one
#' distribution-function slot a family has, and the two need different
#' functions here, so the pair is refused by name.
#'
#' Every quantity is computed as a LOG and never as a complement: the
#' survival from the method of images for small times and from the
#' eigenfunction series for large ones, and each defective function the
#' same way. Against a 700-bit reference over drift -20 to 20, boundary
#' separation 0.3 to 6, relative start 0.001 to 0.999 and normalized
#' time 1e-3 to 10, 3150 points, the worst relative error is 4.3e-14 on
#' the distribution function over both boundaries, 1.2e-12 on either
#' defective one, and on the survival 1.7e-12 where |v| a is at most 1,
#' 7.4e-12 up to 24, 4.6e-11 up to 72 and 8.0e-11 up to 120. These
#' hold on that grid only. The review of 2026-09-24 measured 798 more
#' points at 1200 bits, to |v| a = 250 and a relative start of 0.9999,
#' and there the survival's relative error is 1.7e-11 where |v| a is
#' at most 1, 1.5e-9 at 120 and 1.6e-9 at 250, the worst at relative
#' start 0.9999 and normalized time 0.01. 1.6e-9 is the worst measured
#' anywhere. A survival
#' of exp(-72063) comes back as that log. `RWiener::pwiener()` agrees to
#' 4.5e-11 absolute on the grid the density is pinned on, and all of
#' that difference is RWiener's own error.
#'
#' An interval's mass is formed as a difference of logs, never of
#' probabilities: F(t2) (1 - F(t1) / F(t2)), or the same with the mass
#' still to come after each edge, whichever keeps more digits, or a
#' quadrature of the density where the interval holds a tiny share of
#' the mass on both sides. Against 751 intervals at 1200 bits, 511 of
#' them the review's, with the lower edge at normalized time 1e-3 to 3,
#' |v| a up to 80 and masses down to 1e-99 of F(t1), the worst error
#' of the log mass is 4.9e-11. The probability of a
#' left-censored row, or of an interval, is held at 1e-300 or above,
#' because frmtmb takes it on the probability scale: a row whose true
#' probability is smaller than that contributes a flat barrier rather
#' than its value.
#'
#' **What it recovers, at 30 subjects by 400 trials** (drift 0.4 and 0.9
#' by condition with a subject deviation of 0.35, boundary separation
#' 1.4 with a subject deviation of 0.20 on its log, non-decision time
#' 0.25 s). With a 1 s deadline, 16 percent of trials right-censored,
#' 40 replicates: every fixed effect and both variance components
#' covered on 33 to 38 of 38 fits with an interval. With trials faster
#' than 0.45 s recorded as left-censored and every fifth other trial as
#' interval-censored into its 100 ms bin, 29 and 14 percent of trials,
#' 60 replicates: they covered on 47 to 55 of 55, the lowest being the
#' drift intercept at 85.5 percent, Wilson interval 73.8 to 92.4.
#' Every fit converged. The arms draw from the same seeds, so their
#' coverages are not independent of each other.
#'
#' The plain model only. Under `variability` the distribution function
#' is the same series averaged over the per-trial parameters, and the
#' drift average has no closed form there, so the pair is refused.
#'
#' @section Contaminant trials:
#' `contaminant = TRUE` mixes a second process into every row: a
#' response time uniform over a window, with a boundary that is a coin
#' flip. Its share is the distributional parameter `lambda`, on a logit
#' link, which takes a formula like any other:
#'
#' ```
#' frm(bf(rt | dec(response) ~ condition, lambda ~ 1, bias = 0.5),
#'     family = wiener(contaminant = TRUE, contaminant_range = c(0, 3)),
#'     data = dat)
#' ```
#'
#' It needs no [frmtmb::mixture()]: the uniform has no free parameter,
#' so it is one more term in the density rather than a second family. A
#' right-censored row enters through both parts' survivals; a left- or
#' interval-censored row through the diffusion's defective function at
#' its boundary and half of the uniform's, the coin flip landing there.
#'
#' **The window.** Give it as `contaminant_range =`: the task's
#' response window, such as `c(0, deadline)`. Without it, a model that
#' declares its deadline through `trunc(ub = )` uses the fastest
#' response to that deadline, and any other model is refused, with the
#' class `frmtmb_eam_contaminant_range_error`. A response outside the
#' window is refused with the same class, since the contaminant has
#' density zero there. The refusal of a missing window is the
#' measured reason: without a deadline bounding both processes, the
#' slowest diffusion trial sets the top of the observed range, the
#' window grows with the sample, and `lambda` collapses. On one subject
#' of 4000 trials, 25 seeds, true share 0.0495, contaminants uniform on
#' 0 to 5 s: the observed range gave `lambda` 0.0212 and covered on 2 of
#' 25, the true window 0.0498 and 25 of 25. With a 5 s deadline declared
#' by `trunc(ub = 5)`, the window from the fastest response to the
#' deadline gave 0.0508 against a recorded share of 0.0496, and the one
#' to the slowest response 0.0518, both covering on 25 of 25; the
#' deadline is the top used, being the one that does not move with the
#' sample (`dev/phase3b-cont-window.R`). The review of 2026-09-24
#' measured the same pattern on its own design: 0.0050 and 5 of 25
#' without a deadline, 0.0545 and 22 of 25 with one.
#'
#' Where the practice comes from. Ratcliff and Tuerlinckx (2002) model
#' the contaminant as uniform over the range of the observed response
#' times, and the DMAT toolbox follows them; HDDM uses a uniform of
#' fixed density instead. None of the three was checked against its
#' source or text for this page; the attributions are the review's of
#' 2026-09-24, recorded as such.
#'
#' **What it recovers, at 30 subjects by 400 trials.** Drift 0.4 and
#' 0.9 by condition with a subject deviation of 0.35, boundary
#' separation 1.4 with a subject deviation of 0.20 on its log,
#' non-decision time 0.25 s, 5 percent contaminants uniform on 0.1 to
#' 5 s, `max_ndt = 0.5`. With `contaminant_range = c(0.1, 5)`, over 40
#' replicates, `lambda` came back at 0.0499 and its Wald interval
#' covered on 39 of 40; every fixed effect and both variance components
#' covered on 34 to 40 of 40, and every fit converged. The same draws
#' with a 3 s deadline declared by `trunc(ub = 3)` and the default
#' window, 44 replicates: `lambda` 0.0306 against a recorded share of
#' 0.0307, covering on 35 of 36. On data with no contaminant and the
#' window given, 30 replicates: see the collapse case below.
#'
#' **A guess faster than the non-decision time.** The non-decision time
#' is bounded by the fastest response, so a contaminant faster than the
#' true non-decision time pulls the bound under it and the fit cannot
#' reach the truth. On the design above at the default bound, 20
#' replicates on the first round's build, the non-decision time came
#' back at 0.111 s against 0.25 and the drift effect covered on 0 of
#' 20. Under the contaminant such a row
#' has a likelihood, so `max_ndt` above the fastest response is allowed
#' with `contaminant = TRUE`: give a bound you know the non-decision time
#' is under.
#'
#' **The collapse case.** Maximum likelihood has no prior to hold
#' `lambda` off its edge. When the data give it nothing, it runs down
#' the logit: the estimate is a logit near -20, its standard error is in
#' the thousands, and the log-likelihood is the plain family's.
#' [frmtmb::diagnose()] reports it as a parameter at the end of its
#' link. On the design above with no contaminant and the window
#' `c(0, 5)` given, 30 replicates: `lambda` reached the edge on 25,
#' `diagnose()` named it on 27, all 25 at the edge among them, and at
#' the edge the log-likelihood was the plain family's to within
#' 8.9e-07 on every one.
#'
#' At `lambda = 0` the density is the plain family's: bit for bit on
#' every row whose diffusion density is at least the contaminant's, and
#' to the rounding of one log-sum anchor on the others.
#'
#' @section Accuracy:
#' The density is the Navarro and Fuss (2009) pair of series, both
#' evaluated at a fixed truncation and combined with a smooth weight,
#' because an automatic-differentiation tape cannot choose between them
#' on a parameter. It agrees with `RWiener::dwiener()` to better than
#' 1e-12 relative on the log scale over normalized times from 1e-3 to
#' 50. See `vignette("ddm")`.
#'
#' @param max_ndt Upper bound for the non-decision time, in the units of
#'   the response, applied to every row. `NULL`, the default, takes the
#'   fastest response of each row's `ndt_group()`, or of the whole data
#'   set when the model has no `ndt_group()`. It cannot be combined with
#'   `ndt_group()`. Give it when a component of a [frmtmb::mixture()]
#'   needs the bound up front.
#' @param variability Which across-trial variability parameters to
#'   estimate: any of `"sv"` (drift rate), `"sz"` (start point) and
#'   `"st"` (non-decision time). The default estimates none, which is
#'   the plain Wiener model.
#' @param nodes Gauss-Legendre node counts for the two quadratures, as a
#'   named vector. Only the entries for the `variability` parameters in
#'   use are read. The defaults are measured rather than chosen: `sz`
#'   reaches 1e-10 in the log density by 7 nodes everywhere it was
#'   probed, and `st` needs more because a response time below
#'   `ndt + st / 2` cuts the range and the integrand turns on sharply at
#'   the cut. See `vignette("ddm")` for the measurement.
#' @param allow_unreachable Permit a `max_ndt` above the smallest
#'   response time. Only correct inside a mixture, where another
#'   component carries the rows the Wiener density cannot reach.
#' @param link Link for the drift rate. Identity by default, and there
#'   is rarely a reason to change it: the drift rate is signed.
#' @param contaminant Mix a uniform contaminant into the density, with
#'   its mixing proportion estimated as the distributional parameter
#'   `lambda` on a logit link. `FALSE`, the default, is the plain model.
#'   See Contaminant trials.
#' @param contaminant_range The response times the uniform contaminant
#'   spreads over, as `c(lower, upper)` in the units of the response:
#'   the task's response window. `NULL`, the default, is allowed only
#'   when the model declares a deadline with `trunc(ub = )`, and is then
#'   the fastest response to that deadline; otherwise it is refused. See
#'   Contaminant trials for why.
#'
#' @return A `frmtmb_family`.
#'
#' @references
#' Navarro, D. J. and Fuss, I. G. (2009). Fast and accurate calculations
#' for first-passage times in Wiener diffusion models. *Journal of
#' Mathematical Psychology*, 53(4), 222-230.
#'
#' Ratcliff, R. and Tuerlinckx, F. (2002). Estimating parameters of the
#' diffusion model: approaches to dealing with contaminant reaction
#' times and parameter variability. *Psychonomic Bulletin & Review*,
#' 9(3), 438-481.
#'
#' @examples
#' set.seed(1)
#' dat <- ddm_simulate(300, mu = 0.8, bs = 1.4, ndt = 0.3, bias = 0.5)
#' fit <- frm(bf(rt | dec(upper) ~ 1, bias = 0.5),
#'            family = wiener(), data = dat)
#' fixef(fit)
#'
#' @export
wiener <- function(max_ndt = NULL, variability = character(0),
                   nodes = c(sz = 7L, st = 21L),
                   allow_unreachable = FALSE, link = "identity",
                   contaminant = FALSE, contaminant_range = NULL) {
  if (!is.logical(contaminant) || length(contaminant) != 1L ||
      is.na(contaminant)) {
    frm_stop("wiener(): `contaminant` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(contaminant_range)) {
    if (!contaminant) {
      frm_stop("wiener(): `contaminant_range` is the range of the ",
               "contaminant, and there is none without contaminant = TRUE.",
               call. = FALSE)
    }
    if (!is.numeric(contaminant_range) || length(contaminant_range) != 2L ||
        any(!is.finite(contaminant_range)) || contaminant_range[1L] < 0 ||
        !(contaminant_range[2L] > contaminant_range[1L])) {
      frm_stop("wiener(): `contaminant_range` must be two finite response ",
               "times, the first at least 0 and below the second.",
               call. = FALSE)
    }
  }
  if (!is.null(max_ndt)) {
    if (!is.numeric(max_ndt) || length(max_ndt) != 1L ||
        !is.finite(max_ndt) || max_ndt <= 0) {
      frm_stop("wiener(): `max_ndt` must be one positive finite number, ",
               "or NULL to take it from the data.", call. = FALSE)
    }
  }
  if (!is.logical(allow_unreachable) || length(allow_unreachable) != 1L ||
      is.na(allow_unreachable)) {
    frm_stop("wiener(): `allow_unreachable` must be TRUE or FALSE.",
             call. = FALSE)
  }
  if (contaminant && allow_unreachable) {
    # allow_unreachable exists for a component of frmtmb::mixture(), and
    # a mixture never finalizes its components, so the contaminant's
    # range would never be set; and the pair would be two contaminants
    # for one set of fast guesses
    frm_stop("wiener(): `contaminant = TRUE` and `allow_unreachable = ",
             "TRUE` cannot be combined. allow_unreachable is for a ",
             "component of mixture(), where the other component is the ",
             "contaminant; contaminant = TRUE is the contaminant inside ",
             "the family, and needs no mixture.", call. = FALSE)
  }
  cfg <- list(max_ndt = max_ndt, link = link,
              allow_unreachable = isTRUE(allow_unreachable),
              variability = ddm_check_variability(variability),
              nodes = ddm_check_nodes(nodes),
              contaminant = isTRUE(contaminant),
              contaminant_range = if (is.null(contaminant_range)) NULL else
                as.numeric(contaminant_range))
  # `ndt` is a fraction of a bound the data settles, so the family
  # carries a bound from the moment it is built: the one `max_ndt`
  # names, or a refusing placeholder. See ddm_ndt_preinstall().
  ddm_ndt_preinstall(ddm_family(cfg, delta = 1e-9), max_ndt, "wiener")
}

#' The variability names, in the order the dpars are declared.
#'
#' Canonical order rather than the user's, so that two spellings of the
#' same model give the same parameter vector and the same summary.
#'
#' `what` carries the family name into the refusal, because [wiener()]
#' and [wiener_gng()] offer the same three parameters under the same
#' argument and a user should be told which of the two refused.
#'
#' @noRd
ddm_check_variability <- function(variability, what = "wiener") {
  known <- c("sv", "sz", "st")
  if (is.null(variability)) variability <- character(0)
  if (!is.character(variability) || anyNA(variability) ||
      !all(variability %in% known) || anyDuplicated(variability)) {
    frm_stop(what, "(): `variability` names the across-trial variability ",
             "parameters to estimate, as a character vector with no ",
             "repeats, drawn from \"sv\" (drift rate), \"sz\" (start ",
             "point) and \"st\" (non-decision time).", call. = FALSE)
  }
  known[known %in% variability]
}

#' Merge a user's node counts over the defaults.
#'
#' @noRd
ddm_check_nodes <- function(nodes) {
  out <- c(sz = 7L, st = 21L)
  if (is.null(nodes) || !length(nodes)) return(out)
  if (!is.numeric(nodes) || is.null(names(nodes)) ||
      !all(names(nodes) %in% names(out)) || anyDuplicated(names(nodes)) ||
      any(!is.finite(nodes)) || any(nodes < 1) ||
      any(nodes != round(nodes))) {
    frm_stop("wiener(): `nodes` gives the Gauss-Legendre node counts as a ",
             "named vector of whole numbers at least 1, with names drawn ",
             "from \"sz\" and \"st\".", call. = FALSE)
  }
  out[names(nodes)] <- as.integer(nodes)
  out
}

`%||%` <- function(x, y) if (is.null(x)) y else x

#' Build the family object, given whatever the data has settled.
#'
#' Called twice: once by [wiener()], before any data exists, and once
#' more from `family_finalize()` with the unreachable-row margin the
#' response determines. The second call is what the fit actually uses.
#'
#' Every density, mean and simulator below is written in TIMES.
#' `ddm_ndt_install()` wraps them so that each receives `ndt` and `st`
#' multiplied by the row's own bound; nothing here knows that the
#' estimated quantity is a fraction.
#'
#' @noRd
ddm_family <- function(cfg, delta, crange = c(NA_real_, NA_real_)) {
  vv <- cfg$variability
  nd <- ddm_nodes(vv, cfg$nodes)
  st_on <- "st" %in% vv
  dpars <- c("mu", "bs", "ndt", "bias", vv)
  # Both `ndt` links are PLACEHOLDERS. ddm_ndt_install() replaces them
  # once the bound is known: with a scaled logit carrying that bound
  # when there is one number, and with a plain logit on a fraction when
  # ndt_group() makes it per row. See the head of ddm-shared.R.
  links <- list(mu = cfg$link, bs = "log", ndt = "logit", bias = "logit")
  if ("sv" %in% vv) links$sv <- "log"
  # a width on the same (0, 1) scale as bias, so the logit is the
  # scaled logit its own support asks for
  if ("sz" %in% vv) links$sz <- "logit"
  # and a duration in the response's units, bounded by the same
  # quantity that bounds the non-decision time it is centered on: a log
  # link here lets the optimizer walk out to a width no response time
  # could have come from, where every row's range is cut and the
  # surface is flat
  if (st_on) links$st <- "logit"

  lpdf <- if (!length(vv) && !cfg$allow_unreachable) {
    # The plain Wiener density, untouched: no variability parameter
    # exists, so there is nothing to average over and no quadrature to
    # pay for, and the bounded link guarantees every row is reachable,
    # so there is nothing to hold off the singularity either.
    function(y, dpars, aterms) {
      ddm_lpdf_both(y - dpars[["ndt"]], dpars[["mu"]], dpars[["bs"]],
                    dpars[["bias"]], ddm_indicator(aterms))
    }
  } else if (!length(vv)) {
    # The same density with the decision time held at `delta`. Only
    # reached when the user has declared that some rows are unreachable,
    # which is the mixture case: there the Wiener component's own
    # likelihood for a fast guess is zero, and it has to be a zero the
    # log-sum-exp can differentiate. -Inf is not: it exponentiates to
    # zero correctly and then contributes NaN to every gradient.
    function(y, dpars, aterms) {
      ddm_lpdf_both(ddm_floor(y - dpars[["ndt"]] - delta, delta),
                    dpars[["mu"]], dpars[["bs"]], dpars[["bias"]],
                    ddm_indicator(aterms))
    }
  } else {
    function(y, dpars, aterms) {
      ddm_lpdf_var(y, dpars[["mu"]], dpars[["bs"]], dpars[["bias"]],
                   dpars[["ndt"]],
                   if (is.null(dpars[["sv"]])) 0 else dpars[["sv"]],
                   if (is.null(dpars[["sz"]])) 0 else dpars[["sz"]],
                   if (is.null(dpars[["st"]])) 0 else dpars[["st"]],
                   ddm_indicator(aterms), nd, st_on, delta)
    }
  }

  init <- list(
    # A drift of zero is the honest starting guess: the sign of the
    # drift is what the data are there to tell us, and a start with
    # the wrong sign costs more than a start at the middle.
    mu = function(y, aterms) 0,
    bs = function(y, aterms) 1.5,
    # a placeholder; ddm_ndt_install() sets the real one, which is
    # half the fastest response or half the group's own bound
    ndt = function(y, aterms) 0.5 * min(y),
    bias = function(y, aterms) 0.5)
  # Small starts for the variability parameters, because zero is on the
  # boundary of every one of their links and a large start makes the
  # first quadrature straddle a range the data cannot support.
  if ("sv" %in% vv) init$sv <- function(y, aterms) 0.3
  if ("sz" %in% vv) init$sz <- function(y, aterms) 0.05
  # a placeholder, as ndt's is
  if (st_on) init$st <- function(y, aterms) 0.1 * min(y)

  # The distribution function, for cens() and trunc(). The plain model
  # only: under variability it is the same series averaged over the
  # three per-trial distributions, and the drift average is not in
  # closed form there (see the variability section of wiener-cdf.R).
  # A variability model is refused by name, in ddm_check_response(). The
  # slots are still filled for it, with the refusal, because frmtmb
  # checks that a censored model's family HAS a distribution function
  # before it reads the family's own validator, and its message there
  # would name neither this family nor the reason.
  lcdf <- function(q, dpars, aterms) ddm_stop_cens_var()
  lccdf <- function(q, dpars, aterms) ddm_stop_cens_var()
  if (!length(vv)) {
    # Over both boundaries, except on a LEFT- or INTERVAL-censored row:
    # that trial reached a boundary and dec() says which, so its
    # probability is that boundary's defective distribution function.
    # A right-censored row reached none and is scored from lccdf, over
    # both. The censoring code is data, so the choice is a mask.
    lcdf_rows <- function(q, dpars, aterms) {
      lF <- ddm_rt_lcdf2(q - dpars[["ndt"]], dpars[["mu"]], dpars[["bs"]],
                         dpars[["bias"]])[["lF"]]
      cen <- aterms[["cens"]]
      kb <- if (is.null(cen)) 0 else as.numeric(cen == -1 | cen == 2)
      if (!any(kb == 1)) return(exp(lF))
      lb <- ddm_rt_lcdf_b(q - dpars[["ndt"]], dpars[["mu"]], dpars[["bs"]],
                          dpars[["bias"]], ddm_indicator(aterms))
      # frmtmb takes this slot on the probability scale and forms
      # log(F) for a left-censored row and log(F(y2) - F(y)) for an
      # interval. Two things break there, and each stopped every fit of
      # the first recovery arm with a NaN gradient: F underflows as a
      # decision time goes to zero, and late in the distribution F(y2)
      # and F(y), each computed to 1e-13, differ by less than their own
      # error, so the difference came out NEGATIVE. So F is floored at
      # 1e-300, and at an interval's UPPER edge, recognized as frmtmb's
      # own cens_y2 vector, which is data, this returns F(y) PLUS the
      # interval's mass computed directly in log space
      # (ddm_rt_linterval_b()), with F(y) the same bits the lower-edge
      # call returns. frmtmb's difference is then that mass to the
      # rounding of F(y), never negative, and held at 4 ulp of F(y) or
      # above, a finite barrier where the mass is below what the scale
      # holds. The exact fix is a log-difference slot in frmtmb.
      Fq <- ddm_floor(exp(kb * lb + (1 - kb) * lF), 1e-300)
      y2 <- aterms[["cens_y2"]]
      ylo <- aterms[["cens_ylo"]]
      if (is.null(y2) || is.null(ylo) ||
          !identical(as.numeric(q), as.numeric(y2))) {
        return(Fq)
      }
      ki <- as.numeric(cen == 2)
      nd <- dpars[["ndt"]]
      l1 <- ddm_rt_lcdf_b(ylo - nd, dpars[["mu"]], dpars[["bs"]],
                          dpars[["bias"]], ddm_indicator(aterms))
      F1 <- ddm_floor(exp(l1), 1e-300)
      lD <- ddm_rt_linterval_b(ylo - nd, q - nd, dpars[["mu"]],
                               dpars[["bs"]], dpars[["bias"]],
                               ddm_indicator(aterms))
      Fi <- F1 + ddm_floor(exp(lD), 4 * .Machine$double.eps * F1)
      ki * Fi + (1 - ki) * Fq
    }
    lccdf_rows <- function(q, dpars, aterms) {
      ddm_rt_lcdf2(q - dpars[["ndt"]], dpars[["mu"]], dpars[["bs"]],
                   dpars[["bias"]])[["lS"]]
    }
    # frmtmb calls both slots on EVERY row and reads back only the
    # censored ones. Each is a few thousand tape operations per row, so
    # on 12,000 rows the tape did not fit in memory: every fit of the
    # left- and interval-censored recovery arm, and many of the
    # right-censored one, ended in std::bad_alloc. So each slot works on
    # the rows frmtmb will read, which the censoring code, as data,
    # names, and returns a neutral 1 (or log 1) on the others. Under
    # trunc() every row's distribution function is the normalizer, and
    # all rows are computed.
    #
    # On those rows every one is scored at its boundary, so the marginal
    # function is not formed at all, and at an interval's upper edge
    # only the interval rows are. The lower edge and the upper-edge call
    # form F(y) with the same arithmetic on the same values, so the two
    # agree to the bit and frmtmb's difference is the directly computed
    # mass (see lcdf_rows above for why the mass is computed directly).
    lcdf_edge <- function(q, dpars, aterms) {
      ddm_floor(exp(ddm_rt_lcdf_b(q - dpars[["ndt"]], dpars[["mu"]],
                                  dpars[["bs"]], dpars[["bias"]],
                                  ddm_indicator(aterms))), 1e-300)
    }
    lcdf_upper <- function(q, dpars, aterms) {
      nd <- dpars[["ndt"]]
      ylo <- aterms[["cens_ylo"]]
      up <- ddm_indicator(aterms)
      l1 <- ddm_rt_lcdf_b(ylo - nd, dpars[["mu"]], dpars[["bs"]],
                          dpars[["bias"]], up)
      F1 <- ddm_floor(exp(l1), 1e-300)
      lD <- ddm_rt_linterval_b(ylo - nd, q - nd, dpars[["mu"]],
                               dpars[["bs"]], dpars[["bias"]], up, l1)
      F1 + ddm_floor(exp(lD), 4 * .Machine$double.eps * F1)
    }
    lcdf <- function(q, dpars, aterms) {
      cen <- aterms[["cens"]]
      if (is.null(cen) || length(cen) != length(q) ||
          !is.null(aterms[["trunc_lb"]]) || !is.null(aterms[["trunc_ub"]])) {
        return(lcdf_rows(q, dpars, aterms))
      }
      y2 <- aterms[["cens_y2"]]
      if (!is.null(y2) && !is.null(aterms[["cens_ylo"]]) &&
          identical(as.numeric(q), as.numeric(y2))) {
        return(ddm_on_rows(lcdf_upper, q, dpars, aterms, which(cen == 2), 1))
      }
      ddm_on_rows(lcdf_edge, q, dpars, aterms, which(cen == -1 | cen == 2), 1)
    }
    lccdf <- function(q, dpars, aterms) {
      cen <- aterms[["cens"]]
      if (is.null(cen) || length(cen) != length(q)) {
        return(lccdf_rows(q, dpars, aterms))
      }
      ddm_on_rows(lccdf_rows, q, dpars, aterms, which(cen == 1), 0)
    }
  }

  mean_fn <- if (!length(vv)) {
    function(dpars, aterms) ddm_mean_rt(dpars, aterms)
  } else {
    # The closed-form conditional mean is the mean of the WRONG model
    # once the parameters vary between trials, and a post-fit method
    # that quietly returns a number for the wrong model is worse than
    # one that refuses.
    gh <- ddm_gauss_hermite(21L)
    function(dpars, aterms) ddm_mean_rt_var(dpars, aterms, nd, gh)
  }
  sim <- if (!length(vv)) {
    function(dpars, aterms, n) ddm_sim_rt(dpars, aterms, n)
  } else {
    function(dpars, aterms, n) ddm_sim_rt_var(dpars, aterms, n, nd)
  }

  if (isTRUE(cfg$contaminant)) {
    dpars <- c(dpars, "lambda")
    links$lambda <- "logit"
    # Small and inside the link: zero is the edge, and a start there
    # would make the first gradient step decide whether contamination
    # exists at all.
    init$lambda <- function(y, aterms) 0.05
    pb <- ddm_boundary_prob_fn(vv, nd)
    lpdf <- ddm_cont_lpdf(lpdf, crange)
    if (!is.null(lcdf)) {
      lcdf <- ddm_cont_lcdf(lcdf, crange)
      lccdf <- ddm_cont_lccdf(lccdf, crange)
    }
    mean_fn <- ddm_cont_mean(mean_fn, pb, crange)
    sim <- ddm_cont_sim(sim, pb, crange)
  }

  fam <- frmtmb::custom_family(
    "wiener",
    accepts_aterms = ddm_accepts[["wiener"]],
    dpars = dpars,
    links = links,
    lpdf = lpdf,
    lcdf = lcdf,
    lccdf = lccdf,
    # The boundary a trial ended at reaches the density as dec() or as
    # vint1, and either will do, so the requirement is declared as the
    # choice it is rather than checked by hand after the frame is built.
    required_aterms = list(c("dec", "vint1")),
    # Either will do, and only one of them: ddm_indicator() reads dec
    # and falls back to vint1, so a model supplying both used to fit
    # with the second column silently unread, even when the two
    # contradicted each other. The pair order is that precedence.
    exclusive_aterms = list(c("dec", "vint1")),
    family_finalize = function(fam, y, aterms) {
      ddm_finalize(fam, cfg, y, aterms)
    },
    # one closure per family object, so the validator knows whether
    # this model has a distribution function to censor with
    valid_y = function(y, aterms) {
      ddm_check_response(y, aterms)
      if (length(vv) && (!is.null(aterms[["cens"]]) ||
                         !is.null(aterms[["trunc_lb"]]) ||
                         !is.null(aterms[["trunc_ub"]]))) {
        ddm_stop_cens_var()
      }
      cen <- aterms[["cens"]]
      if (!is.null(cen) && any(cen == -1 | cen == 2) &&
          (!is.null(aterms[["trunc_lb"]]) ||
           !is.null(aterms[["trunc_ub"]]))) {
        ddm_stop_cens_trunc()
      }
      invisible(NULL)
    },
    init_dpars = init,
    type = "continuous",
    post = list(mean_fn = mean_fn),
    sim = sim,
    sim_refusal = NULL)
  if (isTRUE(cfg$contaminant)) fam[["contaminant_range"]] <- crange
  # An interval's lower edge, which lcdf needs at the upper-edge call to
  # form the interval's mass directly; see lcdf above. Only when a row
  # is interval-censored, so every other model's data is what it was.
  if (!length(vv)) {
    fam[["aterm_data"]] <- function(y, aterms) {
      cen <- aterms[["cens"]]
      if (is.null(cen) || !any(cen == 2)) return(list())
      list(cens_ylo = as.numeric(y))
    }
  }
  fam
}

#' Fill in everything the family could not know until it had the data.
#'
#' The non-decision-time bound and the margin by which an unreachable
#' row is held off the singularity are both properties of the response.
#' `family_finalize()` is the slot frmtmb provides for exactly this: it
#' runs once at frame assembly, after the response is validated and
#' before any link is used, and whatever it returns is the family the
#' rest of the fit sees. Before that slot existed this package wrote the
#' bound into an environment the link closures read at call time, which
#' worked only for as long as the undocumented call order held.
#'
#' @noRd
ddm_finalize <- function(fam, cfg, y, aterms) {
  lo <- min(y)
  sp <- ddm_ndt_spec(y, aterms, cfg$max_ndt, "wiener")
  if (!is.null(cfg$max_ndt) && sp$ub > lo && !cfg$allow_unreachable &&
      !isTRUE(cfg$contaminant)) {
    # min(y) itself is allowed, and is the default: the logit never
    # reaches 1 at a finite linear predictor, so ndt < the bound stays
    # strict. Anything above min(y) does admit parameter values with no
    # likelihood.
    frm_stop("wiener: max_ndt = ", format(sp$ub), " is above the smallest ",
             "response time (", format(lo), "). The density is zero at ",
             "and below the non-decision time, so a bound above min(rt) ",
             "admits parameter values with no likelihood. In a mixture ",
             "the other component covers those trials, and ",
             "allow_unreachable = TRUE says so.", call. = FALSE)
  }
  # A settled bound is KEPT rather than re-derived. The family is
  # rebuilt from `cfg` here, so the install's own guard cannot see the
  # carried bound and has to be handed it: without this a leave-one-out
  # refit that dropped the fastest trial would be a refit of a
  # DIFFERENT model, which is what the review measured at 0.334522
  # against 0.345458.
  keep <- ddm_ndt_keep(fam)
  # The contaminant's range is a property of the fitted data for the
  # same reason, and is kept the same way.
  crange <- fam[["contaminant_range"]]
  if (isTRUE(cfg$contaminant) && (is.null(crange) || anyNA(crange))) {
    crange <- cfg$contaminant_range %||% ddm_cont_default_range(y, aterms)
    if (!(crange[2L] > crange[1L])) {
      frm_stop("wiener(contaminant = TRUE): every response time is ",
               format(crange[1L]), ", so the observed range the uniform ",
               "contaminant spreads over has zero width.", call. = FALSE)
    }
  }
  # A response outside the window has contaminant density zero, and a
  # row the Wiener part cannot reach then has no likelihood at all; a
  # window narrower than the data is a mistake in the window.
  if (isTRUE(cfg$contaminant)) ddm_cont_check_range(y, aterms, crange)
  # Under the contaminant a row faster than the non-decision time HAS a
  # likelihood, the contaminant's, so a `max_ndt` above the fastest
  # response is a model rather than a mistake. The Wiener part then
  # needs the density that holds such a row off the singularity, which
  # is the one allow_unreachable selects. Only then: at the default
  # bound every row is reachable and the plain density is kept, which
  # is what makes the family at lambda = 0 the plain one.
  ub <- keep[["ub"]] %||% sp$ub
  cfg2 <- cfg
  if (isTRUE(cfg$contaminant) && !is.null(cfg$max_ndt) && ub > lo) {
    cfg2$allow_unreachable <- TRUE
  }
  ddm_ndt_install(ddm_family(cfg2, delta = 1e-9 * lo, crange = crange),
                  ub,
                  if (is.null(keep)) sp$floors else keep[["floors"]],
                  "wiener",
                  sizes = if (is.null(keep)) sp$sizes else keep[["sizes"]])
}

#' The decision indicator, under whichever spelling supplied it.
#'
#' `dec()` is the term this package contributes to frmtmb's addition-term
#' registry and the one every reference on the model uses; `vint()` is
#' the general-purpose route that was the only one available before that
#' registry existed, and it keeps working.
#'
#' @noRd
ddm_indicator <- function(aterms) {
  up <- aterms[["dec"]]
  if (is.null(up)) up <- aterms[["vint1"]]
  up
}

#' Response and decision-indicator validation.
#'
#' The missing-indicator refusal is NOT here. It is declared, as
#' `required_aterms = list(c("dec", "vint1"))`, which frmtmb reads as
#' "either of these two spellings will do" and enforces during frame
#' assembly, before the frame is built rather than after. What is left
#' here is the part a declaration cannot express: whether the values in
#' the term are the two a boundary indicator may take.
#'
#' @noRd
ddm_check_response <- function(y, aterms) {
  if (any(!is.finite(y)) || any(y <= 0)) {
    frm_stop("wiener: the response must be a strictly positive, finite ",
             "response time.", call. = FALSE)
  }
  ddm_check_units(y, "wiener")
  up <- ddm_indicator(aterms)
  if (any(!is.finite(up)) || any(up != 0 & up != 1)) {
    frm_stop("wiener: the decision indicator must be 0 (lower boundary) ",
             "or 1 (upper boundary). dec() coerces a factor or a character ",
             "vector for you, taking its second level as the upper ",
             "boundary; vint() does not, so recode it with ",
             "as.integer(decision == \"upper\").", call. = FALSE)
  }
  invisible(NULL)
}
