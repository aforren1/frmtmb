# eLife 56938: what the paper actually fits

Status: working notes. Written incrementally; sections marked TODO are still being pinned.

## Identity

Shinn M, Lam NH, Murray JD (2020). *A flexible framework for simulating and fitting
generalized drift-diffusion models*. eLife 9:e56938.

This is a **software/methods paper**, not a single scientific claim. It introduces PyDDM
(<https://github.com/murraylab/PyDDM>) and demonstrates it by fitting a family of
drift-diffusion models to two published perceptual decision datasets. So the maintainer's
question "can we fit the model(s)" is really "can frmtmb express the generalized DDM
(GDDM) class", with the paper's five named variants as the concrete test set.

## Scientific setting

Random-dot-motion discrimination. On each trial a subject views a moving dot field at some
motion coherence `C` and reports the perceived direction; both the **choice** and the
**reaction time** are recorded. The modeling claim is that the standard DDM (constant
drift, constant bounds) is too rigid, and that a generalized DDM with leaky integration,
time-collapsing bounds and a nonlinear coherence-to-drift mapping fits monkey data
substantially better.

## Data structure

Two datasets, both fit **per subject**, not hierarchically.

**Roitman & Shadlen (2002)** - 2 rhesus monkeys, random dot motion, saccadic report.
Six coherence levels. Distributed with PyDDM as `doc/downloads/roitman_rts.csv`
(columns: `monkey, rt, coh, correct, trgchoice`; 6149 data rows). Coherences present are
0, 0.032, 0.064, 0.128, 0.256, 0.512. The paper's own script keeps monkey 1 only and
trims `0.1 < rt < 1.65` s "for compatibility with" Ratcliff & McKoon-era reanalyses.
Original source: <https://shadlenlab.columbia.edu/resources/RoitmanDataCode.html>.

**Evans & Hawkins (2019)**, *Cognition* - human random dot motion, no-feedback-delay
condition only, pooled across subjects. Deposited at <https://osf.io/2vnam/>.

Per-trial observation is therefore a **(choice, RT) pair** with one covariate (coherence),
i.e. a defective bivariate density: one density over RT for each of the two absorbing
boundaries, together integrating to less than 1 when the process can fail to terminate
within `T_dur`.

## Response coding

The paper's default is **accuracy coding** (`correct` in {0,1}, upper bound = correct
response), which is why drift is `mu0 * C` with no sign flip. PyDDM also ships a
stimulus-coded variant (`roitman_shadlen_stimulus_coding.py`) where the upper bound is a
fixed direction and drift changes sign with the stimulus. Both codings matter for the
frmtmb mapping: accuracy coding is the one that maps cleanly onto a two-boundary family
with a single non-negative drift.

## The models

### M1. 3-parameter DDM (paper Eq. 14) - the baseline

    dx = mu dt + sigma dW,   sigma = 1
    mu(x, t, C) = mu0 * C            (constant in x and t)
    B(t)        = B0                 (constant)
    X0          = delta(x)           (unbiased start)

Free: `mu0`, `B0`, `t_nd`. This is exactly the textbook Wiener first-passage-time model
with `z = 0.5`, no across-trial variability.

### M2. 8-parameter DDM (paper Eq. 16) - drift free per coherence

Same as M1 but `mu(x, t, C) = mu_j`, one independent drift per coherence level j
(6 levels). Free: 6 x `mu_j`, `B0`, `t_nd`. Unlike M1 and M5 this one carries **no lapse
overlay**: `p*_i(t) = p_i(t - t_nd)`. The paper fits it under HDDM's conventions,
including a fixed outlier probability `p_outlier = 0.05`.

### M3. 11-parameter "full DDM" (paper Eq. 15) - across-trial variability

Same as M2 plus the three Ratcliff across-trial variability terms:

    mu(x,t,...) ~ Normal(mu_j, s_v^2)                     (drift varies trial to trial)
    X0(x,...)   = 1/(2 s_z + 1) for -s_z <= x <= s_z      (start point varies)
    p*_i(t)     = 1/(2 s_T + 1) * sum over k in [-s_T, s_T] of p_i(t - t_nd - k)

Free: 6 x `mu_j`, `B0`, `t_nd`, `s_v`, `s_z`, `s_T`; `p_outlier = 0.05` fixed. The
likelihood is now a triple integral of the Wiener density over the three nuisance
distributions. The paper says outright that this model "does not fall into the GDDM
framework because the across-trial variability in drift rate is a random variable", and
fits it with HDDM rather than PyDDM.

### M4. 18-parameter DDM (paper Eq. 17) - saturated per condition

All three DDM parameters free per coherence level: `mu_j`, `B_j`, `t_nd,j` for each of
the 6 coherences. The paper obtains it via EZ-diffusion moment matching rather than full
likelihood, but as a *model* it is just M1 with every parameter indexed by condition.

### M5. 6-parameter GDDM (paper Eq. 13) - the paper's headline model

    mu(x, t, C) = mu0 * (C / C_max)^alpha - l * x   (leaky/unstable integration)
    sigma       = 1
    X0          = delta(x)
    B(t)        = B0 * exp(-t / tau)                (exponentially collapsing bound)
    p*_i(t)     = 0.95 * p_i(t - t_nd) + 0.025      (uniform contaminant mixture)

Free: `mu0`, `alpha` (coherence nonlinearity), `l` (leak; negative = leaky, positive =
unstable), `B0`, `tau`, `t_nd`. `C_max` is fixed, not fitted: 0.512 for Roitman &
Shadlen, 0.4 for Evans & Hawkins.

**The published Eq. 13 has a typesetting error.** It prints the drift term as
`mu0 (C - C_max)^alpha` with a genuine minus sign, and that error is in the JATS source
so it propagates identically to the HTML, the PDF and the PMC copy. It must be a
division. The paper's *own Appendix 1* gives the class it used:

    class DriftCoherenceLeak(ddm.models.Drift):
        required_parameters = ["driftcoh", "leak", "power", "maxcoh"]
        def get_drift(self, x, conditions, **kwargs):
            return self.driftcoh * (conditions["coh"]/self.maxcoh)**self.power \
                   + self.leak * x

so the drift is `mu0 * (C / C_max)^alpha + leak * x`, i.e. the paper's `l` is `-leak`.
The subtraction reading is also nonsense on its face: it is non-positive for every
`C <= C_max`, gives exactly zero drift at maximum coherence, and is not real-valued for
non-integer `alpha`. The bioRxiv preprint has no `alpha` at all (a 5-parameter GDDM with
plain `mu0 * C`), so both the nonlinearity and its typo arrived at revision. Anything
reimplementing this model from the printed equation will silently fit the wrong thing.

### What the paper does not report

**There is no table of fitted parameter values, for either dataset.** Nor is there a
single numeric BIC, log-likelihood or MSE anywhere in the body, the figure legends, the
tables, the appendix or the back matter. The model comparison exists only as unlabeled
bar plots (Figure 2a-c for monkey N, supplements for monkey B and the human data), and
the JATS XML carries no `source-data` elements, so the plotted values were never
deposited. A prototype therefore **cannot** be checked against the paper's own numbers.
The nearest available target is PyDDM's tutorial fit, below, which is a *similar but not
identical* model. That constrains what "recovery" can mean here and is the reason the
prototype below leans on simulation-based recovery rather than replication.

The tutorial version of this model that ships with PyDDM
(`doc/downloads/roitman_shadlen.py`) drops `alpha` and uses a Poisson (exponential)
mixture instead of the uniform one, so it has 5 free parameters. Its **published fitted
values for monkey 1** are given in the script as a copy-pasteable pre-fit model:

    driftcoh = 10.49091,  leak = -0.482,  B = 1.811,  tau = 1.992,
    nondectime = 0.211,   pmixturecoef = 0.02 (fixed),  rate = 1 (fixed)

with `dx = 0.01, dt = 0.01, T_dur = 2`. Note PyDDM parameterizes the collapsing bound as
`B * exp(-tau * t)` (see `pyddm/models/bound.py::BoundCollapsingExponential`), so its
`tau` is a *rate*, the reciprocal of the paper's time constant. These are the numbers a
prototype should try to reproduce.

## The likelihood, pinned from the source

The prose is less useful here than `pyddm/models/loss.py`. `LossLikelihood.loss()` does,
per condition group `k`:

    loglikelihood += sum(log(sols[k].pdf("_top")[idx_top]))
    loglikelihood += sum(log(sols[k].pdf("_bottom")[idx_bottom]))

where `idx` are the RT values rounded to the nearest `dt` bin. Two consequences that
decide the frmtmb mapping:

1. **The likelihood factorizes over trials.** It is a plain sum of per-row log densities.
   There is no cross-trial coupling, no sequential dependence, no latent state carried
   between trials. This is *not* `frmtmb_structure()` territory.
2. **The expensive object is shared within a condition.** The defective density is
   obtained once per unique condition (here, per coherence) by numerically solving the
   Fokker-Planck equation on a `t x x` grid (Crank-Nicolson, tridiagonal solves), then
   *looked up* at each trial's RT bin. Six solves serve all 6149 trials.

Note the paper itself states **no** loss function. There are 19 numbered display
equations and not one of them is a likelihood, a BIC or an MSE; those exist only as prose
and as code. The definitions above are therefore taken from PyDDM v0.5.0, the release
contemporary with the August 2020 paper. In that release undecided trials *do* contribute
(`loglikelihood += log(prob_undecided) * n_undecided`); current master comments that line
out as "not a valid way to incorporate undecided trials". The point is moot for the
Roitman fits, where `T_dur = 2` and RTs are trimmed below 1.65 s, so nothing is undecided. Contaminant mass is instead handled by the mixture overlay.

The contaminant overlay, from `pyddm/models/overlay.py::OverlayUniformMixture`:

    p_top(t)    <- p_top(t)    * (1 - c) + 0.5 * c * norm / n_t
    p_bottom(t) <- p_bottom(t) * (1 - c) + 0.5 * c * norm / n_t

with `norm = sum(p_top) + sum(p_bottom)`, i.e. a mixture with a uniform-over-`[0, T_dur]`
lapse distribution split evenly across the two responses. `OverlayPoissonMixture` is the
same construction with an exponential rather than uniform contaminant.

Non-decision time is a pure **shift** of the first-passage density, `p*(t) = p(t - t_nd)`,
applied as a post-hoc overlay rather than inside the diffusion.

## Fitting procedure

- Full-distribution **maximum likelihood** on the joint (choice, RT) distribution. No
  priors, no MCMC, no hierarchy. Per subject.
- Default optimizer: **differential evolution** (global). Nelder-Mead, BFGS, basin
  hopping and hill climbing are also offered.
- Model comparison by BIC (`log(n) * k - 2 * loglik`), raw negative log likelihood, and
  mean squared error against the empirical RT histograms.
- Densities come from a numerical **Fokker-Planck** solve (Crank-Nicolson preferred,
  backward Euler fallback) except when PyDDM detects that the model is a constant-drift,
  constant-bound DDM, in which case it substitutes the closed-form analytic solution
  (`pyddm/analytic.py`).

That last detail is the crux of the whole feasibility question: **M1, M2 and M4 take the
analytic route; M3 and M5 do not.**

## Licensing note for anything fetched

PyDDM is MIT licensed. `roitman_rts.csv` redistributed in its docs is derived from the
Roitman & Shadlen (2002) public release. Copies used for probing live only in the session
scratchpad and are not committed to this repository.
