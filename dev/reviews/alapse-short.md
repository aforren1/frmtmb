# Short review of lane alapse, fix round 1, 2026-11-17

Scope: the claim of fix round 1 (`dev/alapse-findings.md`, "Fix round
1"): code 0 without convergence warnings only where the maximum is proved
global, by convexity (fixed quantile, constant sigma) or by the search
over the quantile (estimated constant quantile). Worktree
`frmtmb-wt-alapse`, library `wt-alapse-lib`, base `rellib-r22`. Scripts
`dev/reviews/alapse-short-*.R`, runner `dev/reviews/alapse-short-run.sh`,
logs `dev/reviews/alapse-review-out/alapse-short-*.txt` (gitignored).
Every number below is pasted from those logs.

The library is the worktree's source: every function of `R/ald-fixed.R`,
`R/fit.R`, `R/confint.R`, `R/se-check.R` and `R/ald-laplace.R` deparses
equal to the installed namespace (`alapse-short-sync.txt`, 0 differ in
each).

## Verdict: FIX

1. (HIGH) The quantile search and its restart ignore the optimizer's
   box. With a bound on the quantile's or sigma's intercept
   (`set_prior(lb =, ub =)`), the fit moves outside the bound and
   reports code 0, "the highest maximum", and no convergence warning.
   rellib-r22 stays inside the bound and warns. A silent wrong result
   that 0.83.0 did not give.
2. (MEDIUM) `frmtmb_ald_local` and code 1 on fits that the lane's own
   proof covers: a constant sigma on a link other than log. The fit
   is the global maximum (equal to the log-link fit's logLik to 8e-11)
   and still warns.
3. (LOW, file) The same false warning on a model of two responses that
   share no parameter; the warning's remedy ("Fix the quantile and drop
   the scale model ... whose maximum is then the highest") does not
   apply there.

Everything else that the brief asked for holds (sections 1 to 6).

## 1. The convexity proof

Derivation. The asymmetric Laplace log density (`R/families.R`
`fam_asym_laplace()`, and brms) is
`log p + log(1 - p) - log sigma - rho_p((y - mu) / sigma)`, with
`rho_p(u) = u (p - [u < 0])`. With case weights `w_i >= 0` and an offset
`o_i` in mu, the negative log-likelihood at a fixed p is
`W log sigma + sum w_i rho_p(y_i - o_i - x_i b) / sigma + const`.
`rho_p` is positively homogeneous, so with `t = 1 / sigma`, `a = b / sigma`
it is `-W log t + sum w_i rho_p(t (y_i - o_i) - x_i a)`: `-log t` is
convex, and each `rho_p` of an affine map is convex. The map
`(log sigma, b) -> (t, a)` is a bijection of the parameter space, so a
local minimum is global. The proof holds:

- with case weights (non-negative multipliers of convex terms);
- with an offset in mu (it moves into the data, `y - o`);
- with the `zi` part of `zero_inflated_asym_laplace()`: brms and frmtmb
  read only `zi` on a zero row and `log(1 - zi) + ALD` on the others, so
  the objective separates into the ALD part on the nonzero rows and a
  Bernoulli part in the zi coefficients. Both of frmtmb's zi links
  (`dpar_links_unit`: logit, identity) give a concave Bernoulli
  log-likelihood;
- with a finite box on log sigma (an interval in t is convex), which the
  code also guards by `inside(p0)`.

It does not hold with censoring or truncation (a tail probability term),
a bound or fixed value on a location coefficient, a prior density (the
log prior is not known concave in `(t, a)`), or a sigma or quantile model.

`ald_global_case()` against that: censoring, truncation, a bound or a
fixed location coefficient, me(), autocor and a non-identity mu link
stop the proof earlier (`ald_lp_blocks()`), so no code 0 comes from it;
a prior density on any parameter, more than one response, and a sigma or
quantile model give "local". zi parameters are read (`zpos`); any other
parameter (a mo() simplex, for example) gives "local" by the `known`
test. So nothing the proof does not cover is called "convex". It
excludes more than it must (finding 2 and 3).

Measured (`alapse-short-fits.txt`): fixed quantile 0.3 with no extra
term, with weights, with an offset (rq objective relative difference 0),
zi with a fixed and an estimated quantile, zi with
`link_zi = "identity"`: each code 0, global, one warning
(`frmtmb_ald_no_se`).

## 2. The search over the quantile

`S(p) = min_b sum w_i rho_p(r_i)`: for fixed b each term is linear in p,
so S is a minimum of linear functions, concave and piecewise linear.
Sigma profiles in closed form (`sigma = S / W`), giving
`W (log p + log(1 - p) - log(S(p) / W) - 1)`, a decreasing function of
S plus a term in p alone. So a lower bound on S gives an upper bound on
the profile. On `[a, b]` a concave S lies above its chord; on the end
intervals the chord through 0 is lower still (`S >= 0`), so the bound
stays valid. I checked the closed forms: the end-interval bounds
`W (log(b / S_b) + log W - 1)` and `W (log((1 - a) / S_a) + log W - 1)`
are the suprema at p -> 0 and p -> 1 of the chord bound; the interior
stationary condition of `log p + log(1 - p) - log(c0 + c1 p)` is
`c1 p^2 + 2 c0 p - c0 = 0`, as coded, with the ends as candidates. Each
`S(p)` comes from a vertex that passes the subgradient test, so it is
exact. The bound object is right, and the proof is valid where the box
does not bind (finding 1).

Brute force (`alapse-short-search.R`, `alapse-short-search.txt`): design
B (100 rows, y = 1 + 0.5 x + N(0, 1), seed 2000 + s), seeds 109, 22,
140, 2, 3, 4, 12, 20, 101, 102 and fresh 601 to 610. Four variants:
plain, weights U(0.2, 3), ties (y rounded to integers, x to 0.1), and
both. The global maximum comes from `rq.wfit(method = "br")` on a 1999-point
grid of p (step 5e-4) and, exactly, from the weighted quantile process
maximized on each interval:

    SUMMARY plain n 20 | code0 20 | global 20 | ald_local 0 | conv-warn 0 | gap_grid > 1e-6: 0 (max -9.91e-08) | gap_exact > 1e-6: 0 (max 2.84e-14) | global & gap_exact > 1e-6: 0
    SUMMARY w     n 20 | code0 20 | global 20 | ald_local 0 | conv-warn 0 | gap_grid > 1e-6: 0 (max -1.7e-09) | gap_exact > 1e-6: 0 (max 8.44e-11) | global & gap_exact > 1e-6: 0
    SUMMARY ties  n 20 | code0 20 | global 20 | ald_local 0 | conv-warn 0 | gap_grid > 1e-6: 0 (max 2.42e-12) | gap_exact > 1e-6: 0 (max 1.55e-11) | global & gap_exact > 1e-6: 0
    SUMMARY wties n 20 | code0 20 | global 20 | ald_local 0 | conv-warn 0 | gap_grid > 1e-6: 0 (max -4.4e-08) | gap_exact > 1e-6: 0 (max 1.52e-11) | global & gap_exact > 1e-6: 0
    plain 109 code 0 global TRUE ll -134.505740 exact -134.505740 p 0.7247 p_exact 0.7247 evals 34 moved 0.4197

80 of 80 fits are at the global maximum (largest shortfall 8.44e-11),
and the search took 31 to 37 linear programs.

## 3. Tied data

`alapse-short-ties.R` (`alapse-short-ties.txt`).

- Part 1: 10 integer-y (rpois, x in 0:3) and 10 rounded-y (one
  decimal, x in 0:4) data sets, n = 60, quantiles 0.25, 0.5, 0.75 and
  0.9, unit and U(0.5, 2) weights: 160 problems. Plain `ald_vertex()`
  returns NULL on 102 of them. `ald_vertex_try()` and
  `ald_vertex_ties()` certify 160 of 160, objective against rq "br"
  within 1.98e-16 relative.
- Part 2, adversarial: `ald_vertex_ties()` evaluated with
  `ald_tie_eps` raised (a copy of the function in an environment that
  redefines the constant; no package code changed), so the
  perturbed optimum's basis is often not optimal on the data:

      PART2 eps 0.001  | perturbed basis not optimal on the data: 0, certified 0 | optimal: 160, certified 160 | ...
      PART2 eps 0.05   | perturbed basis not optimal on the data: 53, certified 0 | optimal: 107, certified 97 | worst certified rel 1.69e-16
      PART2 eps 0.3    | perturbed basis not optimal on the data: 155, certified 0 | optimal: 5, certified 1 | worst certified rel 0
      PART2 eps 1      | perturbed basis not optimal on the data: 157, certified 0 | optimal: 3, certified 0 | worst certified rel 0

  Of 365 non-optimal bases, the proof certified none. Every certified b
  is at rq's objective. The proof reads only the data as they are, as
  the code comment says.
- Part 3, ties below the perturbation (integer y plus k 1e-11): 37 of
  60 certified, worst relative objective 2e-16. The other 23 keep the
  optimizer's warnings, which is safe.
- Part 4: frm() with a fixed quantile (0.25, 0.5, 0.9) on the part-1
  data: 60 of 60 code 0, global, `frmtmb_ald_no_se` only, objective
  against rq at most 0 relative (min -1.5e-16).

## 4. frmtmb_ald_local

`alapse-short-fits.txt`:

    sigma ~ z estq                     code 1 | cert ok TRUE  global FALSE | warn frmtmb_ald_local,frmtmb_ald_no_se
    sigma ~ z fixq                     code 1 | ... | warn frmtmb_ald_local,frmtmb_ald_no_se
    quantile ~ z                       code 1 | ... | warn frmtmb_ald_local,frmtmb_ald_no_se
    multi-response (mvbind)            code 1 | ... | warn frmtmb_ald_local,frmtmb_ald_no_se
    multi-response (bf + bf)           code 1 | ... | warn frmtmb_ald_local,frmtmb_ald_no_se
    prior density on sigma, fixq       code 1 | ... | warn frmtmb_ald_local,frmtmb_ald_no_se
    prior density on quantile, estq    code 1 | ... | warn frmtmb_ald_local,frmtmb_ald_no_se
    prior density on b, fixq           code 1 | cert ok FALSE | warn frmtmb_conv_code,frmtmb_conv_gradient,frmtmb_ald_no_se

A prior density on b stops the proof itself ("a prior density on the
coefficients"), so the optimizer's own warnings stay. That is correct.
On the proved-global fits, `frmtmb_ald_local` fired 0 times in 80 (the
search fits), 60 (part 4) and 8 (section 1).

check_se: `frmtmb_ald_local` is a convergence warning, emitted by
`check_convergence()` and not under `check_se`. With
`check_se = "ignore"`, the warnings are `frmtmb_ald_local` only. With
`"stop"`, `frmtmb_ald_local` is emitted, then the error
`frmtmb_ald_no_se`. On a global fit with `"ignore"`: no warning.
This is consistent with the generic convergence classes.

## 5. confint

On the B109 fit, `method = "profile"` and `"uniroot"` stop with class
`frmtmb_ald_no_se` for `x`, `(Intercept)`, `sigma_(Intercept)` and
`quantile_(Intercept)`. On the zi fit (seed 9), `zi_x` profile 0.51849
to 1.41492 and uniroot 0.51959 to 1.41514, against glm()'s profile
0.51832 to 1.41516. `zi_(Intercept)` is the same. `parm = c("zi_x",
"x")` stops for the whole call (it names `x`), which is acceptable.

## 6. Tests on the lane library

`alapse-short-tf-*.txt`, one file per process:

    RESULT test-ald-fixed.R pass=161 fail=0 err=0 skip=0 warn=0
    RESULT test-se-limit.R pass=380 fail=0 err=0 skip=0 warn=0

## Finding 1 (HIGH): the search leaves the box

`alapse-short-bounds.R lane|base` (`alapse-short-bounds-lane.txt`,
`-base.txt`), design B seed 109 (global p 0.7247, logit 0.968, log
sigma -1.267; local p 0.4197, logit -0.324, log sigma -1.041):

    lane: quantile Int in [-2, 0.3]  code 0 global TRUE  | p 0.7247 sigma 0.2817 logLik -134.505741 | par 1.5938 0.3995 -1.2669 0.9681 | box [-Inf -Inf -Inf -2] .. [Inf Inf Inf 0.3] | warn frmtmb_ald_no_se
    lane: sigma Int in [-1.1, 2]     code 0 global TRUE  | ... | par 1.5938 0.3995 -1.2669 0.9681 | box [-Inf -Inf -1.1 -Inf] .. [Inf Inf 2 Inf] | warn frmtmb_ald_no_se
    base: quantile Int in [-2, 0.3]  code 1 global FALSE | p 0.4195 ... | par 0.7258 0.4053 -1.0417 -0.3249 | warn frmtmb_conv_code,frmtmb_conv_gradient
    base: sigma Int in [-1.1, 2]     code 1 global FALSE | p 0.4195 ... | par 0.7258 0.4053 -1.0417 -0.3250 | warn frmtmb_conv_code,frmtmb_conv_gradient

The formulas are `bf(y ~ x, quantile ~ 1)` and `bf(y ~ x, sigma ~ 1)`
with `set_prior("", class = "Intercept", dpar = ..., lb =, ub =)`, and
`restarts = 0` gives the same result. Cause: `ald_quantile_search()`
searches all of (0, 1) and profiles sigma without a bound. In
`ald_polish()`, the restart sets `start` from the search, and
`ald_local()` checks `inside()` only on its line-search candidates, so a
start outside the box whose Newton gain is about 0 settles where it
starts. Neither `proved(run2)` nor the final `f1 <= f0` test reads the
box. A `uniform()` bound on the same coefficient (with
`frmtmb_prior_unbounded`) did not move. The fit kept p 0.4197 with
`frmtmb_ald_local` ("the run from the highest quantile found did not
settle"), which is safe.

Required: in `ald_global_case()`, give "search" only when the sigma and
quantile coefficients have no finite bound (or restrict the search to
the box's p range and keep "local" when a sigma bound exists); in
`ald_polish()`, accept `run2` only when `inside(run2$par)`. Add the two
cases above as tests. They fail on this build.

## Finding 2 (MEDIUM): a false local warning on a non-log sigma link

`alapse-short-over.txt`:

    log link     : logLik -304.39157849 code 0 warn frmtmb_ald_no_se
    softplus link: logLik -304.39157849 code 1 warn frmtmb_ald_local,frmtmb_ald_no_se
    diff -8.14e-11

`ald_global_case()` requires `link == "log"` for sigma. For one constant
sigma, any link is a bijection of its coefficient onto sigma > 0, so
the convexity proof (and the search, which works in sigma and p
directly) holds whatever the link. The fit is at the global maximum,
and it reports code 1 and a warning that says it may not be. By the
round's rule ("warns on a correct fit is worse than none"), drop the
link condition for sigma (and for the quantile in the search, which is
over p), with the box condition of finding 1.

## Finding 3 (LOW, file)

Two asym_laplace responses with a fixed quantile each, sharing no
parameter: joint logLik -577.34707996 is the sum of the single fits'
logLik (-301.51649181 + -275.83058815, difference -1.64e-11), yet the
fit has code 1 and `frmtmb_ald_local`. Each response separates, so each
can be proved alone. The warning's remedy ("Fix the quantile and drop
the scale model ... whose maximum is then the highest") is wrong for
this case, and for a prior density. File it, or word the remedy by the
`local` reason.
