# More than two choices in gddm(): feasibility probe

Date: 2026-09-07. Lane `wt-gddm`, worktree
`C:/Users/adf44/source/r/frmtmb-wt-gddm`, branch `wt-gddm`, base c18253e
(frmtmb 0.53.0, frmtmb.eam 0.4.0). RTMB 1.9, numDeriv 2016.8.1.1,
R 4.6.1, Windows 11. Probe scripts in `dev/gddm-nchoice/`; they run from
the worktree root as `Rscript dev/gddm-nchoice/probeN-*.R`.

**No package code was changed.** Nothing under any `R/`, nothing in any
`DESCRIPTION`. Every prototype is a standalone script; four of them
reach into `frmtmb.eam:::` for the shipped 1-D solver, the `gd_b3`
B-spline and `ddm_floor`, which is the idiom
`tests/testthat/test-gddm-solver.R` already uses.

## Verdict

**Do not extend `gddm()` to more than two choices.** The question was
asked of `gddm()`, and `gddm()` exists for state-dependent drift and
moving boundaries. The only route that keeps those is a grid solve on an
(n-1)-dimensional domain. Measured at three alternatives, and at the
grid it needs to reach even **8 to 17 times worse accuracy** than the
shipped one-dimensional family, that route costs **106.6 s of tape build
against 8.56 s and 483 ms per gradient against 6.7 ms**: 12 times the
build and 72 times the gradient, for an answer an order of magnitude
worse. At four alternatives it needs about **20 minutes of tape build
and tens of seconds per gradient, per condition**. That is not a family
anyone can fit.

Three findings qualify that, and the third is the one worth arguing
about.

1. **Route 1, the race, does not cover the gap, and this is provable
   rather than arguable.** A genuine multi-alternative diffusion is a
   RELATIVE-evidence model: adding the same constant to every drift
   changes nothing at all. Measured over 60000 trials, its choice
   proportions and median response time are identical to every printed
   digit at drift +0, +1 and +3. `rdm()` moves its median from 0.609 to
   0.396 over the same shift. The two families are not
   reparameterizations of each other.

   Fitted rather than reasoned about: over 20000 trials from the
   three-alternative diffusion, `rdm(3)` and the diffusion have the same
   six parameters and `rdm(3)` lands **91.0 log-likelihood units** below
   it. That is 30 times the roughly p/2 = 3 a correctly specified model
   gets for free, so the misspecification is not in doubt. But it is
   **0.0046 nats per trial**, so a realistic 500-trial session would
   expect a gap of about 2.3, which is smaller than that same free 3.
   The race is provably the wrong model and you need tens of thousands
   of trials from one subject to notice. That is the whole case for
   doing nothing, and it is section 6.

2. **Route 2, the (n-1)-dimensional Fokker-Planck solve, works, is
   differentiable, and is not worth its price.** It tapes, its gradient
   matches numDeriv to 9e-8 at worst over five parameters, it conserves
   mass to 1.5e-4, and it needed
   one non-obvious correction that a plausible implementation gets wrong
   and does not notice (section 4.2). Its cost is set by the node count,
   which is `N^(n-1)`, and the measured constant is 33 to 53 microseconds
   of tape build per node-step, **flat in the dimension**.

3. **Route 3 exists, is EXACT, and reaches exactly three alternatives
   and nothing else.** For three alternatives the continuation region is
   an equilateral triangle, which is a reflection-group fundamental
   domain, so the method of images gives a closed-form first-passage
   density: a finite alternating sum of Gaussians times a difference of
   two normal distribution functions, with a Girsanov tilt for the
   drift. It tapes to SECOND order (Hessian correct to 7.6e-7), builds a
   tape 4.3 times faster than the shipped `gddm()`, and has no grid
   error of any kind. At n = 2 it reduces exactly to `wiener()`, so it
   is the honest n-choice generalization of the analytic family in a way
   the race is not.

   It is still probably not worth shipping, and section 6 says why: it
   is three alternatives only, one stopping rule only, equal thresholds
   only, constant drift only and fixed boundaries only. It is therefore
   NOT a `gddm()` at all. It could take none of `gddm()`'s drift,
   boundary or starting-point components. It would be a new family next
   to `wiener()`, and the stopping rule it is stuck with is the one that
   happens to be solvable rather than the one the literature writes down
   most often.

**Recommendation.** Ship nothing. Keep the existing refusal in
`gd_check_response()`, which already points a third response level at
`lba()`. One word of it is wrong and should be corrected when some other
lane is editing that file anyway: it says multi-alternative choice
"needs racing accumulators rather than one accumulator between two
walls". Section 3 shows that is false. It needs a higher-dimensional
absorbing domain, and a race is one particular choice of domain.

If the three-alternative closed form is ever wanted, section 5 is the
recipe and `dev/gddm-nchoice/images.R` is working code. It is about two
days as a standalone family, and it is not `gddm()`.

## 1. What the shipped family costs, for scale

`probe0-baseline.R`. Every number later is a multiple of these.

The shipped `gddm()` at its default grid (`dt = 0.01`, `ny = 201`), over
`t_max = 2` and 2000 trials:

| conditions | tape build | one `fn` | one gradient |
| --- | --- | --- | --- |
| 1 | 8.56 s | under 5 ms | 6.7 ms |
| 2 | 17.53 s | 3.3 ms | 10.0 ms |
| 4 | 37.41 s | 6.7 ms | 13.3 ms |
| 8 | 109.21 s | 3.3 ms | 53.3 ms |

Accuracy, against the closed-form Wiener density it reduces to at
constant drift and fixed boundaries: worst absolute error in the log
density past 0.2 s is **0.00176**, which reproduces the 0.0018 in
`gddm-plots-findings.md`. Unrenormalized mass over [0, 2] is 0.996831.

One solve alone, in plain doubles, is 0.71 s at `dt = 0.01`, `ny = 201`
(`probe5-cost.R`).

Keep one thing in view for every comparison below. The shipped solver
runs its Thomas sweep as a scalar R loop of length `ny` inside every
step, and builds its tape at **213 microseconds per node-step**
(8.56 s over 201 by 200). The prototypes of sections 4 and 5 vectorize
the same recurrence across grid lines and build at **33 to 53
microseconds per node-step** (section 4.4). The prototypes are
therefore about five times better implemented per node than the family
they are being measured against, so every cost below is flattering to
route 2, not the reverse.

## 2. The geometry, and why three is special

`dev/gddm-nchoice/common.R`.

Take n accumulators, `dx_i = mu_i dt + dW_i`, unit independent noise.
Only relative evidence decides, so project onto the sum-zero plane. With
`U` an orthonormal basis of that plane (n by n-1), `z = U'x` has
`Cov(dz) = U'U = I`: **the reduced process is isotropic Brownian motion
in n-1 dimensions**, with drift `U'mu`. That is Roxin's reduction
(*J. Math. Neurosci.* 9:5, 2019, which gives the reduction and states
that first-passage statistics there need simulation or a numerical
method, with no closed form). The isotropy is the whole reason anything
below works.

The stopping rule: alternative k wins when its evidence beats the mean
of all n by c, that is `<v_k, z> >= c` with `v_k = U'(e_k - 1/n)`. The n
vectors `v_k` have equal length and equal pairwise angles, so the
continuation region is a **regular (n-1)-simplex** centred on the
origin.

At n = 2 the rule is `x_1 - x_2 = +/- 2c`, so `z = (x_1-x_2)/sqrt(2)` is
a unit-diffusion Wiener process with drift `(mu_1-mu_2)/sqrt(2)`
absorbed at `+/- c sqrt(2)`. In `wiener()`'s parameterization that is
`bs = 2 sqrt(2) c`, `v = (mu_1-mu_2)/sqrt(2)`, `bias = 0.5`. The
one-dimensional case is `wiener()` exactly, which is the check the whole
construction has to pass.

At n = 3 the simplex is an equilateral triangle. **The group generated
by reflections in the three sides of an equilateral triangle tiles the
plane** (it is the affine Weyl group of A2 and the triangle is its
alcove). That, and only that, is what makes section 5 possible.

At n = 4 the simplex is a regular tetrahedron, which does **not** tile
three-space by reflections, and whose Dirichlet spectrum is not known in
closed form either. The alcove of affine A3 is a tetrahedron but not a
regular one. So the closed form stops at three for a structural reason,
not for want of effort.

Geometry checks printed by `probe1-images.R`: the three unit normals
have unit length, pairwise inner products of exactly -0.5 (120 degrees),
and sum to zero.

## 3. Route 1: what the race cannot express

`probe4-race.R`, measurement 1.

`lba(n)` and `rdm(n)` race n accumulators to ABSOLUTE thresholds and
factor the likelihood as one accumulator's density times the others'
survivals. The multi-alternative diffusion above stops on RELATIVE
evidence. The two differ in a way that is directly observable and needs
no fitting to see.

Add the same constant to every drift. The relative-evidence model cannot
notice: `mu_k - mean(mu)` is unchanged, so the joint distribution of
(choice, time) is identical. 60000 trials per row:

| model | choice proportions | median rt |
| --- | --- | --- |
| triangle, drift + 0 | 0.6193 0.2999 0.0808 | 0.5930 |
| triangle, drift + 1 | 0.6193 0.2999 0.0808 | 0.5930 |
| triangle, drift + 3 | 0.6193 0.2999 0.0808 | 0.5930 |
| `rdm()`, drift + 0 | 0.4268 0.3305 0.2427 | 0.6094 |
| `rdm()`, drift + 1 | 0.4058 0.3293 0.2649 | 0.4888 |
| `rdm()`, drift + 3 | 0.3842 0.3324 0.2834 | 0.3958 |

The triangle rows are identical to every digit printed, which they must
be. The point is that `rdm()`'s are not: a shift of +3 takes 35 percent
off its median response time. A stimulus manipulation that raises the
overall evidence for every alternative while holding the contrasts fixed
is therefore a direct test between the two families, and neither can
imitate the other on it.

So the race is a DIFFERENT model, not a cheaper parameterization of the
same one. The refusal in `gd_check_response()` is wrong where it says a
multi-alternative diffusion "needs racing accumulators". It needs a
higher-dimensional absorbing domain. The race is the special case where
the domain is a corner box and the noise is independent, which is
exactly why it factorizes into a density times survivals.

### 3.1 How big the misfit is when you fit it

`probe4-race.R`, measurement 2. 20000 trials from the three-alternative
diffusion with `mu = (0.9, 0.1, -1.0)`, `c = 0.7071`, `ndt = 0.25`.
Every model here has exactly six free parameters, so the log-likelihoods
compare directly and AIC adds nothing.

| model | logLik | fit time | gap |
| --- | --- | --- | --- |
| triangle at the TRUTH | -14778.307 | - | 1.62 |
| triangle at its MLE | -14776.687 | 289 s | 0 |
| `rdm(3)` at its MLE | -14867.677 | 44 s | **91.0** |
| `lba(3)` at its MLE | -15136.934 | 18 s | **360.2** |

The triangle model's own MLE beats the truth by 1.62, against the
roughly p/2 = 3 a correctly specified six-parameter model gets for free,
which is the check that the closed form and the simulator agree about
what model they are.

`rdm(3)` is 91.0 below. That is 30 times the free 3, so the
misspecification is not a sampling accident. `lba(3)` is 360.2 below,
but its fit reported false convergence and a maximum absolute gradient
of 4.56, so read 360.2 as a bound rather than as its best.

`rdm()` gets as close as it does by degenerating: its fitted `v3` is
-18.64 on the log link, so the third accumulator's drift rate is about
8e-9 and it never finishes on its own, and its fitted `A` is 0.026, so
the start-point range collapses. It reproduces the data by turning
itself into something that is not a race.

One caveat that runs the safe way. The data are drawn from a 250-tile
density and scored by a 118-tile one, so the roughly 0.5 percent of
trials past the 118-tile horizon (section 5.3 puts it at 2.52 s of
decision time at this r) are scored by a truncated series. That can only
push the triangle model's log-likelihood DOWN, so 91.0 understates the
gap rather than inflating it. The 1.62 between the truth and the MLE,
which uses the same 118 tiles for both, says the effect is small.

**And 91.0 over 20000 trials is 0.0046 nats per trial.** A 500-trial
session, which is a normal amount of data from one subject in a
three-alternative task, would expect a gap of about 2.3, smaller than
the 3 a correctly specified model collects for free. So the misfit is a
large-sample fact. It is real, it is not detectable in one subject's
data, and section 6 treats that as the answer.

## 4. Route 2: the (n-1)-dimensional Fokker-Planck solve

`dev/gddm-nchoice/pde2d.R`, `probe3-pde.R`, `probe5-cost.R`.

### 4.1 The scheme

Written for n = 3. Coordinates `s_k = <v_k, z>`, keeping `s1` and `s2`
(`s3 = -s1-s2`): the continuation region is then
`{s1 < c, s2 < c, s1+s2 > -c}`, a right triangle two of whose sides are
grid lines and whose third is a grid diagonal, so no side is
staircased. The price is a noise covariance with a cross term,
`Sigma = [[2/3,-1/3],[-1/3,2/3]]`, because the `s_k` are not orthogonal.
Whitening to make the noise isotropic turns the triangle equilateral and
puts two of its sides off the grid, which is worse.

Douglas ADI: one explicit stage carrying the whole operator including
the cross term, then two implicit corrections, one per direction, each a
set of independent tridiagonal systems.

The lines are ragged, because the domain is half a square. Rather than
sweep ragged lines, the solve runs over the FULL square and the rows
outside the triangle carry the identity, `(lo, di, up) = (0, 1, 0)`, with
a zero right-hand side. The Thomas sweep then leaves them at zero and
they decouple exactly, so every line has the same length and the whole
recurrence **vectorizes across lines**. That is what keeps the tape
build from being one R-level scalar iteration per node per step. It
solves about twice as many nodes as the triangle has, which is counted
in the cost.

The start is spread with the same cubic B-spline `gddm()` uses for its
non-decision-time shift, for the same reason: `c` is a parameter, so the
fractional grid index of the start point moves with a parameter, and
only the WEIGHTS may move, never an index.

### 4.2 The correction a plausible implementation gets wrong

The cross-derivative stencil at a node one layer inside the DIAGONAL
wall reads `(i-1, j-1)`, which is a layer OUTSIDE the domain. Masking it
to zero is the obvious thing and it is wrong: Dirichlet pins p on the
wall, not beyond it, and the smooth continuation across a wall whose
reflection preserves `Sigma` is ODD. This wall's reflection does
preserve `Sigma`, checked algebraically. Reflecting `(i-1, j-1)` across
`i + j = N` gives `(N-j+1, N-i+1)`, which for a layer-one node is the
node itself, so the ghost value is exactly `-p` at the node.

The two axis-aligned walls have no ghost at all, because every interior
node has i and j at least 2 and at most N-1, so their diagonal stencil
neighbours are either interior or exactly on a wall where p is genuinely
zero. Only the oblique wall has the problem, and only it is wrong.

Measured at `t_max = 2` and 400 steps, against the closed form of
section 5. Errors are the worst absolute error in the log density over
t in [0.15, 2].

| | mass | wall 1 | wall 2 | wall 3 |
| --- | --- | --- | --- | --- |
| N = 81, ghost masked to zero | 0.97938 | 0.0085 | 0.0063 | **0.2856** |
| N = 81, odd ghost | **1.000153** | 0.0139 | 0.0160 | 0.0305 |

Masked, the diagonal wall's flux is about 24 percent low and **does not
converge under refinement**: the error runs 0.300, 0.313, 0.334, 0.375
as N goes 24, 36, 54, 81, while the two grid-aligned walls converge
fine. Total mass sits at 0.979 instead of 1.

A three-fold symmetry check makes the same point without needing any
reference: at zero drift from a centred start the three fluxes must be
equal, and masked they come out 0.3317, 0.3317, 0.2554. Corrected they
come out 0.3317, 0.3317, 0.3315 at N = 36.

This is the kind of thing that decides whether route 2 is a week or a
month of work. Two of the three walls converge whatever you do; the
third is silently wrong until someone works out the reflection.

### 4.3 Accuracy and cost, corrected

| N | nt | in-domain nodes | one solve | mass | wall 1 | wall 2 | wall 3 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 24 | 400 | 253 | 0.17 s | 0.996590 | 0.0860 | 0.0984 | 0.1802 |
| 36 | 400 | 595 | 0.28 s | 0.999373 | 0.0439 | 0.0499 | 0.0974 |
| 54 | 400 | 1378 | 0.42 s | 1.000126 | 0.0235 | 0.0267 | 0.0540 |
| 81 | 400 | 3160 | 0.83 s | 1.000153 | 0.0139 | 0.0160 | 0.0305 |

Convergence is between first and second order in N, not the clean second
order the one-dimensional Crank-Nicolson gets, because the cross term is
explicit. Against the shipped family's 0.00176: at 3160 nodes and twice
the shipped time resolution, **route 2 is 8 to 9 times less accurate on
the grid-aligned walls and 17 times less accurate on the oblique one**.

Tape, 2000 trials, one condition. The last two rows are the grids whose
accuracy is quoted in the table above, so they are the ones the
comparison has to be made on.

| N | nt | in-domain nodes | tape build | one gradient |
| --- | --- | --- | --- | --- |
| 24 | 100 | 253 | 1.27 s | under 5 ms |
| 36 | 100 | 595 | 2.23 s | 6.7 ms |
| 36 | 200 | 595 | 6.76 s | 26.7 ms |
| 54 | 200 | 1378 | 10.96 s | 46.7 ms |
| **54** | **400** | **1378** | **42.44 s** | **213 ms** |
| **81** | **400** | **3160** | **106.64 s** | **483 ms** |

Against the shipped family's 8.56 s and 6.7 ms: to reach an answer that
is still 8 to 17 times worse, route 2 wants **12 times the tape build
and 72 times the gradient**. And section 1 applies: the shipped family
would itself get about five times faster if its Thomas sweep were
vectorized the way this prototype's is, so the honest ratio is worse
still.

The gradient is right. Against numDeriv at N = 24, nt = 100, 200 trials,
the five parameters agree to 2.1e-11, 1.3e-9, 2.4e-9, 4.3e-8 and 9.1e-8
relative.

### 4.4 What four alternatives would cost

`probe5-cost.R` runs the SAME Douglas ADI in one, two and three
dimensions on a cube, at matched resolution per axis, so the exponent in
the dimension is measured rather than assumed. The cube is not the
regular tetrahedron four alternatives actually needs; the cost of the
scheme is set by the node count and the number of directional sweeps,
both of which the cube reproduces. (The cube is also not a fiction: it
is the domain of the n-accumulator diffusion with ABSOLUTE thresholds
and correlated noise, which is the race with its independence assumption
removed.)

| d | nodes/axis | nt | nodes | tape build | one gradient | build per node-step |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 51 | 100 | 51 | 6.28 s | 3.3 ms | 1231 us |
| 1 | 51 | 200 | 51 | 8.33 s | under 5 ms | 817 us |
| 2 | 25 | 100 | 625 | 3.31 s | 16.7 ms | 53.0 us |
| 2 | 37 | 100 | 1369 | 5.89 s | 30.0 ms | 43.0 us |
| 2 | 37 | 200 | 1369 | 11.16 s | 100 ms | 40.8 us |
| 3 | 13 | 100 | 2197 | 7.50 s | 160 ms | 34.1 us |
| 3 | 17 | 100 | 4913 | 16.36 s | 166.7 ms | 33.3 us |
| 3 | 21 | 100 | 9261 | 33.53 s | 1.48 s | 36.2 us |

The build cost per node-step is **flat in the dimension** at 33 to 53
microseconds once the grid is big enough for per-call overhead to stop
mattering. The d = 1 rows are the overhead-dominated end, and they are
why the shipped family's 8.56 s buys so little: at 51 nodes almost all
of it is R-level call overhead rather than tape.

So the cost of this route is entirely the node count, and the node count
is `N^(n-1)`.

Extrapolating at the measured rates to four alternatives at the
resolution section 4.3 needs (53 nodes per axis, 200 steps):
`53^3 = 148877` padded nodes, `2.98e7` node-steps, so **about 20 minutes
of tape build** at 40 us per node-step and, from the d = 3, n = 21 row's
1.6 us per node-step, **tens of seconds per gradient**, per condition. A
fit with 50 gradient evaluations over four conditions is hours of
gradients on top of hours of tape building. Five alternatives is a
four-dimensional grid and needs no arithmetic to rule out.

## 5. Route 3: the exact three-alternative density

`dev/gddm-nchoice/images.R`, `probe1-images.R`, `probe2-tape.R`.

### 5.1 The construction

The group generated by reflections in the three sides of an equilateral
triangle tiles the plane and the triangle is its fundamental domain, so
the Dirichlet heat kernel is the ALTERNATING sum of free heat kernels
over the orbit of the source:

    q(z, t | z0) = sum_g det(g) G_t(z - g z0),  G_t(w) = e^{-|w|^2/2t} / 2 pi t

Drift is a Girsanov tilt, `p = exp(<a, z - z0> - |a|^2 t/2) q`, exact
because q vanishes on the boundary and so the tilt does not disturb the
boundary condition. The defective density for choice k is the outward
flux through side k, `-(1/2) int_side d_n p`, and the integral along the
side is Gaussian, so it closes:

    f_k(t) = e^{-<a,z0> + r beta_k - beta_k^2 t / 2} / (2 t sqrt(2 pi t))
             * sum_g det(g) (r - h_g) e^{-(r - h_g)^2 / 2t + alpha_k b_g}
               * [ Phi((L/2 - b_g - alpha_k t) / sqrt(t))
                   - Phi((-L/2 - b_g - alpha_k t) / sqrt(t)) ]

with r the inradius, `L = 2 sqrt(3) r` the side length, `beta_k` and
`alpha_k` the normal and tangential components of the drift at side k,
and `h_g`, `b_g` the normal and tangential coordinates of image g. Every
quantity is elementary.

**The group elements are data.** `R_k(z) = (I - 2 n_k n_k') z + 2 r n_k`,
so every element is `g(z) = A_g z + r w_g` with `A_g` and `w_g`
independent of every parameter. The boundary r, the start z0 and the
drift a all enter the closed form smoothly. That is what makes it
tapeable, and it is the same discipline `gddm()` follows when it refuses
to let an integer index depend on a parameter.

### 5.2 It is right

Five checks, four of them against something that does not go through the
series at all.

**The kernel vanishes on the walls.** Worst `|q|` over nine points
across the three sides, at t = 0.3 and t = 1.5: **6.18e-17**.

**The flux is the derivative of the survival.** Integrating q over the
triangle by grid quadrature and differencing:

| t | S(t) | -dS/dt | sum_k f_k(t) | rel |
| --- | --- | --- | --- | --- |
| 0.10 | 0.977176 | 0.869667 | 0.871600 | 2.2e-03 |
| 0.30 | 0.652814 | 1.707951 | 1.708281 | 1.9e-04 |
| 0.60 | 0.280336 | 0.813903 | 0.813805 | 1.2e-04 |
| 1.00 | 0.087235 | 0.255088 | 0.255051 | 1.4e-04 |
| 1.50 | 0.020216 | 0.059128 | 0.059120 | 1.4e-04 |

The residual is the quadrature's, not the series'.

**The mean first-passage time matches the torsion function.** At zero
drift the mean first-passage time solves `(1/2) Lap u = -1` with u = 0 on
the walls, and on an equilateral triangle that has a closed form: the
three side functions `l_k = r - <n_k, z>` are affine with unit gradients
meeting at 120 degrees and summing to 3r, so `Lap(l1 l2 l3) = -3r`
exactly, and `u = 2 l1 l2 l3 / 3r` is the answer.

| case | mass | P1 | P2 | P3 | mean | exact mean |
| --- | --- | --- | --- | --- | --- | --- |
| zero drift, centre, r = 0.866 | 0.999590 | 0.332859 | 0.332859 | 0.333872 | 0.4860 | 0.4999904 |
| zero drift, offset, r = 0.866 | 0.999587 | 0.399673 | 0.192976 | 0.406938 | 0.4482 | 0.4622536 |
| drift (1.2, -0.4), r = 0.866 | **1.000000000** | 0.532267 | 0.079671 | 0.388062 | 0.4262 | - |
| drift (2.5, 1.0), r = 0.866 | **1.000000000** | 0.944691 | 0.016126 | 0.039183 | 0.2250 | - |
| zero drift, centre, r = 1.6 | 1.000000328 | 0.333333 | 0.333333 | 0.333333 | 1.7066670 | **1.7066667** |

Read the last row: exactly 1/3 per choice, and the mean agrees with the
closed form to **2e-7 relative**. Read the first two rows as the
truncation failure of section 5.3, not as an error in the construction.
At r = 0.866 a zero-drift case has a long tail in the dimensionless time
`t / r^2`, and a 400-tile table does not reach the end of it. The drift
rows have no such tail and come out at exactly 1.000000000.

**Symmetry.** At zero drift from the centre the three choice
probabilities must be 1/3 each and the three densities identical at
every t. Measured directly on the flux, before any integration:
identical to eight decimals at t = 0.1, 0.3 and 1.0.

**Against simulation.** Euler-Maruyama with a Brownian-bridge absorption
correction. Plain Euler misses excursions between monitoring points and
reports first passages late by O(sqrt(dt)); conditional on the two
endpoints the normal component is a Brownian bridge whose probability of
touching a flat wall is `exp(-2 d0 d1 / dt)` in closed form, and
absorbing with that probability removes the leading bias and lets the
step be 25 times coarser. `mu = (0.9, 0.1, -1.0)`, `c = 0.7071`, 100000
trials, step 1e-3:

| choice | series P | sim P | se | z |
| --- | --- | --- | --- | --- |
| 1 | 0.616589 | 0.617120 | 0.001537 | 0.35 |
| 2 | 0.303962 | 0.303650 | 0.001454 | -0.21 |
| 3 | 0.079448 | 0.079230 | 0.000854 | -0.26 |

The per-choice response-time quantiles, three choices by the 0.1, 0.3,
0.5, 0.7 and 0.9 quantiles, have |z| at most 1.70 over the fifteen
comparisons and a mean |z| of 0.80. Refining the step to 2.5e-4 does not
move them systematically. Series and simulation agree within Monte Carlo
error everywhere.

### 5.3 Where it fails, and it does fail

The image sum is an alternating sum of Gaussians. At decision time t the
images that still matter reach out to about `6 sqrt(t)`, so the number
of tiles a horizon needs grows like the area of that disc over the area
of the triangle: **linearly in the dimensionless time `tau = t / r^2`**.
Past the truncation the sum is not merely inaccurate, it is garbage,
because the terms that would have cancelled are missing. That is what
spoils the two zero-drift rows in section 5.2.

Measured: the largest t at which an m-tile sum is within 1e-8 of the
400-tile sum.

| tiles | r = 0.600 | r = 0.866 | r = 1.600 | tau = t / r^2 |
| --- | --- | --- | --- | --- |
| 22 | 0.170 | 0.355 | 1.207 | 0.47 |
| 46 | 0.497 | 1.036 | 3.523 | 1.38 |
| 76 | 0.889 | 1.853 | 6.301 | 2.47 |
| 118 | 1.207 | 2.516 | 8.556 | 3.35 |
| 166 | 1.590 | 3.314 | 11.269 | 4.41 |
| 250 | 2.295 | 4.784 | 16.268 | 6.37 |
| 340 | 2.844 | 5.927 | 20.154 | 7.89 |

The three r columns collapse onto one `tau` to three significant figures
(0.472, 0.473, 0.471 in the first row), which confirms the scaling law
rather than assuming it, and `tau_max` is about `tiles / 43`.

What that costs in practice: from the S(t) table above the survival
falls about as `exp(-1.6 tau)`, so 250 tiles reaches down to a survival
of about 2e-5, five decades below the density's peak. That is a MILDER
limit than the shipped `gddm()`'s own, which at its first grid node is
1.0e13 times the closed-form truth (`gddm-plots-findings.md`). But it is
a real one, and it is the mirror image of the limit `wiener()` already
lives with: there the small-time series fails at large t, and Navarro
and Fuss pair it with a large-time series.

**The pair exists here too and was NOT built.** Poisson summation over
the translation lattice turns the image sum into an eigenfunction
expansion `sum_K c_K exp(-|K|^2 t / 2)` over the dual lattice, which
converges fastest exactly where the images converge slowest. This probe
did not implement it. What a follow-up would have to settle: the term
count at a given tolerance; whether the crossover between the two series
can be chosen without comparing on a parameter, which is the trap
`gddm()` already documents and solves with a smooth weight; and whether
the side integral of a complex plane wave stays well conditioned.

### 5.4 Cost

`probe2-tape.R`. One condition, no grid anywhere.

| trials | tiles | tape build | one `fn` | one gradient |
| --- | --- | --- | --- | --- |
| 500 | 46 | 0.80 s | 12 ms | 28 ms |
| 500 | 118 | 0.94 s | 16 ms | 34 ms |
| 500 | 250 | 1.61 s | 12 ms | 28 ms |
| 2000 | 46 | 0.84 s | 22 ms | 28 ms |
| 2000 | 118 | **1.97 s** | 28 ms | **50 ms** |
| 2000 | 250 | 3.44 s | 50 ms | 152 ms |
| 8000 | 46 | 1.63 s | 24 ms | 60 ms |
| 8000 | 118 | 5.61 s | 134 ms | 244 ms |
| 8000 | 250 | 12.51 s | 170 ms | 360 ms |

Against the shipped `gddm()`'s 8.56 s and 6.7 ms at 2000 trials: **the
tape builds 4.3 times faster and the gradient runs 7.5 times slower**.

That trade has a structural cause worth stating plainly. **The closed
form costs per TRIAL; the grid costs per CONDITION.** `gddm()` solves
one PDE per distinct parameter setting and every trial in that condition
then reads a gather, so its cost is flat in the trial count. The image
series evaluates a sum over images at each trial's own response time. In
a design with few conditions and many trials the grid amortizes and the
closed form does not; with many conditions it is the other way round.

The obvious fix, evaluating the closed form on a shared time grid per
condition and interpolating, buys the grid's amortization back at the
price of reintroducing an interpolation error, and would then be
strictly better than route 2 at n = 3 on both cost and accuracy. It was
not prototyped.

The truncation shows up in the likelihood as well as the density: at
2000 trials the negative log-likelihood is 1424.218 at 46 tiles against
1424.226 at 118 and 250, so 46 tiles is already wrong in the third
decimal of the objective.

Note what disappears along with the grid. The density is evaluated at
`t - ndt` directly, so the non-decision time needs no B-spline
convolution and its derivative is exact. There is no density floor, no
renormalization to undo a mass loss that depends on the parameters, and
no leading-edge region where an implicit scheme puts mass the truth does
not have. Three of `gddm()`'s standing compromises simply do not arise.

### 5.5 The tape is clean

At 200 trials and 118 tiles, against numDeriv (Richardson):

| parameter | RTMB | numDeriv | rel |
| --- | --- | --- | --- |
| mu12[1] | -15.72192925 | -15.72186003 | 4.4e-06 |
| mu12[2] | -10.38520595 | -10.38513697 | 6.6e-06 |
| logc | -10.00311757 | -10.00289917 | 2.2e-05 |
| z0[1] | -27.97887553 | -27.97957590 | 2.5e-05 |
| z0[2] | -27.88175819 | -27.88120601 | 2.0e-05 |
| ndt | -46.23786005 | -46.23788762 | 6.0e-07 |

The residual is numDeriv's, not RTMB's. Second order tapes too: RTMB's
Hessian matches `numDeriv::jacobian` of the gradient to **7.63e-07**
relative, so Laplace over random effects and `sdreport()` would both
work, which is what a family in this package has to support.

And it fits. 4000 trials at 118 tiles, `nlminb` from the truth,
converged in 51 iterations and 34.3 s, with `sdreport()` returning
finite standard errors for all six parameters:

| parameter | estimate | truth | se | z |
| --- | --- | --- | --- | --- |
| mu12[1] | 0.83695 | 0.90000 | 0.03198 | -1.97 |
| mu12[2] | 0.13421 | 0.10000 | 0.03394 | 1.01 |
| logc | -0.35612 | -0.34658 | 0.00827 | -1.15 |
| z0[1] | -0.01102 | 0.00000 | 0.01091 | -1.01 |
| z0[2] | 0.01956 | 0.00000 | 0.01622 | 1.21 |
| ndt | 0.25226 | 0.25000 | 0.00286 | 0.79 |

Every parameter is within two standard errors of the truth, including
both components of the start-point bias, which is the one the
two-dimensional geometry adds and which has no counterpart in
`wiener()`. `nlminb` reported one `NA/NaN function evaluation` on the
way, at a trial step the truncated series could not represent; it
recovered and converged, which is the same behaviour `gddm()`'s density
floor exists to produce.

## 6. Why route 3 still probably should not be built

The construction is exact, cheap and correct. The case against it is
about what it reaches, not whether it works.

**Three alternatives, and structurally no more.** Section 2. A user with
a four-choice task gets nothing at all, and there is no path from three
to four.

**One stopping rule, chosen because it is solvable.** "Beat the mean of
the others by c" gives a simplex. The rule more often written down,
"beat the BEST competitor by c", gives a hexagon at n = 3, and a regular
hexagon is not a planar reflection domain: the only ones are the
rectangle and the 45-45-90, 30-60-90 and equilateral triangles. So that
rule has no closed form. At n = 2 every one of these rules coincides,
which is why the ambiguity has never had to be faced. Shipping the
solvable rule and calling it the multi-alternative diffusion would be a
claim the geometry had made, not the psychology.

**Equal thresholds only.** Unequal thresholds per alternative, which is
how a response bias would be modelled, make the triangle scalene and
destroy the reflection group. Start-point bias is fine and is a free
two-vector; threshold bias is not available at all.

**It is not a `gddm()`.** Constant drift only, fixed boundaries only. It
could accept none of `gddm_drift_leak()`, `gddm_drift_coherence()`,
`gddm_bound_exponential()`, `gddm_bound_linear()` or
`gddm_start_uniform()`. Collapsing boundaries and state-dependent drift
are exactly what `gddm()` exists for, and they are exactly what breaks
the images. This would be a new family sitting next to `wiener()`, with
its own component-free surface, its own documentation, and a name that
would have to explain why it takes none of `gddm()`'s arguments.

**Demand.** `lba(n)` and `rdm(n)` already fit n-alternative designs.
Section 3's separating manipulation, raising overall evidence while
holding contrasts fixed, is not a design that appears in the standard
multi-alternative literature. The models the race cannot express are
real, and section 3 proves it, but they are not models this package has
been asked for.

## 7. Complex AD: available, and not what this problem needs

RTMB 1.9 carries `adcomplex` with `fft`, `Re`, `Im`, `Mod`, `Arg`,
`Conj` and a complex `solve`, and nothing in frmtmb uses it. Priced here
and not used, for three reasons.

**`fft` has nothing to convolve.** The natural use would be the
non-decision-time convolution, and section 5.4 removes it: the closed
form is evaluated at `t - ndt` directly with an exact derivative, which
is strictly better than any convolution. The other candidate, a spectral
solve of the Fokker-Planck equation, wants a periodic domain, and an
absorbing boundary is the opposite of periodic.

**The spectral object that does exist is real.** The eigenfunction
series of section 5.3 is naturally written as an alternating sum of
complex plane waves over the six-element point group, but those sums
collapse to real trigonometric combinations, and a complex multiply is
two real tape nodes either way. Complex arithmetic would be a
convenience in writing it, not a capability that unlocks it.

**A transform in TIME does not close.** The Laplace transform of the
first-passage density is elementary in one dimension, which is where
`wiener()`'s two series come from. On the triangle the resolvent of the
Dirichlet Laplacian is not elementary, and inverting it numerically
means a contour integral whose accuracy would have to be established
from scratch, on the tape, at every parameter value the optimizer
visits.

If `adcomplex` gets a first use in this package it should be where the
frequency domain is the natural home of the model, which is what
`dev/frequency-domain-todo.md` is about, and not here.

## 8. What was NOT done

- **No package code was changed.** No file under any `R/`, no
  `DESCRIPTION`, no `NAMESPACE`, no test, no vignette, no `NEWS.md`.
  `git status` in the worktree shows only this document and
  `dev/gddm-nchoice/`.
- **The dual eigenfunction series of section 5.3 was not implemented.**
  The large-t failure of the image series is measured and its scaling
  law confirmed, but the series that would repair it is described, not
  built. That is the single largest gap in this probe, and it is the one
  that would have to close before route 3 could ship.
- **The regular tetrahedron was not built.** Section 4.4's
  four-alternative cost comes from the same ADI scheme on a cube at
  matched resolution and step count, plus the measured flat cost per
  node-step. Tetrahedral geometry would change the node count by a
  constant factor of about six and would not change the exponent, which
  is what the extrapolation turns on.
- **No fit to real data.** Every number here is against a closed form, a
  simulation, or another prototype.
- **The observable-terms comparison in `probe4-race.R` measurement 3
  produced no table.** It reads fitted coefficients back through the
  family's links to resimulate, and `rdm()`'s degenerate optimum
  (`v3` at -18.64 on the log link) is not a point that route survives.
  The log-likelihood gap and the per-trial nats in section 3.1 carry the
  argument instead, and they need no back-transformation. The script
  still prints the fitted coefficients, which is where the degeneracy
  was found.
- **The "beat the best competitor" hexagon was not solved.** It has no
  closed form (section 6), and route 2 would handle it at route 2's
  price; a second PDE domain was not worth building to say so.
- **`gd_check_response()`'s wording was not corrected**, though section
  3 shows one clause of it is wrong. Changing it is a package edit and
  this lane makes none.
- **No package was installed into any shared library.** The probes read
  frmtmb 0.53.0 and frmtmb.eam 0.4.0 from the user library and installed
  nothing.

## Probe scripts

| file | what |
| --- | --- |
| `dev/gddm-nchoice/common.R` | the n-alternative geometry, the Helmert basis, the bridge-corrected simulator |
| `dev/gddm-nchoice/images.R` | route 3: the reflection-group unfolding and the closed-form density |
| `dev/gddm-nchoice/pde2d.R` | route 2: the two-dimensional Douglas ADI solve on the triangle |
| `dev/gddm-nchoice/probe0-baseline.R` | what the shipped two-choice `gddm()` costs and how accurate it is |
| `dev/gddm-nchoice/probe1-images.R` | route 3 correctness: walls, survival, torsion function, symmetry, simulation |
| `dev/gddm-nchoice/probe2-tape.R` | route 3 on the tape: truncation horizon, gradient, Hessian, cost, recovery |
| `dev/gddm-nchoice/probe3-pde.R` | route 2: convergence, the ghost correction, symmetry, cost, gradient |
| `dev/gddm-nchoice/probe4-race.R` | route 1: the separating manipulation and the race fits |
| `dev/gddm-nchoice/probe5-cost.R` | cost per node-step in one, two and three dimensions |
