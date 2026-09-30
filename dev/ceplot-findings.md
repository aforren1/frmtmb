# Lane wt-ceplot: plot(), conditional_effects() and small post-fit methods

Base: 57c25589 (frmtmb 0.66.0, frmtmb.sample 0.14.0), measured against
the read-only build `C:/Users/adf44/source/r/rellib-r4`. Lane library:
`C:/Users/adf44/source/r/wt-ceplot-lib` (frmtmb and frmtmb.sample from
this worktree; the other six extensions load from rellib-r4 on top of
the lane's frmtmb). brms 2.23.0's source for every function ported
here was dumped from the installed namespace into
`dev/ceplot-brmssrc/` by `dev/ceplot-brmssrc.R` and
`dev/ceplot-brmssrc2.R`; the helpers `use_alias()`, `limit_chars()`,
`replace_args<-`, `make_point_frame()` and `extract_pars()` were read
with `dev/ceplot-brmssrc3.R` and `dev/ceplot-brmssrc4.R`.

## 1. What changed

### 1.1 plot() of a conditional-effects object (rows :154, :162, :164, :169)

`plot.frmtmb_conditional_effects()` now has exactly the formals of
brms's `plot.brms_conditional_effects()`, in brms's order: `x, ncol,
points, rug, mean, jitter_width, stype, line_args, cat_args,
errorbar_args, surface_args, spaghetti_args, point_args, rug_args,
facet_args, theme, ask, plot, ...`. A test pins the list against brms's
own formals when brms is installed and against a literal copy when it
is not. What each does here:

| argument | brms | frmtmb |
|---|---|---|
| `ncol` | facet columns | facet columns, as before |
| `points` | observations over the curves | as before; with a `conditional_smooths()` result it is brms's error "Argument 'points' is invalid for objects returned by 'conditional_smooths'." (it was a message) |
| `rug` | `geom_rug()` of the first predictor's observed values | `graphics::rug()` of the same values, from the point frame `points` uses, for a numeric first predictor |
| `mean` | with spaghetti, `FALSE` drops the estimate line | the same; without spaghetti the line is always drawn, as in brms |
| `jitter_width` | deprecated, with a warning, then used as the jitter width | the same warning, then the horizontal spread of the observations |
| `stype` | `"contour"` or `"raster"` for a surface | colored contour lines, or `graphics::image()` |
| `*_args` | passed to the ggplot2 layer | translated to base graphics; see below |
| `theme` | a ggplot2 theme, or brms's error | `NULL`; a ggplot2 theme is accepted with a warning that it is not applied; anything else is brms's error |
| `ask` | prompt after the first page | the same, brms's page rule |
| `plot` | `FALSE` returns the ggplot objects undrawn | `FALSE` returns the plot objects undrawn (1.2) |
| `...` | `do_plot` is the alias of `plot`, with a warning | the same; any other name is an error that names it, except base graphical parameters, which `plot()` methods tolerate |

The layer lists: `colour`/`color`/`col`, `fill`, `alpha`,
`linewidth`/`lwd`, `linetype`/`lty`, `shape`/`pch` and `size` (point
size for points, line width for lines) are translated. In `line_args`,
`fill` and `alpha` color the band, as they color the ribbon of brms's
`geom_smooth()`. `errorbar_args$width` is the cap width in axis units
(ggplot2's meaning); `point_args$width` spreads the observations by up
to that much on either side, deterministically (no RNG, so a replot is
identical and the caller's seed is untouched). `rug_args$sides` takes
`"b"` and `"t"`; `"l"` and `"r"` draw nothing, because the rug maps
the x values only, which is what ggplot2 draws for brms's layer too.
`surface_args$bins` is the number of contour levels. `facet_args` takes
`ncol`, `nrow` and `scales` (`"free"` and `"free_y"` give each panel
its own vertical range). brms's own refusals of a malformed list are
kept verbatim: "Argument 'line_args' must be named." and "Argument(s)
mapping cannot be replaced.". Any element with no counterpart is named
in ONE warning ("plot() ignores line_args$stat, cat_args$position:
..."), because a setting the caller asked for and did not get would
otherwise pass in silence.

A call that sets any layer list, `rug` or `mean = FALSE` draws a
several-condition display with the base panel grid rather than
tinyplot, because the tinyplot layer has no slot for those settings.

Rendered for a visual check by `dev/ceplot-plots.R` into
`dev/ceplot-log/plot-*.png` (defaults; every layer argument; contour
and raster; spaghetti with and without the mean; a free-scale facet
row; a hypothesis page with a cut title).

### 1.2 What `plot = FALSE` returns

brms returns an invisible named list of ggplot objects, one per
effect; `plot = FALSE` skips drawing them. frmtmb does not depend on
ggplot2 and draws with base graphics, so it has no ggplot object to
return. The base-graphics equivalent of "an object that draws when it
is printed" is returned instead: an invisible named list of
`frmtmb_ce_plot` objects, one per effect, each holding the effect's
data frame and the settings of the call, with `print()` and `plot()`
methods that draw it. Returned whether or not the call draws, as brms
returns its list. The alternative, `grDevices::recordPlot()` objects,
was not taken: a recorded plot has to be drawn to be recorded, is not
safe to serialize across R versions, and holds the pixels of one
device rather than the display. Neither can be extended with `+` as a
ggplot object can, so `expect_ggplot(plot(me, plot = FALSE)[[1]])` in
brms's suite still does not hold (rows :154, :162, :164, :169, :391,
:396 stay not held, now on the class alone; section 5).

### 1.3 plot() of a hypothesis (rows :391, :396)

`plot.frmtmb_hypothesis()` has brms's formals `x, nvariables, N,
ignore_prior, chars, colors, theme, ask, plot, ...`. Hypotheses are
stacked `nvariables` (default 5) to a page in one column, as brms's
`facet_wrap(ncol = 1)` lays them out; they were one per page. `N` is
brms's alias with brms's warning. `chars` cuts a long title exactly as
brms's `limit_chars()` does (four brms outputs pinned in a test,
`dev/ceplot-limitchars.R`). `colors` must be two colors, brms's error
otherwise; the first fills the histogram. `ignore_prior` is accepted
and changes nothing: a fit has no prior, and frmtmb.sample does not
draw from the prior, so `prior_samples` is all `NA`; brms draws the
posterior alone in that case too (`length(unique(samples$Type)) == 1`
after `na.omit()`). `plot = FALSE` returns a list of `frmtmb_hyp_plot`
objects, one per page, which draw when printed.

### 1.4 conditional_effects(): invalid effects (row :205, and :203)

brms's rule, read from `conditional_effects.brmsfit()`: an effect is
valid when every variable of it is a predictor of some distributional
or nonlinear parameter or a grouping variable (`get_all_effects(...,
comb_all = TRUE)`), and a pair names two different variables. One
invalid effect among valid ones is dropped with brms's warning; with
none valid, brms's error. frmtmb now does the same, with brms's
quoting: "Some specified effects are invalid for this model:
'xx'\nValid effects are (combinations of): 'x', 'f'". The valid set is
every variable the model reads outside the response (`ce_locked_vars()`,
the `cs()` variables and the grouping variables). A variable present in
the data but not in the model is invalid, as in brms; it drew a flat
curve before. The order-three check still comes first, as in brms.
brms 2.23.0's answers on the same calls: `dev/ceplot-log/brms-mi.txt`.

### 1.5 conditional_effects(): the ordinal warning (row :217)

brms warns "Predictions are treated as continuous variables in
'conditional_effects' by default which is likely invalid for ordinal
families. Please set 'categorical' to TRUE." whenever it draws the
expected category number of an ordinal fit, which is its default.
frmtmb draws that summary only under `categorical = FALSE` (its
default is the per-category display, the one brms's warning asks for;
a deliberate divergence since 0.60, `?conditional_effects`, "Ordinal
responses"). The warning is now given there, in `ce_display_kind()`,
which the fit and draws methods share, worded for the case: "...
with categorical = FALSE, which is likely invalid ...". The default
call and `dpar = "mu"` give none, and a test holds both. brms warns once
per effect; frmtmb once per call. Row :217 asserts the warning on the
DEFAULT call, so it does not hold. The user decided on 2026-09-30 to
keep frmtmb's per-category default, the display brms's own warning
recommends, as a deliberate divergence from brms; the row is
re-verdicted at consolidation.

### 1.6 conditional_effects() of mi(x, idx = ) (backlog)

brms 2.23.0 does not answer this display: on the same model every call
stops with "Index of response 'x' contains duplicated values."
(`dev/ceplot-brms-mi.R`), because the grid holds `index()` at one
value in every row. 0.66.0 stopped inside the prediction with "Could
not match all indices". It is now refused by name with the reason
(`ce_mi_idx_check()`): "conditional_effects() cannot display response
'y': its predictor reads mi(x, idx = g1), so each row reads x on the
row whose index(g2) equals its g1. The display grid holds both at one
reference value, so every grid row names the same row and the index
cannot be matched. brms 2.23.0 stops on the same grid ...". A response
that does not read the term (`resp = "x"` of the same model) and
`mi(x)` without `idx` are not refused, and a test holds both. Making
the display work would need a rule for which row of `x` a synthetic
grid row reads, which neither package has.

### 1.7 New levels: crossed terms and mm(by = ) (backlog)

Read first: `dev/postfit2-findings.md` sections 7c to 7e. The
bootstrap and draws bands place a new level by moving grouping columns
to a placeholder level of the term. For crossed
`(1 | g) + (1 | h) + (1 | g:h)` with `g` and `h` at observed levels
never seen together, every placeholder of `g:h` moves `g` or `h`, which
the observed terms read, so the call was refused.

`ce_plan_part()` now falls back, only where no placeholder serves, to
RENAMING instead of moving: the rows keep their columns, and in the
copy of the fit used for those rows alone (`ce_plan_eval()`) the
term's first level is renamed to the label the design gives the rows
(`ce_design_label()`, the grouping expression evaluated on the row,
`"1:1"`). The design matches level labels by string
(`match(gv, bk[["levels"]])` in `pred_design()`), so the rows read that
slot, whose coefficients the draw overwrites; nothing else moves, and
the fit itself is untouched (a test checks its levels afterwards).

Against brms at the same draws, `dev/ceplot-crossed-brms.R` (data seed
49, the P1-M2 case of `dev/postfit2-p2-checks.R`; 200 hand-built draws,
seed 1; brms at `algorithm = "fixed_param"`, one chain per draw,
initialized at the draw; curves at seeds 1 to 3):

```
frmtmb from C:/Users/adf44/source/r/wt-ceplot-lib/frmtmb 
fit sd: g 0.8431 h 1.3912 g:h 0.7179
theta order check (block sd at the estimate): 0.8431 1.3912 0.7179
brms draws 200; r_g:h max diff 1.11e-16; sd_g:h max rel diff 0
observed g = 1, h = 2 (no draw): max diff estimate 1.11e-16 lower 2.22e-16 upper 0
unseen g:h, seed 1: estimate brms 0.7882 1.9584 3.1750 | frmtmb 0.7882 1.9584 3.1750
unseen g:h, seed 1: lower brms -0.5651 0.6918 1.8801 | frmtmb -0.5651 0.6918 1.8801; upper brms 2.3347 3.5997 4.7778 | frmtmb 2.3347 3.5997 4.7778
unseen g:h, seed 1: max diff estimate 0 lower 2.22e-16 upper 4.44e-16
unseen g:h, seed 2: estimate brms 0.7511 2.0092 3.1840 | frmtmb 0.7511 2.0092 3.1840
unseen g:h, seed 2: lower brms -0.7615 0.4867 1.7645 | frmtmb -0.7615 0.4867 1.7645; upper brms 2.3044 3.4632 4.6758 | frmtmb 2.3044 3.4632 4.6758
unseen g:h, seed 2: max diff estimate 4.44e-16 lower 2.22e-16 upper 0
unseen g:h, seed 3: estimate brms 0.7989 2.0225 3.2216 | frmtmb 0.7989 2.0225 3.2216
unseen g:h, seed 3: lower brms -0.5162 0.7059 1.9167 | frmtmb -0.5162 0.7059 1.9167; upper brms 2.1800 3.4764 4.5554 | frmtmb 2.1800 3.4764 4.5554
unseen g:h, seed 3: max diff estimate 0 lower 2.22e-16 upper 0
fit unseen g:h | wald: est 0.8020 1.9986 3.1951 | width 3.4602 3.4398 3.4848
fit observed g:h | wald: est -0.3669 0.8297 2.0263 | width 0.9148 0.8463 1.0237
fit unseen g:h | boot60: est 0.8020 1.9986 3.1951 | width 3.5397 3.5828 3.6507
fit unseen g:h: boot width / wald width 1.0230 1.0416 1.0476
draws unseen g:h: width 2.8997 2.9079 2.8977 (0.8380 0.8454 0.8315 of the Wald width); draws observed g:h: width 0.6587 0.4249 0.6559
done
```

An `mm(g1, g2, by = cbind(f1, f2))` term is one block per by-level. A
new member is now placed at a free level of the block of its own
by-value, where the other by-levels' blocks read it as an observed group
of another by-level and skip it (`by_route_rows()`); two new members
with the same value and by-value share one placeholder, so their weights
add. `ce_mm_in_block()` routes the members per block, and an observed
block of such a term holds only its own observed members. With one
member at an observed group of by-level `a` and the other new in `b`,
each block is read one way only, so the bootstrap answers too (the
same row on a term without `by` stays refused, as before). No brms
comparison is possible: brms 2.23.0 cannot draw a new level of an
`mm(by = )` term at all, and stops with one of three errors
(`dev/ceplot-rev-log/mmby-brms.txt`; this record said "no brms
comparison" for another reason before punch round 1, corrected). The
check is against the Wald band, which
carries each new member's variance from its own block,
`dev/ceplot-mmby.R` (data seed 45, bootstrap seed 5, 200 hand-built
draws):

```
ARM lane frmtmb from C:/Users/adf44/source/r/wt-ceplot-lib/frmtmb 
block sd at the estimate (by-levels a, b): 0.3666 0.7410
both unset, f1 = a, f2 = b | wald: est 0.2738 1.9968 3.7198 | width 1.7606 1.7335 1.7835
both unset, f1 = a, f2 = b | boot60: est 0.2738 1.9968 3.7198 | width 1.4296 1.4040 1.4790 | width / wald 0.8120 0.8099 0.8293
both unset, f1 = a, f2 = b | draws: est 0.3208 2.0285 3.7307 | width 1.5396 1.5200 1.6676 | width / wald 0.8745 0.8769 0.9350
both unset, f1 = f2 = a | wald: est 0.2738 1.9968 3.7198 | width 1.5936 1.5636 1.6188
both unset, f1 = f2 = a | boot60: est 0.2738 1.9968 3.7198 | width 1.4145 1.3570 1.3559 | width / wald 0.8876 0.8679 0.8376
both unset, f1 = f2 = a | draws: est 0.2801 1.9826 3.7358 | width 1.6023 1.5322 1.8370 | width / wald 1.0055 0.9799 1.1348
g1 = 99, g2 = 98, f1 = a, f2 = b | wald: est 0.2738 1.9968 3.7198 | width 1.7606 1.7335 1.7835
g1 = 99, g2 = 98, f1 = a, f2 = b | boot60: est 0.2738 1.9968 3.7198 | width 1.4296 1.4040 1.4790 | width / wald 0.8120 0.8099 0.8293
g1 = 99, g2 = 98, f1 = a, f2 = b | draws: est 0.3208 2.0285 3.7307 | width 1.5396 1.5200 1.6676 | width / wald 0.8745 0.8769 0.9350
g1 = 2 (seen, by a), g2 unset, f2 = b | wald: est 0.5797 2.3027 4.0257 | width 1.5361 1.5065 1.5653
g1 = 2 (seen, by a), g2 unset, f2 = b | boot60: est 0.5797 2.3027 4.0257 | width 1.2750 1.2868 1.4879 | width / wald 0.8300 0.8541 0.9505
g1 = 2 (seen, by a), g2 unset, f2 = b | draws: est 0.5827 2.2935 4.0352 | width 1.6193 1.5286 1.7988 | width / wald 1.0541 1.0147 1.1492
both seen, g1 = 2, g2 = 7 | wald: est 0.4552 2.1782 3.9012 | width 0.4463 0.3292 0.5368
both seen, g1 = 2, g2 = 7 | boot60: est 0.4552 2.1782 3.9012 | width 0.3990 0.3038 0.4667 | width / wald 0.8940 0.9229 0.8694
both seen, g1 = 2, g2 = 7 | draws: est 0.4577 2.1808 3.9258 | width 0.5780 0.2589 0.8139 | width / wald 1.2951 0.7864 1.5162
done
```

The observed-member rows are unchanged from 0.66.0 (the "both seen"
lines match `dev/ceplot-log/mmby-base.txt` exactly).

WITHDRAWN in punch round 1: this record explained the bootstrap's 0.81
to 0.95 at 60 refits as the shrinkage of a standard deviation estimated
from five groups. It was Monte Carlo noise. At 200 refits and seeds 5
to 7 the ratio is 0.83 to 1.11 over the three new-member conditions
(`dev/ceplot-rev-bootratio.R`, `dev/ceplot-rev-log/bootratio.txt`).

### 1.8 fitted(sample_new_levels = "old_levels") (row :314)

brms's `get_new_rdraws()` gives each unseen level one seen level of its
group, `sample(seq_len(nlevels), 1)` once per call, and reads that
level's draws. `predict()` already did this per replicate. `fitted()`
is a Wald summary, so the seen level is substituted in the DESIGN:
`new_level_pick_apply()` maps an unseen row to the chosen level in
`pred_design()` and `mm_newdata_parts()`, so the estimate is that
level's fitted value and `Est.Error` carries that level's conditional
variance from the joint precision, not the block's marginal variance.
The choice is made by `old_level_pick_add()`, which `predict()`'s
`predict_new_level_spec()` now calls too, over every linear predictor
of the response (`fitted_old_levels()`), so the two methods choose the
same level at the same seed. `dev/ceplot-oldlevels.R` (data seed 3,
calls at seeds 1 to 12; brms at 50 draws around the fit):

```
frmtmb from C:/Users/adf44/source/r/wt-ceplot-lib/frmtmb 
fitted: seeds 1..12 read level 9 5 5 8 2 5 10 4 11 11 10 2 (identical() to fitted at it)
Est.Error: gaussian 1.6713 1.6760 | old_levels at seed 1 0.3599 0.3867
predict: seeds 1..12 read level 9 5 5 4 2 5 10 8 11 11 10 2
seen levels (fixef + ranef): -1.404 2.075 1.942 0.153 2.520 -0.570 1.054 0.155 -1.910 3.308 2.291 0.887
predict and fitted pick the same level at 10 of 12 seeds
predict's own choice (predict_new_level_spec): seeds 1..12 level 9 5 5 8 2 5 10 4 11 11 10 2
predict's choice and fitted's agree at 12 of 12 seeds
brms: seeds 1..12 read level 9 5 5 8 2 5 10 4 11 11 10 2
brms and frmtmb fitted pick the same level at 12 of 12 seeds
done
```

(`predict()`'s draw means cannot tell levels 4 and 8 apart, 0.002 apart
against a Monte Carlo error of about 0.007; its own recorded choice is
the comparison.) Row :314 still does not hold: its call also passes
`ndraws = 10`, which a fit refuses (D3).

CORRECTED in punch round 1 (review B1): the choice above was made per
BLOCK, keyed by the block's first coefficient, where brms makes it per
GROUPING FACTOR. The one-block model of this script could not show it.
On `bf(y ~ x + (1 | g), sigma ~ (1 | g))` and `(1 + x || g)` the new
group read different seen groups for different terms, at 0 of 12
seeds, and `predict(old_levels)` had the same per-block rule since
0.66.0. See section 10.

### 1.9 posterior_samples(pars = ) order (row :635)

brms's `extract_pars()` keeps each pattern's matches in `variables()`
order, and brms's `variables()` lists every population-level intercept
first. frmtmb's `variables()` lists each predictor's coefficients
together. `ps_select()` puts the `b_` coefficients in `fixef()`'s order,
which is brms's, before matching, and reuses `draws_extract_pars()`.
(CORRECTED in punch round 1, review m5: `fixef()`'s order is not brms's
on a nonlinear model, where brms keeps a nonlinear parameter's
intercept with its own coefficients. `ps_select()` now applies brms's
rule directly; section 10.)
Without `pars` the columns stay those of `variables()`, so row :634
(`names(draws) == variables(fit1)`) keeps holding: reordering
`variables()` itself is the separate decision of
`dev/shapes-findings.md` section 8. A `pars` that matches nothing
returns `NULL`, as brms's does.

### 1.10 parnames() on a fit (row :995)

The generic moves from frmtmb.sample to frmtmb (`R/confint.R`), with
`parnames.frmtmb_fit()` and brms's warning, and joins frmtmb's owner
table (`parnames = "brms"`). frmtmb.sample drops its own generic and
owner-table entry, imports and re-exports frmtmb's, and keeps
`parnames.frmtmb_draws()`. With brms loaded, brms's generic warns too,
so a call warns twice there, as `parnames()` on draws already did.
(Fixed in punch round 1, review m6: both methods now warn once.)

### 1.11 nsamples(incl_warmup = TRUE) (row :595)

Decision: answer, not refuse. `frm_sample()` keeps the stanfit, and
rstan's record `@sim$n_save` counts the saved iterations with the
warmup (rstan saves the warmup unless `save_warmup = FALSE`). brms's
own count is `n_save[1] * chains` with `incl_warmup = TRUE`, and
`(n_save[1] - warmup2[1]) * chains` without; the draws matrix holds the
post-warmup draws only, but the count needs no stored warmup DRAWS,
only the record, so nothing `frm_sample()` stores changes. With
`save_warmup = FALSE`, `n_save` is the post-warmup count, as it is in
brms. Draws with no stanfit behind them refuse by name. brms's check of
`subset` against the count is added. The ported fixture samples
`iter = 75, warmup = 50`, so row :595's `75` holds.

## 2. Seen to fail on the base build

Every new test, run against rellib-r4 (`dev/ceplot-debug1.R <file>
<pkg> base`), logs `dev/ceplot-log/seen-failing-*.txt`:

```
seen-failing-ce-levels.txt:
  ── 1. Error ('test-ce-levels.R:347:5'): an mm() term with a by variable draws ea
  Error: conditional_effects() cannot draw a new level of the multi-membership term (1 | mm(g1, g2, by = cbind(f1, f2)) [cbind(f1, f2) = 1]) with a by variable on band = "boot" or 
  ── 2. Error ('test-ce-levels.R:385:3'): crossed terms at an unseen combination d
  Error: conditional_effects() cannot draw a new level of the group-level term (1 | g:h) at the grid rows with g = 1, h = 1: to place one it would move a column that another term o
seen-failing-ce-parity.txt:
  ── 1. Error ('test-ce-parity.R:25:3'): one invalid effect among valid ones is dr
  Error: Variable 'xx' is not stored in the model frame; pass the original data via data =
  ── 2. Failure ('test-ce-parity.R:74:3'): the expected category number of an ordi
  Expected `any(grepl(r, seen, fixed = TRUE))` to be TRUE.
  Differences:
  ── 3. Error ('test-ce-parity.R:98:3'): a display of mi(x, idx = ) is refused by 
  Error: Could not match all indices in response 'x' on newdata: 100 value(s) of g1 are not among the values of index(g2) on the rows of newdata that 'x' uses
seen-failing-ce-plot.txt:
  ── 1. Failure ('test-ce-plot.R:69:3'): plot() takes every argument brms's two pl
  Expected `ce_f` to be identical to `brms_ce`.
  Differences:
  ── 2. Failure ('test-ce-plot.R:70:3'): plot() takes every argument brms's two pl
  Expected `hyp_f` to be identical to `brms_hyp`.
  Differences:
  ── 3. Error ('test-ce-plot.R:82:3'): plot = FALSE returns one plot object per ef
  Error: plot() has no argument `plot`. Did you mean `plt`?. It takes: x, ask, points, ncol
  ── 4. Error ('test-ce-plot.R:102:3'): rug draws the observed values of a numeric
  Error: plot() has no argument `plot`. Did you mean `plt`?. It takes: x, ask, points, ncol
  ── 5. Error ('test-ce-plot.R:121:5'): mean = FALSE leaves the estimate out only 
  Error: plot() has no argument `plot` (and 1 more: mean). Did you mean `plt`?. It takes: x, ask, points, ncol
  ── 6. Error ('test-ce-plot.R:133:3'): stype picks contour lines or a raster for 
  Error: plot() has no argument `plot`. Did you mean `plt`?. It takes: x, ask, points, ncol
  ── 7. Error ('test-ce-plot.R:145:3'): layer arguments translate, and the rest ar
  Error: plot() has no argument `rug` (and 8 more: line_args, cat_args, errorbar_args, point_args, rug_args, spaghetti_args, surface_args, facet_args). It takes: x, ask, points, nc
  ── 8. Error ('test-ce-plot.R:175:3'): theme, jitter_width and points follow brms
  Error: plot() has no argument `theme`. It takes: x, ask, points, ncol
  ── 9. Error ('test-ce-plot.R:200:3'): facet_args lays out the panels of several 
  Error in `ce_facet_layout(4L, nrow = 1L)`: unused argument (nrow = 1)
  Backtrace:
  ── 10. Failure ('test-ce-plot.R:215:3'): plot() of a hypothesis takes brms's arg
  Expected `pages(plot(h, ask = FALSE))` to be identical to 1L.
  Differences:
  ── 11. Error ('test-ce-plot.R:216:3'): plot() of a hypothesis takes brms's argum
  Error: plot() has no argument `nvariables` (and 1 more: plot). It takes: x, ask
  ── 12. Error ('test-ce-plot.R:240:3'): a long hypothesis is cut in its title as 
  Error in `hyp_limit_chars(x, 20)`: could not find function "hyp_limit_chars"
  Backtrace:
seen-failing-parity-defects.txt:
  ── 1. Error ('test-brms-parity-defects.R:203:5'): fitted() takes sample_new_leve
  Error: fitted() honors sample_new_levels = "gaussian" only. An unseen level's effect is taken from its block's own estimated covariance, which is what brms's "gaussian" does. "ol
  ── 2. Error ('test-brms-parity-defects.R:231:3'): parnames() on a fit is variabl
  Error in `parnames(fit)`: could not find function "parnames"
  Backtrace:
seen-failing-postfit-draws.txt:
  ── 1. Error ('test-postfit-draws.R:382:5'): on draws, crossed terms at an unseen
  Error: conditional_effects() cannot draw a new level of the group-level term (1 | g:h) at the grid rows with g = 1, h = 1: to place one it would move a column that another term o
  ── 2. Error ('test-postfit-draws.R:420:5'): on draws, an mm() term with a by var
  Error: conditional_effects() cannot draw a new level of the multi-membership term (1 | mm(g1, g2, by = cbind(f1, f2)) [cbind(f1, f2) = 1]) with a by variable on band = "boot" or 
seen-failing-shapes-draws.txt:
  ── 1. Error ('test-brms-shapes-draws.R:87:3'): nsamples() is brms's, with brms's
  Error: nsamples(incl_warmup = TRUE) has nothing to count: frm_sample() discards the warmup rather than storing it, so the object carries post-warmup draws only. ndraws(x) is that
  ── 2. Failure ('test-brms-shapes-draws.R:138:3'): posterior_samples(pars = ) ord
  Expected `names(ps(pars = "^b_"))` to be identical to `brms_b`.
  Differences:
  ── 3. Failure ('test-brms-shapes-draws.R:152:3'): posterior_samples(pars = ) ord
  Expected `ps(pars = "^nothing_")` to be NULL.
  Differences:
seen-failing-summary.txt:
```

The weak form (a missing symbol) is the only one available for
`hyp_limit_chars()` and core's `parnames()`; the behavioural failure of
the title cut would need the drawn text, and `parnames(fit)` had no
method at all. Everything else fails on behaviour: an error where the
answer is now given, a missing warning, a wrong order, a wrong count.

Guard-absent cases, each in a test: valid effects (including a
grouping variable and a pair) give no warning; the default ordinal
display and `dpar = "mu"` give no warning; `resp = "x"` of the `mi(idx)`
model and plain `mi(x)` answer; a seen level under
`sample_new_levels = "old_levels"` is identical to the plain `fitted()`;
`mean = FALSE` without spaghetti draws the same operations as `TRUE`; a
factor predictor gets no rug; the per-panel page counter counts 0, 1, 2
and 3 pages in the same file.

## 3. Tests run

One file per R process, the package attached, `NOT_CRAN=true`
(`dev/ceplot-par.sh`, `dev/ceplot-runtest.R`). Final lane build, all
eight suites (`dev/ceplot-findings-gen.R`):

```
## full2: 329 files, 329 with a RESULT line
lane frmtmb: 195 files, expectations 13764, failed 0, error 0, skipped 163, warning 0, passed 13601
lane frmtmb.coupling: 11 files, expectations 547, failed 0, error 0, skipped 5, warning 0, passed 542
lane frmtmb.eam: 29 files, expectations 1746, failed 0, error 0, skipped 3, warning 0, passed 1743
lane frmtmb.latent: 10 files, expectations 361, failed 0, error 0, skipped 2, warning 0, passed 359
lane frmtmb.learn: 15 files, expectations 442, failed 0, error 0, skipped 13, warning 0, passed 429
lane frmtmb.ode: 11 files, expectations 548, failed 0, error 0, skipped 1, warning 0, passed 547
lane frmtmb.sample: 43 files, expectations 2242, failed 0, error 0, skipped 4, warning 0, passed 2238
lane frmtmb.spline: 15 files, expectations 554, failed 0, error 0, skipped 1, warning 0, passed 553

```

The skips are the gated files and platform skips; no file lacks a
RESULT line. Every log's `lib:` line shows frmtmb loaded from the lane
library, and each other extension from rellib-r4 on top of it.

Gated tier, `FRMTMB_BRMS_FIT_TESTS=true` and
`FRMTMB_DRMTMB_FIT_TESTS=true`, every gated file of core and
frmtmb.sample on the lane, and the two ported methods files on the
base. The five failures in the lane's two methods files are the five
rows of section 5, each a stale "now HOLDS" verdict; the gated tier
skips nothing:

```
## gated2: 21 files, 21 with a RESULT line
base frmtmb: 1 files, expectations 162, failed 0, error 0, skipped 0, warning 0, passed 162
base frmtmb.sample: 1 files, expectations 61, failed 0, error 0, skipped 0, warning 0, passed 61
lane frmtmb: 13 files, expectations 2314, failed 3, error 0, skipped 0, warning 0, passed 2311
lane frmtmb.sample: 6 files, expectations 489, failed 2, error 0, skipped 0, warning 0, passed 487
not clean:
  RESULT lane frmtmb test-brms-suite-methods.R: tests=162 failed=3 error=0 skipped=0 warning=0 passed=159
    fails: conditional_effects has reasonable ouputs
    fails: variables has reasonable ouputs
  RESULT lane frmtmb.sample test-brms-suite-methods.R: tests=61 failed=2 error=0 skipped=0 warning=0 passed=59
    fails: nsamples has reasonable ouputs
    fails: posterior_samples has reasonable outputs

```

## 4. R CMD check --as-cran

```
frmtmb: Status: 1 NOTE
  * checking HTML version of manual ... [20s] NOTE
  Skipping checking math rendering: package 'V8' unavailable
frmtmb.sample: Status: OK
```

The one NOTE is the V8 math-rendering note that lane-rules.md lists as
expected. Both were built and checked in `dev/ceplot-check/<pkg>/` by
`dev/ceplot-check.ps1`, with the lane library first in `R_LIBS`.

## 5. Ported brms-suite rows this change flips

Not edited in `dev/brmsport-verdicts.tsv` (lane rule); the gated files
show each as a stale failure, which is expected.

- frmtmb brmsfit-methods:203: is recorded as a pass in frmtmb's own words (Variable 'Trtc' is not stored in the model frame): STALE OWN-WORDS: brms's assertion holds as written 
- frmtmb brmsfit-methods:205: now HOLDS but is recorded as 'defect' (one invalid effect among valid ones is an error; brms warns 'Some specified effects are invalid' and returns the valid ones (dev/brmsport-probe4.R)). Change its row in dev/brmsport-verdicts.tsv and regenerate. 
- frmtmb brmsfit-methods:995: now HOLDS but is recorded as 'defect' (parnames() is brms's deprecated spelling of variables(). frmtmb.sample answers it on draws, with brms's warning; core has none for a fit, and adding one moves the generic from frmtmb.sample's owner table to core's (R/generic-owners.R)). Change its row in dev/brmsport-verdicts.tsv and regenerate. 
- frmtmb.sample brmsfit-methods:595: now HOLDS but is recorded as 'defect' (nsamples() answers now, as brms's does, but incl_warmup = TRUE is refused: frm_sample() discards the warmup rather than storing it, so there is nothing to count. The gap is dev/brms-api-diff.md (c), 'Blocked, not small'). Change its row in dev/brmsport-verdicts.tsv and regenerate. 
- frmtmb.sample brmsfit-methods:635: now HOLDS but is recorded as 'defect' (posterior_samples(pars = '^b_') answers now and returns the right SET of columns; their ORDER is frmtmb's variables() order, which lists each predictor's coefficients together, where brms lists every intercept first. fixef(), vcov() and summary()$fixed take brms's order at item 2.6f; variables() keeps its own, and reordering it is a separate decision (dev/shapes-findings.md section 8)). Change its row in dev/brmsport-verdicts.tsv and regenerate. 

On the base build both methods files hold every row as
recorded (dev/ceplot-log/gated2.log, the two `base` lines).

## 6. Deliberate divergences from brms, with the reason

- `plot = FALSE` returns base-graphics plot objects, not ggplot
  objects (1.2). frmtmb does not depend on ggplot2.
- A ggplot2 `theme` is accepted and not applied, with a warning; brms
  applies it. Base graphics has no theme object; `par()` is the
  counterpart.
- The ordinal default stays the per-category display, and brms's
  "treated as continuous" warning comes with `categorical = FALSE`
  (1.5). Reason: brms's own warning calls its default display likely
  invalid for ordinal families and asks for `categorical = TRUE`, which
  is frmtmb's default. The user decided this on 2026-09-30 (row :217,
  re-verdicted at consolidation).
- The ordinal warning is given once per call, brms's once per effect.
- The `mi(x, idx = )` refusal names the reason and spares a response
  that does not read the term; brms stops with an index error on every
  response of the model.
- A surface's contour lines default to 10 labeled levels; brms draws 30
  unlabeled ones and a color key, which base `contour()` has no room
  for. `surface_args = list(bins = 30)` gives brms's count.
- Error bars stay arrow-headed bars unless `errorbar_args$width` is
  set; brms's default is `width = 0.3` caps with a dodge of 0.4.
- Observations over a factor predictor keep their small deterministic
  spread (up to 0.08 on either side); brms's `jitter_width` default 0
  draws them stacked. Pre-existing.
- `posterior_samples()` without `pars` keeps `variables()` order, which
  is not brms's intercepts-first order (1.9).
- The default facet layout keeps its roughly square grid,
  `ceiling(sqrt(n))` columns; brms's is `max(floor(sqrt(n)), 3)`.
  Pre-existing and pinned by `test-ce-facets.R`; not this lane's
  question.

## 7. Not done, and defects found but not fixed

- Rows :154, :162, :164, :169, :391 and :396 cannot hold without
  ggplot objects. Their ledger reason should change from "refuses
  brms's arguments" to "returns base-graphics plot objects, not
  ggplot" (a divergence) at consolidation.
- Row :217 does not hold, by the user's decision of 2026-09-30 (a
  deliberate divergence, section 6).
- Row :314 needs `ndraws` on a fit; not done (D3, unrelated).
- `extensions/frmtmb.sample/vignettes/brms-posterior.Rmd` listed
  `posterior_samples()`, `nsamples()` and `parnames()` under "What
  refuses"; they have answered since item 2.6f. The three rows are
  replaced by a sentence saying so.
- `dev/brms-api-diff.md` rows for `nsamples` and `parnames` are stale;
  left for the consolidation, as postfit2 left its rows.

## 8. For the user

- Nothing. Row :217 was decided by the user on 2026-09-30: the
  per-category default stays, as a deliberate divergence from brms.

## 9. Version

Both packages gain exported functions or classes: frmtmb `parnames()`
and the `frmtmb_ce_plot` and `frmtmb_hyp_plot` print and plot methods;
frmtmb.sample re-exports `parnames()`. Both bumps are minor. The
return value of both `plot()` methods changes (BREAKING in NEWS).
frmtmb.sample imports `parnames` from frmtmb, so its frmtmb floor must
rise to the frmtmb version that carries this lane; its DESCRIPTION
still says `frmtmb (>= 0.66.0)`, because no version number is written
in this round.

## 10. Punch round 1

Review `dev/reviews/2026-09-30-ceplot.md`: NOT MERGEABLE on two
blockers. Every item below was checked against the reviewer's own
scripts (`dev/ceplot-rev-*`), rerun on this build into
`dev/ceplot-log/p1-rev-*.txt`, and each fix has a test that a mutant of
it fails (`dev/ceplot-p1-mutants.R`, run by `dev/ceplot-p1-mutrun.R`).

### B1. "old_levels" chooses per grouping factor

`old_level_pick_add()` now keys the choice by the block's grouping
factor and the unseen level's label (`new_level_key()`), and stores the
chosen seen level's LABEL, which `new_level_pick_apply()` (the design,
for `fitted()`) and `predict_new_level_draw()` (the replicate, for
`predict()`) place in each block by matching its levels. Every term of
one grouping factor therefore reads one seen group, as brms's
`get_new_rdraws()`, which runs once per group, gives it. Each design
part now carries `is_new`, the rows that read an unseen level of that
block, so the by-level block a `gr(g, by = f)` row does not read makes
no choice (review m4, the wasted draw), and the candidates are the row's
own by-level's groups, as in brms. The reviewer's
`dev/ceplot-rev-oldlevels.R` and `-oldlevels2.R` on this build:

```
frmtmb from C:/Users/adf44/source/r/wt-ceplot-lib/frmtmb 
(a) frmtmb rows n1,n2,n1 per seed: 9/4/9 5/6/5 5/10/5 8/11/8 2/11/2 5/10/5 10/3/10 4/7/4 11/5/11 11/9/11 10/2/10 2/10/2
(a) brms   rows n1,n2,n1 per seed: 9/4/9 5/6/5 5/10/5 8/11/8 2/11/2 5/10/5 10/3/10 4/7/4 11/5/11 11/9/11 10/2/10 2/10/2
(a) seeds where every row agrees: 12 of 12
(b) frmtmb mu/sigma level per seed: 9/9 5/5 5/5 8/8 2/2 5/5 10/10 4/4 11/11 11/11 10/10 2/2
(b) brms   mu/sigma level per seed: 9/9 5/5 5/5 8/8 2/2 5/5 10/10 4/4 11/11 11/11 10/10 2/2
(b) frmtmb reads ONE seen group for both at 12 of 12 seeds; brms at 12
(b) seed 1 response mean 1.41384; equals a seen group's mean: TRUE
(c) by-level b holds groups 7..12. frmtmb: 7 11 11 9 8 11 8 10 9 9 8 8 | brms: 7 11 11 9 8 11 8 10 9 9 8 8
(c) agree at 12 of 12
(d) (g,h) read. frmtmb: 9,4 5,7 5,2 8,3 2,3 5,2 10,3 4,7 11,5 11,1 10,2 2,2
(d)             brms:   9,4 5,7 5,2 8,3 2,3 5,2 10,3 4,7 11,5 11,1 10,2 2,2
(d) agree at 12 of 12
done
blocks: 2 
seed 1..12: the (x = 0, x = 1) answer equals seen group: 9 5 5 8 2 5 10 4 11 11 10 2 
answers that are one seen group's:  12 of 12
```

`predict(old_levels)` goes through the same key; the reviewer's
no-regression comparison (`dev/ceplot-p1-noreg.R`, a copy writing to
`dev/ceplot-log/`) differs from 0.66.0 in one of 262 results, the
`predict_ol` of `bf(y ~ x + (1 | g), sigma ~ (1 | g))`, which is this
fix:

```
sigma_re   predict_ol    DIFFERS: max diff 1.03
identical 261, differ 1 (errors on base 18)
```

### B2. The rename refuses a row with a grouping variable unset

`ce_plan_part()` calls `ce_plan_refuse()` in the rename branch when any
of the block's variables is `NA` in the row: 0.66.0's refusal of
`y ~ x + trt + (1 | trt:subj)` with nothing set is back, on the
bootstrap and on draws (`dev/ceplot-log/p1-rev-misc.txt`, item 2). With
`trt = "b", subj = "99"` the rename still answers.

### Minors

- m1. `test-ce-levels.R` holds the rename by identity: the display
  minus the same rows predicted with the new `g:h` at no level equals
  the drawn effect, to 64 ulp. `relabel_wrongslot` fails it.
- m2. `test-ce-parity.R` holds a `cs()` variable valid; `valid_nocs`
  fails it.
- m3. `ce_model_vars()` leaves out nonlinear and distributional
  parameter names. `effects = c("a", "x")` drops `a` with brms's
  warning, whose list is `'z', 'x'` (brms: `x z`).
- m4. See B1. The level now agrees with brms at 12 of 12 seeds, where it
  agreed at 3.
- m5. `ps_select()` applies brms's rule itself: every distributional
  parameter's intercept first, a nonlinear parameter's intercept kept
  with its coefficients. brms 2.23.0 measured on four shapes
  (`dev/ceplot-p1-psorder-brms.R`,
  `dev/ceplot-log/p1-psorder-brms.txt`), and the test pins three.
- m6. Both `parnames()` methods warn only when the generic that reached
  them is frmtmb's own; with brms loaded its generic has warned. Tests
  in both packages' `test-generic-collision.R` count the warnings in a
  child process with brms attached (1) and without (1).
- m7. `plot()` names the graphical parameters it ignores in a warning
  (`ce_plot_ignored()`), for both plot methods and the plot methods of
  the plot objects. `frm_check_dots()` gains `.hidden`, so `do_plot` is
  accepted but not listed as an argument. `ask = NULL` is taken as
  `TRUE`, so it is not a breaking change.
- m8. The `method = "predict"` refusal moves into `ce_display_kind()`,
  ahead of the ordinal warning.
- m9. NEWS and `?conditional_effects` say that seed equality with brms
  holds for one new level in a call; the `mm(by = )` ratio is the
  reviewer's 0.83 to 1.11 at 200 refits, the 60-refit explanation is
  withdrawn above, and both say brms 2.23.0 cannot draw an `mm(by = )`
  new level.
- m10. The gated tier below is the whole list of
  `dev/release/run-gated.ps1` (the reviewer's job list,
  `dev/ceplot-p1-jobs-gated.txt`).
- m11. `ce_valid_effects()` reads `ce_model_vars()`, which walks the
  terms itself; a why-comment says why it is kept apart from
  `ce_locked_vars()`.

### Pre-existing defects fixed

- `(1 | g:h)` with integer or character columns. `group_values()`
  evaluates a grouping expression's `:` as the interaction the fit's
  `mkReTrms()` builds (`"<g>:<h>"`), in `pred_design()`,
  `getME("flist")` and the rename label. The reviewer's
  `dev/ceplot-rev-numgrp.R` now gives the integer model's answers equal
  to the factor model's in every row, with no warning
  (`dev/ceplot-log/p1-rev-numgrp.txt`), and `dev/ceplot-rev-misc.R`
  item 1 the same for character columns.
- `poly(x, 2)` and `log(abs(z) + 1)`. The frame keeps the raw variables
  a term reads through a transform (`frame_raw_vars()`, the new frame
  slot `raw_vars`), and the display's grid is built on the model frame
  plus them (`ce_base_frame()`). A variable read only inside `offset()`
  is left out, since lane formrobust owns `offset()` in the display;
  the offset cases of the no-regression run are unchanged (they fail on
  both builds, as before).

### Mutants

Each mutant against the test file that should catch it
(`dev/ceplot-p1-jobs-mut.txt`, logs `dev/ceplot-log/p1mut/`):

```
RESULT lane frmtmb test-ce-levels.R group_values_off: tests=90 failed=0 error=1 skipped=0 warning=1 passed=89
  fails: integer and character grouping columns of g:h predict 
RESULT lane frmtmb test-brms-parity-defects.R ol_offrows: tests=140 failed=6 error=0 skipped=0 warning=0 passed=134
  fails: "old_levels" on gr(by = ) chooses within the row's by-level 
RESULT lane frmtmb test-brms-parity-defects.R ol_perblock: tests=140 failed=13 error=0 skipped=0 warning=0 passed=127
  fails: "old_levels" reads ONE seen group for every term of a factor 
RESULT lane frmtmb test-ce-plot.R plot_ignored_off: tests=68 failed=3 error=0 skipped=0 warning=0 passed=65
  fails: graphical parameters are named when they are ignored 
RESULT lane frmtmb test-ce-parity.R predict_after_warn: tests=30 failed=2 error=0 skipped=0 warning=0 passed=28
  fails: method = "predict" on an ordinal fit stops before it warns 
RESULT lane frmtmb.sample test-brms-shapes-draws.R ps_nl_off: tests=76 failed=2 error=0 skipped=0 warning=0 passed=74
  fails: posterior_samples(pars = ) orders a nonlinear model as brms 
RESULT lane frmtmb test-ce-parity.R rawvars_off: tests=27 failed=0 error=1 skipped=0 warning=0 passed=27
  fails: a variable read through poly() or a transform is an effect 
RESULT lane frmtmb test-ce-levels.R relabel_wrongslot: tests=101 failed=1 error=0 skipped=0 warning=0 passed=100
  fails: a renamed level reads the draw, not a fitted level's effect 
RESULT lane frmtmb test-ce-levels.R rename_na_off: tests=101 failed=1 error=0 skipped=0 warning=0 passed=100
  fails: the rename refuses a row whose grouping variable is unset 
RESULT lane frmtmb test-ce-parity.R valid_nlpars: tests=26 failed=0 error=1 skipped=0 warning=0 passed=26
  fails: a nonlinear parameter's name is not a valid effect 
RESULT lane frmtmb test-ce-parity.R valid_nocs: tests=30 failed=2 error=0 skipped=0 warning=0 passed=28
  fails: a cs() variable is a valid effect 
```

`parnames()`'s single warning has no in-process mutant: its test runs
in a child process, which loads the installed package, not the mutated
namespace. The test fails on the base build (below).

### Seen to fail on the base build

`dev/ceplot-debug1.R <file> <pkg> base`, logs
`dev/ceplot-log/p1-seen-failing-*.txt`. B2, m2 and m8 pass on the base
build, which never had the defect (0.66.0 refused B2 by name, had no
validity check and no ordinal warning); their failing construction is
the mutant above.

```
p1-seen-failing-sample-test-brms-shapes-draws.R.txt:
  ── 1. Error ('test-brms-shapes-draws.R:87:3'): nsamples() is brms's, with brms's
  ── 2. Failure ('test-brms-shapes-draws.R:138:3'): posterior_samples(pars = ) ord
  ── 3. Failure ('test-brms-shapes-draws.R:152:3'): posterior_samples(pars = ) ord
  ── 4. Failure ('test-brms-shapes-draws.R:250:3'): posterior_samples(pars = ) ord
  ── 5. Failure ('test-brms-shapes-draws.R:254:3'): posterior_samples(pars = ) ord
p1-seen-failing-sample-test-generic-collision.R.txt:
  ── 1. Failure ('test-generic-collision.R:62:3'): the package installs the owner 
  ── 2. Failure ('test-generic-collision.R:530:3'): parnames() on draws warns once
p1-seen-failing-test-brms-parity-defects.R.txt:
  ── 1. Error ('test-brms-parity-defects.R:203:5'): fitted() takes sample_new_leve
  ── 2. Error ('test-brms-parity-defects.R:254:5'): "old_levels" reads ONE seen gr
  ── 3. Error ('test-brms-parity-defects.R:288:5'): "old_levels" on gr(by = ) choo
  ── 4. Error ('test-brms-parity-defects.R:303:3'): parnames() on a fit is variabl
p1-seen-failing-test-ce-levels.R.txt:
  ── 1. Error ('test-ce-levels.R:347:5'): an mm() term with a by variable draws ea
  ── 2. Error ('test-ce-levels.R:385:3'): crossed terms at an unseen combination d
  ── 3. Error ('test-ce-levels.R:510:3'): a renamed level reads the draw, not a fi
  ── 4. Error ('test-ce-levels.R:566:5'): integer and character grouping columns o
p1-seen-failing-test-ce-parity.R.txt:
  ── 1. Error ('test-ce-parity.R:25:3'): one invalid effect among valid ones is dr
  ── 2. Failure ('test-ce-parity.R:74:3'): the expected category number of an ordi
  ── 3. Error ('test-ce-parity.R:98:3'): a display of mi(x, idx = ) is refused by 
  ── 4. Error ('test-ce-parity.R:144:3'): a nonlinear parameter's name is not a va
  ── 5. Error ('test-ce-parity.R:181:3'): a variable read through poly() or a tran
p1-seen-failing-test-ce-plot.R.txt:
  ── 1. Failure ('test-ce-plot.R:69:3'): plot() takes every argument brms's two pl
  ── 2. Failure ('test-ce-plot.R:70:3'): plot() takes every argument brms's two pl
  ── 3. Error ('test-ce-plot.R:82:3'): plot = FALSE returns one plot object per ef
  ── 4. Error ('test-ce-plot.R:102:3'): rug draws the observed values of a numeric
  ── 5. Error ('test-ce-plot.R:121:5'): mean = FALSE leaves the estimate out only 
  ── 6. Error ('test-ce-plot.R:133:3'): stype picks contour lines or a raster for 
  ── 7. Error ('test-ce-plot.R:145:3'): layer arguments translate, and the rest ar
  ── 8. Error ('test-ce-plot.R:175:3'): theme, jitter_width and points follow brms
  ── 9. Error ('test-ce-plot.R:200:3'): facet_args lays out the panels of several 
  ── 10. Failure ('test-ce-plot.R:215:3'): plot() of a hypothesis takes brms's arg
  ── 11. Error ('test-ce-plot.R:216:3'): plot() of a hypothesis takes brms's argum
  ── 12. Error ('test-ce-plot.R:240:3'): a long hypothesis is cut in its title as 
  ── 13. Failure ('test-ce-plot.R:252:3'): graphical parameters are named when the
  ── 14. Failure ('test-ce-plot.R:258:3'): graphical parameters are named when the
  ── 15. Error ('test-ce-plot.R:262:3'): graphical parameters are named when they 
p1-seen-failing-test-generic-collision.R.txt:
  ── 1. Failure ('test-generic-collision.R:128:3'): brms loaded but never attached
  ── 2. Failure ('test-generic-collision.R:129:3'): brms loaded but never attached
  ── 3. Failure ('test-generic-collision.R:486:3'): parnames() on a fit warns once
  ── 4. Failure ('test-generic-collision.R:489:3'): parnames() on a fit warns once
  ── 5. Failure ('test-generic-collision.R:490:3'): parnames() on a fit warns once
```

### Tests and checks

All eight suites on the final build, one file per process:

```
## p1full: 329 files, 329 with a RESULT line
lane frmtmb: 195 files, expectations 13827, failed 0, error 0, skipped 163, warning 0, passed 13664
lane frmtmb.coupling: 11 files, expectations 547, failed 0, error 0, skipped 5, warning 0, passed 542
lane frmtmb.eam: 29 files, expectations 1746, failed 0, error 0, skipped 3, warning 0, passed 1743
lane frmtmb.latent: 10 files, expectations 361, failed 0, error 0, skipped 2, warning 0, passed 359
lane frmtmb.learn: 15 files, expectations 442, failed 0, error 0, skipped 13, warning 0, passed 429
lane frmtmb.ode: 11 files, expectations 548, failed 0, error 0, skipped 1, warning 0, passed 547
lane frmtmb.sample: 43 files, expectations 2248, failed 0, error 0, skipped 4, warning 0, passed 2244
lane frmtmb.spline: 15 files, expectations 554, failed 0, error 0, skipped 1, warning 0, passed 553

```

The whole gated list of `dev/release/run-gated.ps1`, with
`FRMTMB_BRMS_FIT_TESTS`, `FRMTMB_DRMTMB_FIT_TESTS` and `FRMTMB_FUZZ`
set:

```
## p1gated: 41 files, 41 with a RESULT line
base frmtmb: 1 files, expectations 162, failed 0, error 0, skipped 0, warning 0, passed 162
base frmtmb.sample: 1 files, expectations 61, failed 0, error 0, skipped 0, warning 0, passed 61
lane frmtmb: 31 files, expectations 2811, failed 3, error 0, skipped 0, warning 0, passed 2808
lane frmtmb.sample: 7 files, expectations 525, failed 2, error 0, skipped 0, warning 0, passed 523
lane frmtmb.learn: 1 files, expectations 71, failed 0, error 0, skipped 0, warning 0, passed 71
not clean:
  RESULT lane frmtmb test-brms-suite-methods.R: tests=162 failed=3 error=0 skipped=0 warning=0 passed=159
    fails: conditional_effects has reasonable ouputs
    fails: variables has reasonable ouputs
  RESULT lane frmtmb.sample test-brms-suite-methods.R: tests=61 failed=2 error=0 skipped=0 warning=0 passed=59
    fails: nsamples has reasonable ouputs
    fails: posterior_samples has reasonable outputs

```

The five gated failures are the five stale "now HOLDS" rows of section
5, confirmed again on this build by `dev/ceplot-flips.R`
(`dev/ceplot-log/flips-lane.txt`, 5 rows: :203, :205, :995, :595,
:635); the gated tier skips nothing.

`R CMD check --as-cran`:

```
frmtmb: Status: 1 NOTE
  * checking HTML version of manual ... [24s] NOTE
  Skipping checking math rendering: package 'V8' unavailable
frmtmb.sample: Status: OK
```

## 11. Punch round 2

Re-check: B2, the pre-existing fixes and m1 to m11 hold; B1 held within
one response and not across two (P1-B1).

### P1-B1. One choice per call, across responses

`fitted()` on a multivariate call now draws the `"old_levels"` choice
once, over the linear predictors of every response it answers
(`fitted_old_levels()` takes a vector of responses), before
`fitted_mv()` hands the fit to each response; a response's own call
reuses a choice already attached. `predict_simulate()` carries one
choice list through `predict_new_level_spec()` for every response
(`pick = `). The key is the grouping factor, the level and, through the
by-level's own block, the by-level. The reviewer's
`dev/ceplot-rev-oldlevels4.R` and `-oldlevels5.R` on this build (brms
2.23.0 at fixed_param, 40 chains of one iteration; data seed 3, calls
at seeds 1 to 12):

```
blocks: y 1 | g y2 1 | g 
fitted  frmtmb y/y2: 9/9 5/5 5/5 8/8 2/2 5/5 10/10 4/4 11/11 11/11 10/10 2/2
fitted  brms   y/y2: 9/9 5/5 5/5 8/8 2/2 5/5 10/10 4/4 11/11 11/11 10/10 2/2
one seen group for both responses: frmtmb 12 of 12, brms 12 of 12
frmtmb's y level equals brms's at 12 of 12
predict_new_level_spec y/y2 per seed (each called at the seed): 9/9 5/5 5/5 8/8 2/2 5/5 10/10 4/4 11/11 11/11 10/10 2/2
predict(mv, old_levels) answers: 5x1x2
done
one predict() call per seed, the choice recorded per response:
  seed 1 : 9 / 9 
  seed 2 : 5 / 5 
  seed 3 : 5 / 5 
  seed 4 : 8 / 8 
  seed 5 : 2 / 2 
  seed 6 : 5 / 5 
  seed 7 : 10 / 10 
  seed 8 : 4 / 4 
  seed 9 : 11 / 11 
  seed 10 : 11 / 11 
  seed 11 : 10 / 10 
  seed 12 : 2 / 2 
```

On the round 1 build the same scripts read 9/4 at seed 1 and one group
for both responses at 0 of 12 seeds (`dev/ceplot-rev-log/p1/`). Tests:
`test-brms-parity-defects.R` holds `fitted()` and `predict()` on the
two-response model to brms's rule, one `sample.int()` per call, at
seeds 1 to 6, and records the choices inside one `predict()` call;
`test-brms-methods.R` (gated) holds `fitted()` against brms at the
same seeds. The divergence entry that section 6 carried for this case
is removed: it was this defect.

### Minors

- `print()` of a plot object takes the graphical parameters `plot()`
  takes and names them in the same warning; it stopped before.
- The hypothesis plot's help says that `col = ` abbreviates `colors = `
  and is that argument, as in brms, so `col = "red"` stops on the
  length of `colors`. A test holds it.
- `?fitted.frmtmb_fit` and `?predict.frmtmb_fit` say that two different
  unseen values of one `mm()` term read two seen levels under
  `"old_levels"`, where brms 2.23.0 gives them one, and that one choice
  serves every response of a multivariate call.

### Tests and check

The core and frmtmb.sample suites on this build, one file per process,
and the gated `test-brms-methods.R`:

```
## p2suite: 238 files, 238 with a RESULT line
lane frmtmb: 195 files, expectations 13863, failed 0, error 0, skipped 164, warning 0, passed 13699
lane frmtmb.sample: 43 files, expectations 2248, failed 0, error 0, skipped 4, warning 0, passed 2244

## p2gated: 1 files, 1 with a RESULT line
lane frmtmb: 1 files, expectations 996, failed 0, error 0, skipped 0, warning 0, passed 996

```

`R CMD check --as-cran`:

```
frmtmb: Status: 1 NOTE
  * checking HTML version of manual ... [20s] NOTE
  Skipping checking math rendering: package 'V8' unavailable
frmtmb.sample: Status: OK
```
