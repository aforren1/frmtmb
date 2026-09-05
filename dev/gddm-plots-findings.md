# gddm vignette figures: working notes

Worktree `frmtmb-wt-gddm-plots`, branch `wt-gddm-plots`, off df9d6bc.
Target: `extensions/frmtmb.ddm/vignettes/gddm.Rmd`, which had eight code
chunks and no figure.

## The four figures

| chunk | section | what it shows |
| --- | --- | --- |
| `fig-bound` | What the generalized model adds | `B(t)` and its mirror for the constant, exponential and linear boundaries, at `bs = 3`, `tau = 1.2`, `kappa = 0.6` over a two-second trial |
| `fig-fit` | Fitting | defective cumulative distributions per coherence, observed steps against the fitted curve, one curve per wall |
| `fig-accuracy` | How accurate, and how expensive | left: solved density against the closed-form Wiener density on both walls, log-log; right: absolute log-density error on the shipped grid and on this vignette's coarser one |
| `fig-grid` | When a row falls off the grid | the fitted upper-wall density near the non-decision time at `dt` 0.04, 0.02 and 0.01 |

Rejected: a signed-response-time histogram with the fitted density laid
over it, which was the first draft of `fig-fit`. The three conditions
differ so much in concentration that one shared density scale squashes
the two easy panels while a per-panel scale stops the row being
comparable, and a bin width fine enough for coherence 0.512 is noise at
coherence 0. Defective cumulative distributions carry the same two
facts, the response times and the choice proportions, on one 0-to-1
scale, and the choice proportion is then the height each curve settles
at.

Also rejected: a figure of tape-build against evaluation cost for
`tridiagonal = "recorded"` and `"atomic"`. It needs two fits, and the
brief was not to fit anything new for a figure.

## Measurements

Chunk times come from a `timeit` knit hook inside one render, so the
figures and the rest of the document are measured under the same machine
load. The machine was shared while this was done and wall-clock renders
varied by nearly a factor of two between runs, which is why the
within-knit share is the number to read.

| render | figures | elapsed |
| --- | --- | --- |
| baseline df9d6bc, three samples | 0 img tags | 37.1, 55.2, 62.3 s |
| with figures, `tinyplot` present | 4 img tags, 4 alt attributes | 49.0, 49.3, 43.0, 49.0 s |
| with figures, `tinyplot` absent | 0 img tags, no errors | 45.0 s |

Within one knit the four figure chunks cost 2.04 to 2.07 s across three
runs, against 40 to 67 s for the whole document: 2.9 to 4.8 percent. The
`fit` chunk is the document. Each figure:

    fig-accuracy  0.74 s   (two solves, dt 0.01/ny 201 and dt 0.02/ny 101)
    fig-fit       0.60 s   (three solves, one per coherence)
    fig-grid      0.53 s   (three solves, one per time step)
    fig-bound     0.20 s   (no solve; three closed-form boundaries)

Every density is solved once and reused. Nothing is fitted for a figure:
`fig-fit` and `fig-grid` read the estimates the `estimates` chunk already
back-transforms.

## Access routes

- Boundary curves are public. `gddm_bound_exponential()$fn(t, p, ctl)`
  returns `list(B, dlogB)`, and `fn` is the seam `gddm_bound_term()`
  documents.
- The solved density is internal: `gd_solve()` and `gd_shift()` behind
  `frmtmb.ddm:::`, which is the idiom the package's own
  `tests/testthat/test-gddm-solver.R` already uses, and
  `ddm_lpdf_both()` for the closed form. There is no public route. The
  family's `lpdf` is reachable only through `frm()`, core frmtmb 0.50.0
  exports no `log_lik`, and `fitted()` and `simulate()` return means and
  draws rather than a density on a grid. The chunk says so in a comment.

## tinyplot 0.7.0, three traps

Each was reproduced on a minimal example before a figure was written
around it.

1. **`log =` plus hand-drawn overlays clip.** After
   `tinyplot(..., type = "n", log = "y")` a following `lines()` is
   clipped to a rectangle much smaller than the panel, so the curve
   simply stops part way across. `par("usr")` is correct and identical
   to what base `plot()` sets, so the coordinates are right and the clip
   region is not. Cure: on a log axis let tinyplot draw every series
   itself through a `by` group. On linear axes `type = "n"` plus
   `lines()` is fine, and `fig-bound` and `fig-fit` use it.
2. **`by` plus tinyplot's own legend destroys `par(mfrow)`.** The first
   panel comes out empty. Cure: `legend = FALSE` and a plain base
   `legend()`. The two then compose, and `abline()` works as well.
3. **`yaxt = "n"` plus a manual `axis(2, ...)` clips the legend.** This
   was isolated with an A/B in one image: identical calls, one with the
   manual axis and one without, and only the manual-axis panel loses the
   right-hand side of its legend. Cure: take tinyplot's own `1e-05` tick
   labels. They are wide, so do not set `par(mar)` either; tinyplot
   overrides it and the axis label then collides with the ticks.

A fourth, smaller: tinyplot drops `NA` rows, so the usual trick of
breaking one series in two with an `NA` does not work. `fig-bound` draws
the mirrored wall as its own `lines()` call instead.

And the theme's title is too large for a third of the page width, so
`fig-fit` names its panels `coh 0.128` rather than `coherence 0.128`;
the full wording runs off the last panel. Titling in the margin with
`mtext()` fits but re-triggers trap 1, clipping the curves at about 0.9
on the vertical axis.

## Prose

No existing claim was contradicted. Two are now pinned by a picture:

- "the two agree to better than 0.01 in the log density" from 0.2 s on.
  Measured on the figure's own case, `mu = 2`, `bs = 2`, `bias = 0.5`:
  the worst log-density error past 0.2 s is 0.0018 on the shipped grid
  and 0.0071 on this vignette's `dt = 0.02`, `ny = 101`. The prose says
  a coarser grid is worse, and it is, by a factor of 4.05.
- "The density at very short decision times ... is much larger than the
  truth." At the first node, 0.01 s, the solver is 1.0e13 times the
  closed form. The error falls below 0.01 at 0.12 s on the shipped grid
  and at 0.16 s on the coarser one.

Numbers used in the new prose, all regenerable: the exponential boundary
closes to `exp(-2 / 1.2)`, about a fifth, over the trial; the lower
wall's fitted and observed plateaus agree to two decimals at 0.50, 0.02
and 0.00; at 0.28 s the `dt = 0.04` density is 2.2 decades above the
`dt = 0.01` one, and the three agree to about a tenth of a decade by
0.36 s.

The sentences around each figure are written to stand without it, since
the gate removes the picture and not the text. The house precedent in
`vignettes/case-studies.Rmd` leaves that prose ungated too.

## Files touched

`extensions/frmtmb.ddm/vignettes/gddm.Rmd`,
`extensions/frmtmb.ddm/DESCRIPTION` (`tinyplot` added to Suggests) and
`extensions/frmtmb.ddm/NEWS.md` (one bullet under the existing 0.2.0
"generalized drift-diffusion model" heading; the repository carries no
tags and `DESCRIPTION` still reads 0.2.0, so that release is the one in
progress and no new heading was opened).
