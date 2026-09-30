# Lane ceplot punch 1: NEWS corrections (B1, m4, m7, m9) and the two
# pre-existing fixes.
f <- "C:/Users/adf44/source/r/frmtmb-wt-ceplot/NEWS.md"
s <- paste(readLines(f), collapse = "\n")
rep1 <- function(s, old, new) {
  n <- lengths(regmatches(s, gregexpr(old, s, fixed = TRUE)))
  if (n != 1L) stop("expected one match, got ", n, ": ", substr(old, 1, 60))
  sub(old, new, s, fixed = TRUE)
}
s <- rep1(s, "  image with contours.
* **`conditional_effects(effects = )`",
"  image with contours. Graphical parameters (`col`, `main`, ...) are
  ignored as before, and now with a warning that names them.
* **`predict(sample_new_levels = \"old_levels\")` chooses one seen level
  per grouping factor**, as brms's `get_new_rdraws()` does, not one per
  term. On `bf(y ~ x + (1 | g), sigma ~ (1 | g))` or `(1 + x || g)` a
  new group took each term's effects from a different seen group, a
  group the data does not have. Under `gr(g, by = f)` the choice is made
  among the row's own by-level only, with no draw for the by-level the
  row does not read.
* **`conditional_effects(effects = )`")
s <- rep1(s, "* `fitted()` takes `sample_new_levels = \"old_levels\"`. Each unseen
  level reads one seen level of its block, chosen at random once per
  call as brms chooses it and as `predict()` chooses it: the row takes
  that level's fitted effect, and `Est.Error` that level's conditional
  variance. At seeds 1 to 12, brms 2.23.0 and frmtmb choose the same
  level (`dev/ceplot-oldlevels.R`).",
"* `fitted()` takes `sample_new_levels = \"old_levels\"`. Each unseen
  level reads one seen level of its grouping factor, chosen at random
  once per call as brms chooses it and as `predict()` chooses it, and
  every term of the factor reads that level: the row takes its fitted
  effects, and `Est.Error` their conditional variance. At seeds 1 to 12
  brms 2.23.0 and frmtmb choose the same level on `(1 | g)` with one or
  two unseen levels, on `(1 | g) + (1 | h)`, on `(1 | g)` in `mu` and
  `sigma`, and within the row's by-level of `gr(g, by = f)`
  (`dev/ceplot-oldlevels.R`, `dev/ceplot-rev-oldlevels.R`).")
s <- rep1(s, "  draws. The bootstrap band has the Wald band's width (1.02 to 1.05 at
  60 refits), and on draws the band equals brms's at the same draws and
  seed to 4.4e-16 (`dev/ceplot-crossed-brms.R`). A new member of an
  `mm()` term with a `by` variable is drawn in its own by-level's block:
  the bootstrap band is 0.81 to 0.95 of the Wald band's width at 60
  refits, about the shrinkage of a standard deviation estimated from
  five groups (`dev/ceplot-mmby.R`).",
"  draws. The bootstrap band has the Wald band's width (0.89 to 1.08 at
  200 refits, seeds 5 to 7). With one new level in the call, the band
  on draws equals brms's at the same draws and seed to 4.4e-16
  (`dev/ceplot-crossed-brms.R`); with two, brms draws level by level
  and frmtmb draw by draw, so the bands agree in law and not number. A
  new member of an `mm()` term with a `by` variable is drawn in its own
  by-level's block: the bootstrap band is 0.83 to 1.11 of the Wald
  band's width at 200 refits. brms 2.23.0 cannot draw such a level at
  all; it stops with one of three errors.")
s <- rep1(s, "  the prediction. brms 2.23.0 stops on the same grid.",
"  the prediction. brms 2.23.0 stops on the same grid.
* A grouping term `(1 | g:h)` on integer or character columns predicts
  on new data. `g:h` was evaluated as R's sequence operator, so every
  `predict()`, `fitted()` and `conditional_effects()` on new data
  stopped with \"non-conformable arrays\" (and a \"numerical expression
  has 3 elements\" warning escaped), or reported observed combinations as
  new levels. It is now the interaction the fit reads.
* `conditional_effects()` varies a variable the model reads through a
  transform, `x` in `poly(x, 2)` or `z` in `log(abs(z) + 1)`, as brms
  does. The model frame keeps only the transformed columns, so naming
  such a variable stopped with \"not stored in the model frame\" and the
  default display found nothing to draw. The fit keeps the raw
  variables beside the model frame for it.
* `conditional_effects()` refuses a bootstrap or draws band by name
  again on a row whose grouping variable is unset where the only way to
  place its new level would be to rename one (`y ~ x + trt + (1 |
  trt:subj)` with nothing set). The first build of the renamed level
  returned `NA` there without a message.")
writeLines(strsplit(s, "\n", fixed = TRUE)[[1]], f)
cat("edited\n")
