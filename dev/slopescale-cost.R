# Lane slopescale, item 1018, brief step 5: what autoscale's random-slope
# plan costs.
# (a) Every fit now plans the random-effect columns (autoscale_plan_z())
#     and looks for a slope it cannot rescale (autoscale_note_unscaled()),
#     whether or not the default engages. Timed directly on the frame of
#     each design, against the whole fit, minimum of 5 rounds of blocks
#     grown past 1.2 s.
# (b) Where the default engages (a slope column spread 1e-4), the fit
#     pays for the pre-fit: timed against autoscale = FALSE on the same
#     data (which stalls) and against the fit at spread 1.
# Usage: Rscript slopescale-cost.R
.libPaths(c("C:/Users/adf44/source/r/wt-slopescale-lib",
            "C:/Users/adf44/source/r/rellib-r21",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages(library(frmtmb))
cat("frmtmb", as.character(packageVersion("frmtmb")),
    find.package("frmtmb"), "\n")
ns <- asNamespace("frmtmb")
blk_time <- function(f, min_block = 1.2, rounds = 5) {
  n <- 1L
  repeat {
    t <- system.time(for (i in seq_len(n)) f())[["elapsed"]]
    if (t >= min_block) break
    n <- n * 2L
  }
  min(vapply(seq_len(rounds), function(r) {
    system.time(for (i in seq_len(n)) f())[["elapsed"]] / n
  }, 0))
}
designs <- list(
  sleepstudy = function() {
    d <- lme4::sleepstudy
    list(fo = Reaction ~ Days + (Days | Subject), fam = gaussian(), d = d)
  },
  pois5k = function() {
    set.seed(1)
    g <- factor(rep(1:100, each = 50))
    x <- rnorm(5000)
    y <- rpois(5000, exp(0.3 + rnorm(100, 0, .4)[g] +
                          rnorm(100, 0, .3)[g] * x))
    list(fo = y ~ x + (1 + x | g), fam = poisson(),
         d = data.frame(y, g, x))
  },
  gau1e5 = function() {
    set.seed(2)
    g <- factor(rep(1:1000, each = 100))
    x <- rnorm(1e5)
    z <- rnorm(1e5)
    y <- 1 + rnorm(1000, 0, .5)[g] + rnorm(1000, 0, .3)[g] * x + 0.2 * z +
      rnorm(1e5)
    list(fo = y ~ x + z + (1 + x + z | g), fam = gaussian(),
         d = data.frame(y, g, x, z))
  })
ctl <- frmtmb_control()
cat("\n(a) per-fit planning cost, design: fit s, plan_z s, note s, ratio\n")
for (nm in names(designs)) {
  ds <- designs[[nm]]()
  f <- suppressWarnings(frm(ds$fo, family = ds$fam, data = ds$d))
  # a block past 1.2 s, as the planning steps; the 1e5-row fit once
  fit_t <- if (nm == "gau1e5") {
    system.time(suppressWarnings(frm(ds$fo, family = ds$fam,
                                     data = ds$d)))[["elapsed"]]
  } else {
    blk_time(function() frm(ds$fo, family = ds$fam, data = ds$d),
             rounds = 3)
  }
  fr <- f$frame
  tz <- blk_time(function() ns$autoscale_plan_z(fr))
  tn <- blk_time(function() ns$autoscale_note_unscaled(fr, NULL, ctl))
  cat(sprintf("  %-10s fit %.4f  plan_z %.2e  note %.2e  (plan_z+note)/fit %.2e\n",
              nm, fit_t, tz, tn, (tz + tn) / fit_t))
}
cat("\n(b) engaged default: s 1e-4 default vs autoscale = FALSE vs spread 1\n")
for (nm in c("sleepstudy", "pois5k")) {
  ds <- designs[[nm]]()
  col <- if (nm == "sleepstudy") "Days" else "x"
  fo <- if (nm == "sleepstudy") Reaction ~ 1 + (Days | Subject) else
    y ~ 1 + (1 + x | g)
  d1 <- ds$d
  dsm <- d1
  dsm[[col]] <- d1[[col]] * 1e-4
  q <- function(...) suppressWarnings(suppressMessages(frm(...)))
  # each a block past 1.2 s, minimum of 5 rounds
  t_def <- blk_time(function() q(fo, family = ds$fam, data = dsm))
  t_off <- blk_time(function() q(fo, family = ds$fam, data = dsm,
                                 control = frmtmb_control(autoscale = FALSE)))
  t_one <- blk_time(function() q(fo, family = ds$fam, data = d1))
  cat(sprintf("  %-10s default %.4f  autoscale=FALSE %.4f  spread 1 %.4f  default/spread1 %.2f\n",
              nm, t_def, t_off, t_one, t_def / t_one))
}
