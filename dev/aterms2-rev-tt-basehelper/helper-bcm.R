# The Bayesian Cognitive Modeling port (test-bcm-*.R and
# vignette("bayesian-cognitive-modeling")).
#
# The families the port needs and the core does not ship live in
# inst/bcm/, so the vignette and these tests fit exactly the same code.
# A test helper could not be the shared source: a vignette has no one
# relative path that reaches tests/testthat/ under both `R CMD build`
# (working directory vignettes/) and pkgdown (working directory the
# package root). inst/ is the one directory an installed package and a
# source build resolve the same way. This is helper-rl.R's arrangement,
# for the same reason.

bcm_sources <- c("binomial-extras.R", "process-trees.R", "similarity.R",
                 "marginal.R")

bcm_source_paths <- vapply(bcm_sources, function(f) {
  system.file("bcm", f, package = "frmtmb")
}, "")

for (p in bcm_source_paths) {
  if (nzchar(p)) source(p, local = TRUE)
}

skip_unless_bcm <- function(which = bcm_sources) {
  ok <- all(nzchar(bcm_source_paths[which]))
  testthat::skip_if_not(ok,
                        paste0("inst/bcm/", paste(which, collapse = ", "),
                               " is not installed"))
}

# The Pearson correlation program of chapter 5, adapted. Two chapters
# use it: Correlation_1 is the model, and ESP's Ability and
# OptionalStopping are the same model on different columns. A test file
# cannot read another test file's definitions, so it lives here.
bcm_corr1_code <- function() {
  paste(
    "data {",
    "  int<lower=0> n;",
    "  array[n] vector[2] x;",
    "}",
    "parameters {",
    "  vector[2] mu;",
    "  vector<lower=0>[2] sigma;",
    "  real<lower=-1, upper=1> r;",
    "}",
    "model {",
    "  matrix[2, 2] T;",
    "  T[1, 1] = square(sigma[1]);",
    "  T[1, 2] = r * sigma[1] * sigma[2];",
    "  T[2, 1] = T[1, 2];",
    "  T[2, 2] = square(sigma[2]);",
    "  target += multi_normal_lpdf(x | mu, T);",
    "}",
    sep = "\n")
}

# The data of Lee and Wagenmakers chapter 12, the psychophysical
# functions of Kuijpers and colleagues as the case study distributes
# them: eight subjects, unequal numbers of stimulus levels, and the
# subject-centered stimulus value the model regresses on.
#
# One long data frame rather than the padded matrices the Stan program
# takes, because that is the shape a formula reads. The padding is what
# the -99 sentinel in the original script exists for.
bcm_psy_data <- function() {
  nstim <- c(27L, 24L, 27L, 22L, 25L, 26L, 28L, 28L)
  xmean <- c(318.888, 311.0417, 284.4444, 301.5909,
             296.2000, 305.7692, 294.6429, 280.3571)
  x <- c(
    200, 220, 240, 260, 270, 275, 280, 285, 290, 295, 300, 305, 310, 330,
    335, 340, 345, 350, 355, 360, 365, 370, 375, 380, 385, 390, 400,
    200, 220, 240, 260, 270, 280, 285, 290, 295, 300, 305, 310, 315, 320,
    325, 330, 340, 350, 355, 360, 365, 370, 380, 400,
    200, 220, 225, 230, 235, 240, 245, 250, 255, 260, 265, 270, 275, 280,
    285, 290, 295, 300, 305, 310, 315, 320, 330, 340, 360, 380, 400,
    200, 220, 240, 260, 270, 275, 280, 285, 290, 295, 300, 305, 310, 315,
    320, 325, 330, 335, 340, 360, 380, 400,
    200, 220, 240, 250, 255, 260, 265, 270, 275, 280, 285, 290, 295, 300,
    305, 310, 315, 320, 325, 330, 335, 340, 360, 380, 400,
    200, 220, 240, 260, 265, 270, 275, 280, 285, 290, 295, 300, 305, 310,
    315, 320, 325, 330, 335, 340, 345, 350, 355, 360, 380, 400,
    180, 200, 210, 220, 230, 235, 240, 245, 250, 255, 260, 265, 270, 305,
    310, 315, 320, 325, 330, 335, 340, 345, 350, 355, 360, 380, 400, 420,
    200, 210, 215, 220, 225, 230, 235, 240, 245, 250, 260, 265, 270, 275,
    280, 285, 290, 295, 300, 305, 310, 315, 320, 330, 340, 360, 380, 400)
  r <- c(
    0, 0, 0, 0, 0, 0, 2, 4, 4, 3, 1, 2, 3, 0, 3, 17, 9, 13, 4, 6, 4, 6,
    5, 8, 2, 2, 8,
    0, 0, 0, 0, 0, 0, 3, 1, 2, 3, 0, 2, 6, 10, 7, 3, 7, 4, 6, 11, 13, 12,
    6, 8,
    0, 2, 1, 3, 2, 8, 2, 4, 2, 3, 2, 1, 3, 11, 9, 8, 6, 8, 10, 14, 8, 6,
    6, 6, 6, 6, 6,
    0, 0, 0, 0, 0, 4, 5, 2, 2, 4, 3, 4, 1, 4, 22, 17, 19, 4, 6, 6, 6, 6,
    0, 0, 0, 0, 1, 4, 1, 1, 3, 5, 4, 3, 3, 4, 4, 8, 10, 19, 9, 16, 2, 6,
    6, 6, 6,
    0, 0, 0, 0, 3, 3, 2, 3, 5, 5, 2, 3, 8, 7, 10, 11, 2, 7, 2, 17, 4, 4,
    2, 6, 6, 6,
    0, 1, 0, 1, 0, 1, 4, 2, 10, 5, 6, 1, 4, 1, 11, 4, 9, 5, 13, 1, 12, 7,
    4, 2, 6, 6, 7, 2,
    0, 0, 0, 3, 1, 6, 4, 11, 2, 4, 2, 2, 4, 5, 9, 7, 13, 8, 8, 8, 6, 6,
    6, 6, 6, 6, 6, 6)
  n <- c(
    6, 6, 6, 6, 6, 4, 16, 19, 17, 9, 9, 11, 5, 5, 6, 24, 13, 16, 4, 9, 6,
    8, 6, 11, 2, 2, 8,
    6, 6, 6, 6, 6, 12, 12, 11, 13, 10, 6, 10, 15, 19, 11, 5, 9, 7, 8, 19,
    16, 12, 7, 8,
    8, 15, 4, 9, 12, 20, 7, 11, 5, 7, 3, 3, 10, 16, 11, 9, 7, 10, 14, 15,
    8, 6, 6, 6, 6, 6, 6,
    6, 6, 6, 6, 13, 13, 15, 11, 12, 13, 12, 6, 3, 12, 30, 24, 21, 4, 6,
    6, 6, 9,
    6, 6, 6, 8, 3, 14, 4, 13, 10, 17, 14, 9, 9, 7, 8, 10, 16, 23, 14, 17,
    2, 6, 6, 6, 6,
    6, 6, 6, 12, 9, 13, 7, 18, 22, 10, 5, 8, 12, 12, 13, 12, 3, 12, 4,
    20, 5, 5, 2, 6, 6, 6,
    2, 9, 2, 7, 5, 4, 18, 10, 26, 8, 20, 2, 7, 5, 12, 6, 13, 10, 16, 4,
    16, 9, 5, 2, 6, 6, 8, 2,
    8, 2, 2, 19, 7, 31, 15, 25, 3, 5, 3, 2, 4, 7, 9, 10, 14, 9, 9, 8, 6,
    6, 6, 6, 6, 6, 6, 6)
  subj <- rep(seq_along(nstim), nstim)
  data.frame(subj = factor(subj), x = x, r = r, n = n,
             xc = x - xmean[subj])
}
