# Review of lane wcheck2, item 4968: the point-mass families near the
# boundary of the rule. For each arm: the check's verdict (frm()), and
# the Laplace objective's profile along the spike dpar's intercept (the
# rest re-optimized), slope per unit toward the limit. Predicted slope:
# rate * E_nz (E_nz: sum of w - 1 over the rows off the point mass),
# less 1 per unit for a prior's Jacobian where one applies.
#   Rscript dev/reviews/wcheck2-review-mass.R <base|lane> <fam>
# fam: hl, hg, zib, zoib, zial
a <- commandArgs(TRUE)
lib_arg <- a[1]
fam_key <- a[2]
source("dev/wcheck2-common.R")
ns <- asNamespace("frmtmb")

set.seed(20261113)
n <- 80
d <- data.frame(obs = factor(seq_len(n)), x = rnorm(n),
                g4 = factor(rep(1:20, each = 4)))
d$y <- 1 + 0.5 * d$x + rnorm(n)
i0 <- seq_len(n) %% 5 == 0
i1 <- seq_len(n) %% 7 == 0 & !i0
d$yz <- ifelse(i0, 0, exp(d$y / 2))
d$yb <- plogis(d$y - 1); d$yb[i0] <- 0
d$yo <- d$yb; d$yo[i1] <- 1
d$ya <- ifelse(i0, 0, d$y)
spec <- switch(fam_key,
  hl = list(fam = hurdle_lognormal(), y = "yz", sp = "sigma", zero = TRUE,
            rate = 1, mass = i0, zi = "hu",
            alt = hurdle_lognormal(link_sigma = "identity")),
  hg = list(fam = hurdle_gamma(), y = "yz", sp = "shape", zero = FALSE,
            rate = 0.5, mass = i0, zi = "hu",
            alt = hurdle_gamma(link_shape = "identity")),
  zib = list(fam = zero_inflated_beta(), y = "yb", sp = "phi", zero = FALSE,
             rate = 0.5, mass = i0, zi = "zi",
             alt = zero_inflated_beta(link_phi = "identity")),
  zoib = list(fam = zero_one_inflated_beta(), y = "yo", sp = "phi",
              zero = FALSE, rate = 0.5, mass = i0 | i1, zi = "zoi",
              alt = zero_one_inflated_beta(link_phi = "identity")),
  zial = list(fam = zero_inflated_asym_laplace(), y = "ya", sp = "sigma",
              zero = TRUE, rate = 1, mass = i0, zi = "zi", alt = NULL))
pm <- spec$mass
nnz <- sum(!pm)
cat(sprintf("%s: %d mass rows, %d others\n", fam_key, sum(pm), nnz))

# index in betad of the spike dpar's intercept, from the frame
spike_idx <- function(fo, fam, w) {
  d$w <- w
  fr <- q(frm(fo, family = fam, data = d, dry_run = "frame"))
  k <- grep(paste0("(^|_|\\|)", spec$sp, "$"), names(fr$linpreds))
  if (!length(k)) k <- grep(spec$sp, names(fr$linpreds))
  fr$linpreds[[k[1]]]$idx[1]
}
arm <- function(lab, w, fo, fam = spec$fam, prior = NULL, jac = 0,
                prof = TRUE) {
  d$w <- w
  E <- sum(w[!pm] - 1)
  lin <- fam_key == "zial"
  pred <- (if (lin) sum(w[!pm]) else spec$rate * E) - jac
  cat(sprintf("--- %s: E_nz %.4g, sum(w) - n %.4g, predicted slope %.4g\n",
              lab, E, sum(w) - n, pred))
  cat("  frm():", outcome(frm(fo, family = fam, data = d, prior = prior)),
      "\n")
  if (!prof) return(invisible())
  ob <- q(frm(fo, family = fam, data = d, prior = prior,
              dry_run = "objective"))
  j <- spike_idx(fo, fam, w)
  p0 <- (if (!is.null(ob$obj)) ob$obj else ob)$par
  jj <- match(j, which(names(p0) == "betad"))
  if (is.na(jj)) jj <- 1L
  g <- c(2, 4, 6, 8, 10, 12)
  gl <- if (spec$zero) -g else g
  v <- profile_disp(ob, gl, j = jj)
  show_profile("  profile, toward the limit", g, v)
}
fo <- as.formula(paste(spec$y, "| weights(w) ~ x + (1 | obs)"))
# near the boundary, both sides
arm("E +0.6, mass 0.5 (sum(w) < n)", ifelse(pm, 0.5, 1 + 0.6 / nnz), fo)
arm("E -0.6, mass 3 (sum(w) > n)", ifelse(pm, 3, 1 - 0.6 / nnz), fo)
# mixed weights off the mass, rescaled to E = +0.4 and -0.4
set.seed(5)
u <- runif(n, 0.3, 1.7)
for (E in c(0.4, -0.4)) {
  w <- u
  w[!pm] <- u[!pm] * (nnz + E) / sum(u[!pm])
  w[pm] <- runif(sum(pm), 0, 4)
  arm(sprintf("mixed weights, E %+.1f, mass random in (0, 4)", E), w, fo)
}
# the mass rows at weight 0: they leave pos and need no effect
w <- ifelse(pm, 0, 1 + 0.6 / nnz)
arm("E +0.6, mass weight 0", w, fo, prof = FALSE)
if (fam_key != "zial") {
  # a formula on the mass probability
  fz <- bf(fo, as.formula(paste(spec$zi, "~ x")))
  arm(sprintf("E +0.6, %s ~ x", spec$zi), ifelse(pm, 0.5, 1 + 0.6 / nnz),
      fz, prof = FALSE)
  arm(sprintf("E -0.6, %s ~ x", spec$zi), ifelse(pm, 3, 1 - 0.6 / nnz), fz,
      prof = FALSE)
  # identity link on the spike dpar
  arm("E +0.6, identity link", ifelse(pm, 0.5, 1 + 0.6 / nnz), fo,
      fam = spec$alt, prof = FALSE)
  arm("E -0.6, identity link", ifelse(pm, 3, 1 - 0.6 / nnz), fo,
      fam = spec$alt, prof = FALSE)
}
# a prior on the spike dpar (no formula: natural scale, log link's
# Jacobian 1 per unit on the MAP route)
pr <- switch(spec$sp,
  sigma = set_prior("normal(0, 1)", class = "sigma"),
  shape = set_prior("gamma(0.01, 0.01)", class = "shape"),
  phi = set_prior("gamma(0.01, 0.01)", class = "phi"))
for (E in c(0.6, 1.6, 3)) {
  arm(sprintf("E %+.1f, prior %s on %s", E, pr[[1]]$prior %||% "", spec$sp),
      ifelse(pm, 0.5, 1 + E / nnz), fo, prior = pr, jac = 1,
      prof = spec$sp == "sigma")
}
# a grouped effect, not observation-level: never refused
fg <- as.formula(paste(spec$y, "| weights(w) ~ x + (1 | g4)"))
arm("E +0.6, (1 | g4)", ifelse(pm, 0.5, 1 + 0.6 / nnz), fg, prof = FALSE)
cat("DONE\n")
