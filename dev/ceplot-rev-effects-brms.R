# Reviewer check (lane ceplot): brms 2.23.0's valid effects for the model
# shapes of dev/ceplot-rev-effects.R, read from
# get_all_effects(brmsterms(), comb_all = TRUE) without a fit: the
# single variables brms accepts as an effect.
#   Rscript dev/ceplot-rev-effects-brms.R > dev/ceplot-rev-log/effects-brms.txt
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressMessages(library(brms))
fm <- list(
  plain = bf(y ~ x + z),
  inter = bf(y ~ x * f),
  smooth = bf(y ~ s(x) + z),
  smooth_by = bf(y ~ s(x, by = f) + f),
  t2 = bf(y ~ t2(x, z)),
  gp = bf(y ~ gp(x) + z),
  gp_by = bf(y ~ gp(x, by = f) + f),
  mo = bf(y ~ mo(o) + x),
  me = bf(y ~ me(xm, sdx) + z),
  mi = bf(y ~ mi(xmi) + w) + bf(xmi | mi() ~ w),
  nl = bf(yp ~ a * exp(b * x), a ~ 1 + z, b ~ 1, nl = TRUE),
  dpar = bf(y ~ x, sigma ~ z),
  reslope = bf(y ~ x + (1 + w | g)),
  offset = bf(yc ~ x + offset(off)),
  weights = bf(y | weights(wt) ~ x),
  cs = bf(yo ~ cs(x) + z, family = acat()),
  poly = bf(y ~ poly(x, 2) + log(abs(z) + 1)),
  grby = bf(y ~ x + (1 | gr(g, by = fg))),
  mm = bf(y ~ x + (1 | mm(g1, g2))),
  gh = bf(y ~ x + (1 | g:h)),
  mix = bf(y ~ 1, mu1 ~ x, mu2 ~ z, family = mixture(gaussian, gaussian)),
  zi = bf(yc ~ x, zi ~ z, family = zero_inflated_poisson()),
  ar = bf(y ~ x + ar(time = t, gr = g)),
  fs = bf(y ~ x + s(z, g, bs = "fs", k = 4)),
  mv = bf(y ~ x) + bf(y2 ~ z),
  lf = bf(y ~ x) + lf(sigma ~ w)
)
for (nm in names(fm)) {
  f <- fm[[nm]]
  r <- tryCatch({
    bt <- brms:::brmsterms(f)
    rsv <- brms:::rsv_vars(bt)
    ae <- brms:::get_all_effects(bt, rsv_vars = rsv, comb_all = TRUE)
    sort(unlist(ae[lengths(ae) == 1]))
  }, error = function(e) paste("ERROR", conditionMessage(e)))
  cat(sprintf("%-10s %s\n", nm, paste(r, collapse = " ")))
}
