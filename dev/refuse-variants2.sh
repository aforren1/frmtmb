#!/usr/bin/env bash
# Second batch of variants, as dev/refuse-variants.sh.
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export FRMTMB_EQUIV_LIB="C:/Users/adf44/source/r/wt-refuse-lib;C:/Users/adf44/source/r/rellib-r16"
export NOT_CRAN=true
cd dev/refuse-variant-out
cat > sds_two2.edit.R <<'E2'
spec$bf <- "bf(y ~ s(x) + s(z, k = 5))"
spec$prior <- "set_prior('student_t(3, 0, 1)', class = 'sds') + set_prior('normal(0, 2)', class = 'sds', coef = 's(z, k = 5)')"
E2
cat > sds_sigma2.edit.R <<'E2'
spec$bf <- "bf(y ~ x, sigma ~ s(z))"
spec$prior <- "set_prior('student_t(3, 0, 1)', class = 'sds', dpar = 'sigma')"
E2
cat > flat_sigma_smooth.edit.R <<'E2'
spec$bf <- "bf(y ~ x, sigma ~ s(z))"
spec$prior <- NULL
E2
cat > cens_mix.edit.R <<'E2'
spec$prior <- NULL
set.seed(7)
d <- spec$data
k <- rbinom(nrow(d), 1, 0.4)
d$y <- ifelse(k == 1, 3 + 0.5 * d$x, -3 - 0.5 * d$x) + rnorm(nrow(d))
d$ce <- 0L
d$ce[d$y > 3.5] <- 1L
d$ce[d$y < -3.5] <- -1L
spec$data <- d
spec$family <- "mixture(gaussian(), gaussian())"
spec$bf <- "bf(y | cens(ce) ~ 1, mu1 ~ x, mu2 ~ x)"
E2
cd ../..
for v in sds_two2 sds_sigma2 flat_sigma_smooth cens_mix; do echo $v; done | \
  xargs -P 4 -I{} bash -c '"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" dev/equiv-variant.R "s() x prior" dev/refuse-variant-out/{}.edit.R > dev/refuse-variant-out/{}.txt 2>&1'
