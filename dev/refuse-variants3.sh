#!/usr/bin/env bash
# Third batch: a censored mixture against its uncensored control.
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export FRMTMB_EQUIV_LIB="C:/Users/adf44/source/r/wt-refuse-lib;C:/Users/adf44/source/r/rellib-r16"
export NOT_CRAN=true
cd dev/refuse-variant-out
cat > mix_control.edit.R <<'E2'
spec$prior <- NULL
set.seed(7)
d <- spec$data
k <- rbinom(nrow(d), 1, 0.4)
d$y <- ifelse(k == 1, 3 + 0.5 * d$x, -3 - 0.5 * d$x) + rnorm(nrow(d))
spec$data <- d
spec$family <- "mixture(gaussian(), gaussian(), order = 'none')"
spec$bf <- "bf(y ~ 1, mu1 ~ x, mu2 ~ x)"
E2
cat > cens_mix_none.edit.R <<'E2'
spec$prior <- NULL
set.seed(7)
d <- spec$data
k <- rbinom(nrow(d), 1, 0.4)
d$y <- ifelse(k == 1, 3 + 0.5 * d$x, -3 - 0.5 * d$x) + rnorm(nrow(d))
d$ce <- 0L
d$ce[d$y > 3.5] <- 1L
d$ce[d$y < -3.5] <- -1L
spec$data <- d
spec$family <- "mixture(gaussian(), gaussian(), order = 'none')"
spec$bf <- "bf(y | cens(ce) ~ 1, mu1 ~ x, mu2 ~ x)"
E2
cd ../..
for v in cens_mix_none; do echo $v; done | \
  xargs -P 2 -I{} bash -c '"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" dev/equiv-variant.R "s() x prior" dev/refuse-variant-out/{}.edit.R > dev/refuse-variant-out/{}.txt 2>&1'
