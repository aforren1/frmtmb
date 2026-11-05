#!/usr/bin/env bash
# L1 identities through the equiv instrument for variants of sampled
# pairs: one edit file per variant, written here, then
# dev/equiv-variant.R measures each on this lane's library. Logs in
# dev/refuse-variant-out/.
export PATH="$PATH:/c/rtools45/usr/bin:/c/rtools45/x86_64-w64-mingw32.static.posix/bin"
export FRMTMB_EQUIV_LIB="C:/Users/adf44/source/r/wt-refuse-lib;C:/Users/adf44/source/r/rellib-r16"
export NOT_CRAN=true
mkdir -p dev/refuse-variant-out
cd dev/refuse-variant-out
cat > sds_t2.edit.R <<'EOF'
spec$bf <- "bf(y ~ t2(x, z))"
spec$prior <- "set_prior('student_t(3, 0, 1)', class = 'sds')"
EOF
cat > sds_t2coef.edit.R <<'EOF'
spec$bf <- "bf(y ~ t2(x, z))"
spec$prior <- "set_prior('student_t(3, 0, 1)', class = 'sds', coef = 't2(x, z)')"
EOF
cat > sds_by.edit.R <<'EOF'
spec$bf <- "bf(y ~ f + s(x, by = f))"
spec$prior <- "set_prior('student_t(3, 0, 1)', class = 'sds', coef = 's(x, by = f)')"
EOF
cat > sds_two.edit.R <<'EOF'
spec$bf <- "bf(y ~ s(x) + s(z, k = 5))"
spec$prior <- "set_prior('student_t(3, 0, 1)', class = 'sds', coef = 's(z, k = 5)')"
EOF
cat > sds_sigma.edit.R <<'EOF'
spec$bf <- "bf(y ~ s(x), sigma ~ s(z))"
spec$prior <- "set_prior('student_t(3, 0, 1)', class = 'sds', dpar = 'sigma')"
EOF
cat > sds_beside_sd.edit.R <<'EOF'
spec$bf <- "bf(y ~ s(x) + (1 | g))"
spec$prior <- "set_prior('student_t(3, 0, 1)', class = 'sds') + set_prior('normal(0, 2)', class = 'sd')"
EOF
cat > skew_se.edit.R <<'EOF'
spec$family <- "skew_normal()"
spec$prior <- NULL
d <- spec$data
d$s <- 0.2 + 0.3 * abs(d$z)
spec$data <- d
spec$bf <- "bf(y | se(s) ~ x)"
EOF
cat > skew_se_sigma.edit.R <<'EOF'
spec$family <- "skew_normal()"
spec$prior <- NULL
d <- spec$data
d$s <- 0.2 + 0.3 * abs(d$z)
spec$data <- d
spec$bf <- "bf(y | se(s, sigma = TRUE) ~ x)"
EOF
cat > mi_se.edit.R <<'EOF'
spec$prior <- NULL
d <- spec$data
d$s <- 0.2 + 0.3 * abs(d$z)
d$y[seq(5, 195, by = 19)] <- NA
spec$data <- d
spec$bf <- "bf(y | mi() + se(s) ~ x)"
EOF
cat > mi_sd_se.edit.R <<'EOF'
spec$prior <- NULL
d <- spec$data
d$s <- 0.2 + 0.3 * abs(d$z)
d$s2 <- 0.2 + 0.3 * abs(d$x)
spec$data <- d
spec$bf <- "bf(y | mi(s2) + se(s) ~ x)"
EOF
cat > id_diag.edit.R <<'EOF'
spec$prior <- NULL
spec$bf <- "bf(y ~ x + (1 + z | gr(g, id = 'p', cor = FALSE)), sigma ~ (1 | gr(g, id = 'p', cor = FALSE)))"
EOF
cat > mi_gev.edit.R <<'EOF'
spec$prior <- NULL
set.seed(5)
d <- spec$data
d$y <- 0.5 + 0.4 * d$x + 0.6 * (rexp(nrow(d))^(-0.1) - 1) / 0.1
d$y[seq(5, 195, by = 13)] <- NA
spec$data <- d
spec$family <- "gen_extreme_value()"
spec$bf <- "bf(y | mi() ~ x, xi ~ 1)"
EOF
cat > mi_lognormal_sd.edit.R <<'EOF'
spec$prior <- NULL
set.seed(6)
d <- spec$data
d$y <- rlnorm(nrow(d), 0.5 + 0.4 * d$x, 0.6)
d$s <- 0.1
spec$data <- d
spec$family <- "lognormal()"
spec$bf <- "bf(y | mi(s) ~ x)"
EOF
cd ../..
for v in sds_t2 sds_t2coef sds_by sds_two sds_sigma sds_beside_sd skew_se skew_se_sigma mi_se mi_sd_se id_diag mi_gev mi_lognormal_sd; do echo $v; done | \
  xargs -P 7 -I{} bash -c '"/c/Program Files/R/R-4.6.1/bin/Rscript.exe" dev/equiv-variant.R "s() x prior" dev/refuse-variant-out/{}.edit.R > dev/refuse-variant-out/{}.txt 2>&1'
