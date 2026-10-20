# Reviewer of lane sampleprior: the derived sd_/cor_ columns on a
# correlated 3x3 block and a nested group; fast path against core's slow
# path per draw; mutations of the fast path's map that the two-point
# check must reject; the theta_ columns' positions against rellib-r9's.
# Data seed 20261118, sampler seed 21.
source("C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/sampleprior-rev-common.R")
ns <- asNamespace("frmtmb.sample")
core <- asNamespace("frmtmb")
set.seed(20261118)
G <- 12
n <- 240
d <- data.frame(x = rnorm(n), z = rnorm(n),
                g = factor(rep(seq_len(G), each = n / G)),
                h = factor(rep(1:4, length.out = n)))
d$g2 <- factor(rep(1:3, length.out = n))
S <- matrix(c(1, .5, -.3, .5, .8, .2, -.3, .2, .6), 3)
B <- MASS::mvrnorm(G, c(0, 0, 0), S)
d$y <- 1 + B[d$g, 1] + (1 + B[d$g, 2]) * d$x + (0.5 + B[d$g, 3]) * d$z +
  rnorm(4)[d$h] + rnorm(G * 3)[as.integer(interaction(d$g, d$g2))] +
  rnorm(n, 0, 0.7)
fit <- frm(y ~ x + z + (1 + x + z | g) + (1 | h) + (1 | g:g2), data = d)
cat("fit variables():", grep("^sd_|^cor_", variables(fit), value = TRUE), "\n")
# lme4 at the ML estimate: the semantics of the sd and the cor
lm4 <- lme4::lmer(y ~ x + z + (1 + x + z | g) + (1 | h) + (1 | g:g2),
                  data = d, REML = FALSE)
vl <- as.data.frame(lme4::VarCorr(lm4))
vo <- core$hyp_vals_only(fit)
ev <- core$hyp_env_vals(fit, vo$vals, vo$comp)
cat("at the ML estimate, frmtmb (hyp_env_vals) vs lme4:\n")
key <- c(sd_g__Intercept = "g.(Intercept).NA", sd_g__x = "g.x.NA",
         sd_g__z = "g.z.NA", cor_g__Intercept__x = "g.(Intercept).x",
         cor_g__Intercept__z = "g.(Intercept).z", cor_g__x__z = "g.x.z",
         sd_h__Intercept = "h.(Intercept).NA",
         "sd_g:g2__Intercept" = "g:g2.(Intercept).NA")
lk <- paste(vl$grp, vl$var1, vl$var2, sep = ".")
for (k in names(key)) {
  cat(sprintf("  %-22s frmtmb %.5f lme4 %.5f\n", k, ev[[k]],
              vl$sdcor[match(key[[k]], lk)]))
}
ds <- suppressMessages(suppressWarnings(frm_sample(
  fit, chains = 2, iter = 600, refresh = 0, seed = 21)))
m <- ds$draws
dn <- grep("^sd_|^cor_", colnames(m), value = TRUE)
cat("derived columns:", dn, "\n")
# the slow path, per draw, against the matrix
th <- grep("^theta_", colnames(m), value = TRUE)
slow <- t(vapply(seq_len(nrow(m)), function(r) {
  v <- vo$vals
  v[vo$comp == "theta"] <- m[r, th]
  unlist(core$hyp_env_vals(fit, v, vo$comp)[dn])
}, numeric(length(dn))))
cat("max |draws - hyp_env_vals per draw| over", nrow(m), "draws:",
    max(abs(m[, dn] - slow)), "\n")
dv <- ns$draws_natural_cols(fit)$derived
cat("fast path in use:", !is.null(dv$fast), " names on it:",
    length(dv$fast$names), "of", length(dv$names), "\n")
# brms's own VarCorr() on these columns, as brms reads a brmsfit's draws
vcds <- tryCatch(VarCorr(ds), error = function(e) conditionMessage(e))
cat("VarCorr(ds)$g$sd Estimate:", if (is.list(vcds)) vcds$g$sd[, 1] else vcds,
    " | column means", colMeans(m[, c("sd_g__Intercept", "sd_g__x",
                                      "sd_g__z")]), "\n")

# mutations of the fast path's map: each must make the check return NULL
fast_src <- paste(deparse(ns$draws_cov_fast, width.cutoff = 500L), collapse = "
")
mutate <- function(from, to, lab) {
  src <- sub(from, to, fast_src, fixed = TRUE)
  src <- strsplit(src, "
", fixed = TRUE)[[1]]
  if (identical(paste(src, collapse = "
"), fast_src)) { cat(lab, ": pattern not found\n"); return() }
  f <- eval(parse(text = src), envir = ns)
  environment(f) <- ns
  r <- f(fit, dv$names, vo, ev)
  cat(sprintf("mutation %-44s -> %s\n", lab,
              if (is.null(r)) "rejected (slow path)" else
                paste("ACCEPTED with", length(r$names), "names")))
}
mutate('b = bi, i = j, j = j', 'b = bi, i = length(tn) + 1L - j, j = length(tn) + 1L - j',
       "sd names reversed within a block")
mutate('b = bi, i = j, j = k', 'b = bi, i = k, j = j', "cor (i, j) transposed")
mutate('i = j, j = k)', 'i = j, j = if (k == 2L && length(tn) == 3L) 3L else if (k == 3L) 2L else k)',
       "cor_Intercept__x and cor_Intercept__z swapped")
mutate('M[e$i, e$j]/(sqrt(M[e$i, e$i]) * sqrt(M[e$j, e$j]))',
       'M[e$i, e$j]/M[e$i, e$i]', "cor as cov / var_i")
mutate('sqrt(M[e$i, e$i])', 'sqrt(M[e$i, e$i]) * (1 + 1e-6 * (M[e$i, e$i] > 2))',
       "sd off by 1e-6 only where var > 2")
# the unmutated function, for the record
f0 <- ns$draws_cov_fast(fit, dv$names, vo, ev)
cat("unmutated: ", if (is.null(f0)) "NULL" else
  paste("accepted,", length(f0$names), "names"), "\n")
saveRDS(list(cols = colnames(m)), "C:/Users/adf44/source/r/frmtmb-wt-sampleprior/dev/sampleprior-rev-log/names-cols-lane.rds")
