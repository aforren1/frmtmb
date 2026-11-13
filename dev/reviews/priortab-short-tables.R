# Short check of lane priortab's fix round: per-row coerced-table check.
# (a) untouched sampling tables, coerced and not, over many models: read,
#     with the same resolved entries as the defaults alone;
# (b) mixed frames.
# PRIORTAB_LIB unset = lane library; PRIORTAB_LIB=none = rellib-r20.
lib <- Sys.getenv("PRIORTAB_LIB", "C:/Users/adf44/source/r/wt-priortab-lib")
if (identical(lib, "none")) lib <- ""
.libPaths(c(if (nzchar(lib)) lib, "C:/Users/adf44/source/r/rellib-r20",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressPackageStartupMessages({library(frmtmb); library(frmtmb.sample)})
cat("frmtmb", dirname(find.package("frmtmb")), "\n")
options(width = 150, digits = 10)
set.seed(1)
n <- 120
d <- data.frame(x = rnorm(n), z = rnorm(n), g = factor(rep(1:10, each = 12)))
d$y <- 1 + d$z * 0.3 + rnorm(10)[d$g] * 0.3 + rnorm(n, sd = 0.3)
d$y2 <- d$y + rnorm(n)
d$Y_a <- d$y; d$Y_b <- d$y2
d$yp <- abs(d$y) + 0.1
d$cnt <- rpois(n, 3)
d$o <- cut(d$y, 4, labels = FALSE)
d$om <- factor(sample(1:4, n, TRUE), ordered = TRUE)
d$t <- rep(1:12, 10)
cases <- list(
  re = list(y ~ z + (1 + z | g), gaussian()),
  student = list(y ~ z + (1 | g), student()),
  sigma_x = list(bf(y ~ z, sigma ~ x + (1 | g)), gaussian()),
  mv = list(bf(mvbind(Y_a, Y_b) ~ z + (1 | g)) + set_rescor(TRUE), gaussian()),
  nl = list(bf(yp ~ a * exp(b * z), a ~ 1 + (1 | g), b ~ 1, nl = TRUE),
            gaussian()),
  cumul = list(o ~ z + (1 | g), cumulative()),
  smooth = list(y ~ s(x, k = 5) + (1 | g), gaussian()),
  mo = list(y ~ mo(om), gaussian()),
  zip = list(bf(cnt ~ z, zi ~ x), zero_inflated_poisson()),
  ar = list(y ~ z + ar(t, g), gaussian()),
  mix = list(y ~ z, mixture(gaussian, gaussian)),
  shifted = list(yp ~ z, shifted_lognormal()),
  negbin = list(cnt ~ z + (1 | g), negbinomial())
)
ent <- function(rp) {
  if (is.null(rp$ri)) return("none")
  vapply(rp$ri$entries, function(e) paste(e$comp, paste(e$idx, collapse = ","),
                                          format_prior_dist_safe(e$dist)), "")
}
format_prior_dist_safe <- function(dd) {
  tryCatch(paste(deparse(dd, width.cutoff = 500L), collapse = ""),
           error = function(e) "?")
}
rs <- function(fit, pr) {
  tryCatch(suppressMessages(frmtmb.sample:::sample_resolve_priors(fit, pr)),
           error = function(e) structure(conditionMessage(e), class = "err"))
}
cat("\n## (a) untouched sampling tables\n")
for (nm in names(cases)) {
  cs <- cases[[nm]]
  fit <- tryCatch(suppressWarnings(frm(cs[[1]], data = d, family = cs[[2]],
                                       dry_run = "objective")),
                  error = function(e) NULL)
  if (is.null(fit)) { cat(nm, ": frm dry run failed\n"); next }
  tab <- tryCatch(default_prior(cs[[1]], data = d, family = cs[[2]],
                                route = "sample"),
                  error = function(e) NULL)
  if (is.null(tab)) { cat(nm, ": no sample table\n"); next }
  r0 <- rs(fit, NULL)
  r1 <- rs(fit, frmtmb:::as_priorlist(tab))
  r2 <- rs(fit, as.brmsprior(as.data.frame(tab)))
  st <- function(r) if (inherits(r, "err")) paste("REFUSED:", substr(r, 1, 160))
                    else if (identical(ent(r), ent(r0))) "same entries"
                    else "DIFFERENT entries"
  cat(sprintf("%-8s table: %s | coerced: %s\n", nm, st(r1), st(r2)))
  if (!inherits(r2, "err") && !identical(ent(r2), ent(r0))) {
    cat("   defaults: ", paste(setdiff(ent(r0), ent(r2)), collapse = " ; "), "\n")
    cat("   coerced : ", paste(setdiff(ent(r2), ent(r0)), collapse = " ; "), "\n")
  }
}
cat("\n## (b) mixed frames, y ~ z + x + (1 | g)\n")
f <- y ~ z + x + (1 | g)
fit <- frm(f, data = d, dry_run = "objective")
ps <- as.data.frame(default_prior(f, data = d, route = "sample"))
pf <- as.data.frame(default_prior(f, data = d))
mix <- function(label, tab) {
  r <- rs(fit, as.brmsprior(tab))
  cat(sprintf("%-60s %s\n", label, if (inherits(r, "err"))
    paste("REFUSED:", sub(".*on the rows ([^.]*)\\..*", "\\1", r)) else
    paste("read;", paste(ent(r), collapse = " ; "))))
}
m1 <- ps; m1$prior[m1$class == "sigma"] <- "normal(0, 1)"
m1$prior[m1$class == "b" & m1$coef == "z"] <- "normal(0, 2)"
m1$prior[m1$class == "Intercept"] <- "(flat)"
mix("M1 sample: sigma, b z densities; Intercept (flat)", m1)
m2 <- m1; m2$prior[m2$class == "Intercept"] <- ""
mix("M2 same, Intercept \"\"", m2)
m3 <- pf; m3$prior[m3$class == "sigma"] <- "student_t(3, 0, 2.5)"
m3$prior[m3$class == "Intercept"] <- "normal(1, 2)"
mix("M3 fit table: sigma, Intercept densities, sd (flat)", m3)
m4 <- m3; m4$prior[m4$class == "sd" & m4$group == ""] <- "student_t(3, 0, 2.5)"
mix("M4 fit table: every default slot given a density", m4)
m5 <- m1; m5$prior[m5$class == "Intercept"] <- "student_t(3, 1, 2.5)"
m5$prior[m5$class == "sd" & m5$group == "g"] <- "(flat)"
mix("M5 sample: sd g row (flat) (brms leaves it blank)", m5)
m6 <- ps; m6$prior[m6$class == "sd" & m6$group == ""] <- "(flat)"
m6$prior[m6$class == "sd" & m6$group == "g" & m6$coef == ""] <- "normal(0, 1)"
mix("M6 sample: sd class (flat), sd g a density", m6)
m7 <- ps; m7$source <- "user"
m7$prior[m7$class == "Intercept"] <- "(flat)"
mix("M7 sample relabelled user, Intercept (flat)", m7)
mix("M8 hand frame: (flat) on b only",
    data.frame(prior = "(flat)", class = "b", coef = "z"))
mix("M9 hand frame: (flat) on sigma only",
    data.frame(prior = "(flat)", class = "sigma"))
cat("\ndone\n")
