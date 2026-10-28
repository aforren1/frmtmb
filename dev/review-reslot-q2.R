# Reviewer short check, item 2. Usage: Rscript dev/review-reslot-q2.R
# brms|before|after. Seeds: data 1, points 99, sampler 1.
#  A. brms: is default_prior() passed back to brm() brms's defaults?
#  B. frm() under default_prior() (fit route) and validate_prior()
#     output equals frm() without / with the plain spec, bitwise.
#  C. frm_sample() under frm()'s route table, unmodified and with one
#     row edited (the brms porting workflow): which densities?
#  D. (after) round trip of a flat run's and a default run's
#     prior_summary(), object and table, on models the lane did not use.
arm <- commandArgs(trailingOnly = TRUE)[1L]
mkdata <- function() {
  set.seed(1)
  G <- 10L; n <- 8L
  d <- data.frame(g = factor(rep(seq_len(G), each = n)), x = rnorm(G * n),
                  z = rnorm(G * n), w = rnorm(G * n))
  d$y <- 1 + 0.5 * d$x + rnorm(G, 0, 0.7)[d$g] + rnorm(G * n)
  d$o <- sample(1:4, G * n, TRUE)
  d$ym <- ifelse(runif(G * n) < 0.4, rnorm(G * n, -2), rnorm(G * n, 2))
  d
}
d <- mkdata()
f1 <- y ~ x + (1 + x | g)
if (identical(arm, "brms")) {
  .libPaths(c("C:/Users/adf44/source/r/rellib-r14",
              "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
  suppressPackageStartupMessages(library(brms))
  tab <- default_prior(bf(f1), data = d)
  print(tab[, c("prior", "class", "coef", "group")])
  s0 <- stancode(bf(f1), data = d)
  s1 <- stancode(bf(f1), data = d, prior = tab)
  cat("A. stancode with default_prior() passed back identical:",
      identical(s0, s1), "\n")
  tab2 <- tab
  tab2$prior[tab2$class == "b" & tab2$coef == "x"] <- "normal(0, 1)"
  s2 <- stancode(bf(f1), data = d, prior = tab2)
  cat(grep("lprior \\+=", strsplit(s2, "\n")[[1]], value = TRUE), sep = "\n")
  # a class-level "" row: (flat)
  s3 <- stancode(bf(f1), data = d, prior = set_prior("", class = "sd"))
  cat("class-level '' on sd leaves an sd term:",
      any(grepl("lprior.*sd_1", strsplit(s3, "\n")[[1]])), "\n")
  quit(save = "no")
}
source("C:/Users/adf44/source/r/frmtmb-wt-reslot/dev/reslot-common.R")
suppressPackageStartupMessages(library(frmtmb.sample))
q <- function(expr) suppressWarnings(suppressMessages(expr))
ents <- function(ri) {
  if (is.null(ri)) return("none")
  paste(vapply(ri$entries, function(e) {
    paste0(e$comp, "[", paste(e$idx, collapse = ","), "]:", e$dist$kind,
           if (isTRUE(e$dist$flat)) "-flat" else "")
  }, ""), collapse = " ")
}

cat("== B. frm() fits under commonly passed tables\n")
same <- function(a, b) identical(a$opt$par, b$opt$par) &&
  identical(a$opt$objective, b$opt$objective)
tabf <- default_prior(f1, data = d)
print(as.data.frame(tabf)[, c("prior", "class", "coef", "group")],
      row.names = FALSE)
f_ml <- q(frm(f1, data = d))
f_tab <- q(frm(f1, data = d, prior = tabf))
cat("  frm(prior = default_prior(fit route)) == frm():", same(f_ml, f_tab),
    " fit$prior NULL:", is.null(f_tab$prior), "\n")
sp <- set_prior("normal(0, 1)", class = "b")
vp <- validate_prior(sp, f1, data = d)
f_sp <- q(frm(f1, data = d, prior = sp))
f_vp <- q(frm(f1, data = d, prior = vp))
cat("  frm(prior = validate_prior(sp)) == frm(prior = sp):", same(f_sp, f_vp),
    "\n")
tabe <- tabf
tabe$prior[tabe$class == "b" & tabe$coef == "x"] <- "normal(0, 1)"
f_te <- q(frm(f1, data = d, prior = tabe))
f_bx <- q(frm(f1, data = d, prior = set_prior("normal(0, 1)", class = "b",
                                               coef = "x")))
cat("  frm(prior = edited table) == frm(prior = b x spec):", same(f_te, f_bx),
    "\n")
cat("  objectives:", format(c(f_ml$opt$objective, f_tab$opt$objective,
                              f_sp$opt$objective, f_vp$opt$objective,
                              f_te$opt$objective, f_bx$opt$objective),
                            digits = 15), "\n")

cat("\n== C. frm_sample() under frm()'s route table\n")
fitd <- q(frm(f1, data = d, dry_run = "objective"))
rs <- function(pr) {
  r <- tryCatch(q(frmtmb.sample:::sample_resolve_priors(fitd, pr)),
                error = function(e) e)
  if (inherits(r, "error")) paste("ERROR", conditionMessage(r)) else ents(r$ri)
}
cat("  defaults:               ", rs(NULL), "\n")
cat("  flat:                   ", rs("flat"), "\n")
cat("  fit-route table:        ", rs(tabf), "\n")
cat("  fit-route table, edited:", rs(tabe), "\n")
tabs <- default_prior(f1, data = d, route = "sample")
cat("  sample-route table:     ", rs(tabs), "\n")
tabse <- tabs
tabse$prior[tabse$class == "b" & tabse$coef == "x"] <- "normal(0, 1)"
cat("  sample table, edited:   ", rs(tabse), "\n")

if (identical(arm, "after")) {
  cat("\n== D. round trips of prior_summary()\n")
  cmp <- function(fit, rp0, rp1, label) {
    fit0 <- fit; fit0$prior <- rp0$effective
    fit1 <- fit; fit1$prior <- rp1$effective
    o0 <- frmtmb.sample:::prior_augmented_obj(fit0, rp0$ri$entries %||% list())
    o1 <- frmtmb.sample:::prior_augmented_obj(fit1, rp1$ri$entries %||% list())
    set.seed(99)
    df <- dg <- 0
    for (j in 1:40) {
      p <- o0$par + stats::rnorm(length(o0$par), 0, 0.7)
      f0 <- o0$fn(p); g0 <- o0$gr(p)
      df <- max(df, abs(o1$fn(p) - f0) / abs(f0))
      dg <- max(dg, max(abs(o1$gr(p) - g0)) / max(abs(g0)))
    }
    cat(sprintf("  %-46s rel fn %.3g  rel grad %.3g\n", label, df, dg))
  }
  models <- list(
    lkj4 = list(f = y ~ x + (1 + x + z + w | g), fam = gaussian()),
    mix = list(f = ym ~ x, fam = mixture(gaussian(), gaussian())),
    mo = list(f = y ~ mo(o) + (1 | g), fam = gaussian()))
  for (nm in names(models)) {
    m <- models[[nm]]
    fit <- q(frm(m$f, data = d, family = m$fam))
    for (pr in list("flat", NULL)) {
      lab <- paste(nm, if (is.null(pr)) "default" else "flat")
      ds <- q(frm_sample(fit, prior = pr, chains = 1, iter = 60, warmup = 30,
                         seed = 1, refresh = 0))
      ps <- prior_summary(ds)
      rp0 <- q(frmtmb.sample:::sample_resolve_priors(fit, pr))
      rp1 <- q(frmtmb.sample:::sample_resolve_priors(fit, ps))
      cmp(fit, rp0, rp1, paste(lab, "object"))
      tb <- as.brmsprior(as.data.frame(ps))
      rp2 <- q(frmtmb.sample:::sample_resolve_priors(fit, tb))
      cmp(fit, rp0, rp2, paste(lab, "as.brmsprior(table)"))
      cat("    mix_theta:", isTRUE(rp0$mix_theta), isTRUE(rp1$mix_theta),
          isTRUE(rp2$mix_theta), " entries:", ents(rp0$ri), "|",
          ents(rp2$ri), "\n")
    }
  }
}
