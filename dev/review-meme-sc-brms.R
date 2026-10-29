# Review of lane meme, short check: brms 2.23.0's default_prior() on the
# short-check shapes, validate_prior() for me(x) with no sdx, and the
# coef-bound refusal for b, sd, meanme, sdme. No compile.
#   Rscript dev/review-meme-sc-brms.R > dev/review-meme-out/sc/brms.txt
.libPaths("C:/Users/adf44/AppData/Local/R/win-library/4.6")
suppressPackageStartupMessages(library(brms))
wt <- "C:/Users/adf44/source/r/frmtmb-wt-meme"
source(file.path(wt, "dev/review-meme-shapes.R"))
source(file.path(wt, "dev/review-meme-sc-shapes.R"))
d <- sc_data()
sh <- sc_shapes(bf, mvbf, set_rescor, cumulative, gp)
out <- list()
for (nm in names(sh)) {
  s <- sh[[nm]]
  t <- tryCatch(as.data.frame(default_prior(s$f, data = d,
                                            family = s$fam %||% gaussian())),
                error = function(e) conditionMessage(e))
  cat("====", nm, "\n")
  print(t)
  out[[nm]] <- t
}
out$nosdx_validate <- as.data.frame(validate_prior(
  set_prior("normal(0, 1)", class = "meanme"), bf(y ~ me(x) + z), data = d))
print(out$nosdx_validate)
cat("==== coef + bound, stancode\n")
f <- bf(y ~ x + w + (1 + x | g) + me(xg, sg, gr = g) + me(wg, sg, gr = g))
ps <- list(
  b_coef_ub = set_prior("normal(0, 1)", class = "b", coef = "w", ub = 0.1),
  sd_coef_ub = set_prior("normal(0, 1)", class = "sd", group = "g",
                         coef = "x", ub = 2),
  sd_coef_lb0 = set_prior("normal(0, 1)", class = "sd", group = "g",
                          coef = "x", lb = 0),
  sd_coef_lb0_ub2 = set_prior("normal(0, 1)", class = "sd", group = "g",
                              coef = "x", lb = 0, ub = 2),
  sd_group_ub = set_prior("normal(0, 1)", class = "sd", group = "g", ub = 2),
  meanme_coef_lb = set_prior("normal(0, 1)", class = "meanme", coef = "mexg",
                             lb = 0),
  sdme_coef_ub = set_prior("exponential(1)", class = "sdme", coef = "mexg",
                           ub = 0.5),
  sdme_coef_lb0 = set_prior("exponential(1)", class = "sdme", coef = "mexg",
                            lb = 0),
  sdme_class_ub = set_prior("exponential(1)", class = "sdme", ub = 3))
for (k in names(ps)) {
  r <- tryCatch({stancode(f, data = d, prior = ps[[k]]); "accepted"},
                error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(sprintf("%-16s %s\n", k, r))
}
fg <- bf(y ~ gp(u) + x)
tg <- default_prior(fg, data = d)
r <- tryCatch({stancode(fg, data = d, prior = tg); "accepted"},
              error = function(e) paste("ERROR:", conditionMessage(e)))
cat(sprintf("%-16s %s\n", "gp table back", r))
saveRDS(out, file.path(wt, "dev/review-meme-out/sc/brms.rds"))
cat("DONE\n")
