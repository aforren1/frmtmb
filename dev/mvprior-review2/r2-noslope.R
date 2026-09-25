# Item 3 over-reach: class "b" where the only "slope" is a special term.
who <- Sys.getenv("R2_WHO", "lane")
Sys.setenv(R2_ARM = if (who == "base") "base" else "lane")
source("C:/Users/adf44/source/r/frmtmb-wt-mvprior/dev/mvprior-review2/r2-prelude.R")
if (who == "brms") suppressPackageStartupMessages(library(brms))
ns <- asNamespace(if (who == "brms") "brms" else "frmtmb")
source(file.path(R2_ROOT, "dev/mvprior-review2/r2-cases.R"))
d <- r2_data()
d$ordf <- d$ord
mods <- list(
  smooth = y ~ s(x),
  mono = yg ~ mo(ordf),
  ranslope = yg ~ 1 + (x | g),
  offset = yg ~ 1 + offset(x),
  noint = yg ~ 0 + x,
  gp = yg ~ gp(x),
  cat_by = yg ~ 1 + (1 | g) + s(z),
  cs = ord ~ cs(x))
names(d)[names(d) == "yg"] <- "yg"; d$y <- d$yg
for (nm in names(mods)) {
  fam <- if (nm == "cs") ns$acat() else gaussian()
  r <- tryCatch({
    v <- as.data.frame(ns$validate_prior(ns$set_prior("normal(0, 1)", class = "b"), mods[[nm]], data = d, family = fam))
    v[is.na(v)] <- ""
    u <- v[v$source == "user", ]
    paste("ACCEPT:", paste(u$class, u$coef, sep = "|", collapse = " ; "))
  }, error = function(e) paste("REFUSE:", substr(gsub("[[:space:]]+", " ", conditionMessage(e)), 1, 220)))
  cat(sprintf("%-9s %-26s %s\n", nm, deparse(mods[[nm]]), r))
}
