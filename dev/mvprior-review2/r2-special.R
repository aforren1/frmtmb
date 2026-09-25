# What class "b" reaches on special-term predictors, lane arm.
source("C:/Users/adf44/source/r/frmtmb-wt-mvprior/dev/mvprior-review2/r2-prelude.R")
source(file.path(R2_ROOT, "dev/mvprior-review2/r2-cases.R"))
d <- r2_data(); d$ordf <- d$ord; d$y <- d$yg
ns <- asNamespace("frmtmb")
for (f in list(yg ~ mo(ordf), y ~ s(x), yg ~ 0 + x, yg ~ z + mo(ordf))) {
  des <- ns$prior_design(f, d, gaussian(), list())
  pt <- des$frame[["par_template"]]
  r <- tryCatch({
    rr <- ns$resolve_priorlist(des, set_prior("normal(0, 0.123)", class = "b"))
    paste(vapply(rr$entries, function(e) paste(ns$par_template_names(pt[[e$comp]], e$comp)[e$idx], collapse = "+"), ""), collapse = " ")
  }, error = function(e) paste("REFUSE:", conditionMessage(e)))
  cat(deparse(f), "->", if (nzchar(r)) r else "(no density)", "\n   pars:",
      unlist(lapply(names(pt), function(c) ns$par_template_names(pt[[c]], c))), "\n")
}
