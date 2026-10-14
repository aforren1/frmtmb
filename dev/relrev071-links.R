# Release reviewer of 0.71.0: every row of the generated brms_mu_links
# table against brms 2.23.0's .family_<name>()$links, and the
# brms_no_mixture list.
.libPaths(c("C:/Users/adf44/source/r/rellib-r9",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
cat("brms", as.character(packageVersion("brms")), "frmtmb",
    as.character(packageVersion("frmtmb")), "\n")
ns <- asNamespace("frmtmb")
bns <- asNamespace("brms")
tab <- get("brms_mu_links", ns)
src <- get("brms_mu_link_source", ns)
bad <- 0
for (f in names(src)) {
  b <- tryCatch(get(paste0(".family_", src[[f]]), bns)()$links,
                error = function(e) NULL)
  if (is.null(b)) { cat("no brms family", f, "\n"); next }
  if (!identical(b, tab[[f]])) {
    bad <- bad + 1
    cat("DIFF", f, ":", tab[[f]], "| brms", b, "\n")
  }
}
cat("families", length(src), "pairs", length(unlist(tab)), "mismatches",
    bad, "\n")
for (f in c("dirichlet", "dirichlet2", "logistic_normal", "frechet",
            "gen_extreme_value", "discrete_weibull")) {
  cat(f, ": brms", get(paste0(".family_", f), bns)()$links,
      "| frmtmb", tab[[f]], "\n")
}
# the links each new frmtmb family accepts for mu, from its constructor
for (f in c("frechet", "gen_extreme_value", "discrete_weibull",
            "logistic_normal", "dirichlet")) {
  fam <- get(f, asNamespace("frmtmb"))
  for (l in tab[[f]]) {
    r <- tryCatch({ fam(link = l); "ok" }, error = function(e)
      paste("ERR", substr(conditionMessage(e), 1, 70)))
    if (r != "ok") cat("  ", f, l, r, "\n")
  }
}
b <- brms::brmsfamily("dirichlet2")
cat("brmsfamily('dirichlet2') link", b$link, "\n")
r <- tryCatch(frmtmb::brmsfamily("dirichlet2", link = "softplus"),
              error = function(e) conditionMessage(e))
cat("frmtmb brmsfamily('dirichlet2', link = 'softplus'):",
    if (is.character(r)) substr(r, 1, 120) else r$link, "\n")
