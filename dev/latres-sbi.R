# brms 2.23.0's table of which family and link has a built-in Stan density
# on the link scale (`sbi_<link>` specials) and which have a joint link:
# stan_eta_transform() reads these two to decide whether `mu` passes through
# the inverse link before the ARMA loop. dev/latres-findings.md quotes it.
.libPaths(c("C:/Users/adf44/source/r/rellib-r16",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
fams <- ls(asNamespace("brms"), pattern = "^\\.family_", all.names = TRUE)
for (f in fams) {
  nm <- sub("^\\.family_", "", f)
  i <- brms:::family_info(nm, "specials")
  li <- brms:::family_info(nm, "links")
  cat(sprintf("%-30s links=%s specials=%s\n", nm, paste(li, collapse = ","),
              paste(i, collapse = ",")))
}
