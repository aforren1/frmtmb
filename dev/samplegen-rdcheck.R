# Every Rd this lane changed, verified by RENDERING it, not by reading
# its source: `%` starts an Rd comment even in verbatim macros, and a
# roxygen backtick span can be executed. Each page is rendered with
# Rd2txt and the lines that must survive are grepped for in the output.
#   Rscript dev/samplegen-rdcheck.R
root <- "C:/Users/adf44/source/r/frmtmb-wt-samplegen"
checks <- list(
  "man/frmtmb-sampling-api.Rd" = c(
    "The borrowed-generic seam", "frm_install_generics(pkgname, owners)",
    "S3method(owner::generic, class)", "bare UseMethod()"),
  "extensions/frmtmb.sample/man/frmtmb-loo-refusals.Rd" = c(
    "bridge_sampler(samples, ...)", "bayes_factor(x1, x2, log = FALSE, ...)",
    "post_prob(x, ..., prior_prob = NULL, model_names = NULL)",
    "x, x1, x2, samples, log, prior_prob, model_names, ..."),
  "extensions/frmtmb.sample/man/draws-diagnostics.Rd" = c("rhat(x, ...)"),
  "extensions/frmtmb.sample/man/frmtmb-draws-refusals.Rd" = c(
    "posterior_samples(x, pars = NA, ...)", "object, x, pars, ..."),
  "extensions/frmtmb.sample/man/posterior_epred.Rd" = c(
    "posterior_linpred(object, transform = FALSE, ...)"))
bad <- 0L
for (f in names(checks)) {
  out <- tempfile()
  tools::Rd2txt(file.path(root, f), out = out,
                options = list(underline_titles = FALSE))
  txt <- gsub("[[:space:]]+", " ", paste(readLines(out, warn = FALSE),
                                         collapse = " "))
  # Rd2txt renders code in quotes; compare without them
  txt <- gsub("['\u2018\u2019]", "", txt)
  for (s in checks[[f]]) {
    ok <- grepl(s, txt, fixed = TRUE)
    if (!ok) bad <- bad + 1L
    cat(sprintf("%-4s %-45s %s\n", if (ok) "ok" else "MISS", basename(f), s))
  }
}
cat("missing:", bad, "\n")
