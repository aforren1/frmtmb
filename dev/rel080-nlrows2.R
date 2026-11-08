# Which basis the fit's tape carries on the three rows of rel080-nlrows.R,
# and whether the probe wins once its runs step without the body's block.
#
#   REV_ARM=<lib> Rscript dev/rel080-nlrows2.R
src <- readLines("dev/review-nlstart-fresh2.R")
stop_at <- grep("^cat\\(\"ROWHEAD", src)
eval(parse(text = src[1:(stop_at - 1L)]))
rows <- list(c("gaussian", "gomp", 22), c("poisson", "two_sat_rev", 33),
             c("bernoulli", "mm", 11))
for (rw in rows) {
  fam <- rw[1]; body <- rw[2]; seed <- as.numeric(rw[3])
  d <- gen(fam, body, seed * 1000 + match(body, names(B)))
  fm <- form(fam, body)
  r <- withr::with_options(list(frmtmb.nl_sign_probe = FALSE), suppressWarnings(
    suppressMessages(frm(fm, family = F[[fam]], data = d))))
  ob <- r$obj$env$.frmtmb_opt_basis
  cat(fam, body, seed, "| basis:", if (is.null(ob)) "none" else
        paste(vapply(ob$blocks, `[[`, "", "kind"), collapse = ","),
      "nl", isTRUE(ob$nl), "| ll", format(as.numeric(logLik(r)), digits = 9), "\n")
}
cat("NLROWS2 DONE\n")
