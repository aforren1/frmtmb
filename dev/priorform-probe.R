# Constructions for the priorform lane's items, run against either the
# reference build of the base commit or the lane's own library, so the
# before and after come from one script.
#   Rscript dev/priorform-probe.R ref
#   Rscript dev/priorform-probe.R lane
args <- commandArgs(trailingOnly = TRUE)
which_lib <- if (length(args)) args[1] else "ref"
lib <- switch(which_lib,
  ref = "C:/Users/adf44/source/r/rellib-r3",
  lane = "C:/Users/adf44/source/r/priorform-lib")
.libPaths(c(lib, "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(frmtmb))
cat("frmtmb", as.character(packageVersion("frmtmb")), "from",
    find.package("frmtmb"), "\n\n")

show <- function(label, expr) {
  cat("---", label, "---\n")
  out <- tryCatch(expr, error = function(e) paste("ERROR:",
                                                   conditionMessage(e)))
  print(out)
  cat("\n")
}

set.seed(20260916)
d <- data.frame(y = rnorm(60), x = rnorm(60), z = rnorm(60),
                g = factor(rep(1:6, 10)))

cat("### item 4: bf() idempotence\n")
show("identical(form, bf(form)), dpar", {
  form <- bf(y ~ a, sigma ~ 1)
  identical(form, bf(form))
})
show("identical(form, bf(form)), nl", {
  form <- bf(y ~ a, sigma ~ 1, a ~ x, nl = TRUE)
  identical(form, bf(form))
})

cat("### item 6: set_prior() vector arguments\n")
show("set_prior('normal(0, 2)', class = c('b', 'sd'))", {
  bp <- set_prior("normal(0, 2)", class = c("b", "sd"))
  list(n = length(bp), prior = bp$prior, class = bp$class)
})

cat("### item 7: the rest of the prior surface\n")
for (nm in c("default_prior", "validate_prior", "empty_prior",
             "as.brmsprior")) {
  cat(sprintf("  %-16s exported: %s\n", nm,
              nm %in% getNamespaceExports("frmtmb")))
}
show("print(set_prior('normal(0,1)'))",
     utils::capture.output(print(set_prior("normal(0,1)"))))

cat("### item 12: duplicated group-level effects\n")
show("(1 | g) + (x | g)", {
  frm(y ~ x + (1 | g) + (x | g), d, dry_run = "frame")
  "accepted"
})

cat("### item 13: x * cs(g)\n")
show("y ~ x * cs(g), cumulative", {
  dd <- d
  dd$yo <- factor(cut(dd$y, 3), ordered = TRUE)
  frm(yo ~ x * cs(g), dd, family = cumulative(), dry_run = "frame")
  "accepted"
})

cat("### item 14: bf(y ~ ~ x)\n")
show("bf(y ~ ~ x)", {
  b <- bf(y ~ ~x)
  "accepted"
})
