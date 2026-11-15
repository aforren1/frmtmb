# The release review's F3: brmsfit-methods:169 plots the spaghetti display
# that line 167 builds. Line 167 is refused on a maximum-likelihood fit
# (row 167, cannot transfer), and hygiene's helper fix (backlog 3042) now
# marks line 169 as reading that stale name, so the row no longer
# measures a plot. It is re-adjudicated as row 167 is; the ggplot
# divergence stays measured at rows 162 and 164.
f <- "dev/brmsport-verdicts-manual-fit.tsv"
x <- readLines(f)
i <- grep("^brmsfit-methods:169\t", x)
stopifnot(length(i) == 1)
x[i] <- paste("brmsfit-methods:169", "frmtmb", "cannot transfer", "no-draws",
  paste("plot() of the spaghetti display of brmsfit-methods:167, whose call",
        "is refused on a maximum-likelihood fit (row 167: one curve per draw",
        "needs draws, and its Wald band is refused by name), so the line",
        "reads a stale name and plots nothing; the base-graphics plot",
        "objects against is(object, 'ggplot') are measured at",
        "brmsfit-methods:162 and :164"), sep = "\t")
con <- file(f, "wb"); writeLines(x, con, sep = "\r\n"); close(con)
cat("row 169 re-adjudicated\n")
