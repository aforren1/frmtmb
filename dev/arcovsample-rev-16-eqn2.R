# REVIEW re-check, detail: the two-argument \eqn in man/frmtmb-autocor.Rd
# through Rd2HTML and Rd2latex as well as Rd2txt.
#
#   Rscript dev/arcovsample-rev-16-eqn2.R

f <- "C:/Users/adf44/source/r/frmtmb-wt-arcovsample/man/frmtmb-autocor.Rd"

h <- capture.output(tools::Rd2HTML(f, out = stdout(), package = "frmtmb"))
k <- grep("max", h, fixed = TRUE)
cat("Rd2HTML: ", length(h), " lines; lines mentioning 'max':\n", sep = "")
for (i in k) cat("  [", i, "] ", trimws(h[i]), "\n", sep = "")

l <- capture.output(tools::Rd2latex(f, out = stdout()))
cat("\nRd2latex: ", length(l), " lines parsed\n", sep = "")
k2 <- grep("max", l, fixed = TRUE)
for (i in k2) cat("  [", i, "] ", trimws(l[i]), "\n", sep = "")

cat("\nRd parse warnings, if any:\n")
w <- tryCatch({ tools::checkRd(f); "checkRd returned no condition" },
              warning = function(cnd) paste("WARNING:",
                                            conditionMessage(cnd)),
              error = function(cnd) paste("ERROR:",
                                          conditionMessage(cnd)))
print(w)
cat("DONE\n")
