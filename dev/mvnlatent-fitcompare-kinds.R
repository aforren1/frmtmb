# Which quantities of a dev/mvnlatent-fitcompare.R save are values and
# which are refusals (character), so "identical" is read correctly.
x <- readRDS(commandArgs(trailingOnly = TRUE)[1])
q <- c("se", "fixef", "predict", "fitted", "probs", "simulate")
for (k in q) {
  ch <- vapply(x, function(f) is.character(f[[k]]), TRUE)
  cat(k, ": value in", sum(!ch), "fits; refusal in", sum(ch), "\n")
  if (any(ch)) print(unique(vapply(x[ch], function(f) f[[k]], "")))
}
cat("convergence codes:", paste(names(table(vapply(x, `[[`, 0, "conv"))),
                                 table(vapply(x, `[[`, 0, "conv")),
                                 collapse = ", "), "\n")
se <- unlist(lapply(x, `[[`, "se"))
cat("SE entries", length(se), "finite", sum(is.finite(as.numeric(se))), "\n")
