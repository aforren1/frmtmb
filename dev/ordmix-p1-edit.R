# Punch round 1: two compat.R row texts (m6, the cs() x cumulative reason).
p <- "C:/Users/adf44/source/r/frmtmb-wt-ordmix/R/compat.R"
s <- readLines(p)
old1 <- "to at most 2.5 ulp over 20 shapes (dev/ordmix-lpcheck.R);"
new1 <- paste0("to at most 2.5 ulp over the 20 shapes of dev/ordmix-lpcheck.R. ",
  "A probit component whose latent distance from a threshold passes ",
  "about 38 has a NaN density where brms's is finite (the review's ",
  "three-component probit, sratio and acat mixture, one perturbed ",
  "point; dev/test-backlog.md);")
old2 <- "Refused: category-specific effects are not identified under the cumulative parameterization."
new2 <- paste0("Refused: cs() moves each threshold of a row by its own ",
  "amount, so under the cumulative parameterization a row's thresholds ",
  "can cross, and the category between two crossed thresholds then has ",
  "a negative probability. brms 2.23.0 fits it, with a warning that the ",
  "effects are experimental; hurdle_cumulative() takes it as brms does.")
i1 <- grep(old1, s, fixed = TRUE); i2 <- grep(old2, s, fixed = TRUE)
stopifnot(length(i1) == 1, length(i2) == 1)
s[i1] <- sub(old1, new1, s[i1], fixed = TRUE)
s[i2] <- sub(old2, new2, s[i2], fixed = TRUE)
writeLines(s, p)
cat("done\n")
