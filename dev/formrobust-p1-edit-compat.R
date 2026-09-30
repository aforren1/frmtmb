# Punch round 1, m3: the compat fitted row names the divergence. Record.
p <- "C:/Users/adf44/source/r/frmtmb-wt-formrobust/R/compat.R"
x <- paste(readLines(p), collapse = "\n")
old <- "a response that is NA, or absent, is filled with its expected value."
stopifnot(lengths(regmatches(x, gregexpr(old, x, fixed = TRUE))) == 1L)
x <- sub(old, paste0("a response that is NA, or absent, is filled with its ",
  "expected value. That diverges from brms on purpose: brms fills it with ",
  "a draw per posterior draw, so its fitted() carries the fill's spread, ",
  "and a maximum likelihood fit has no draws to propagate the fill ",
  "through. frmtmb.sample's posterior_epred() fills with draws, as brms ",
  "does."), x, fixed = TRUE)
con <- file(p, "wb"); writeLines(strsplit(x, "\n")[[1]], con, sep = "\r\n")
close(con)
