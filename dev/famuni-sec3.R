setwd("C:/Users/adf44/source/r/frmtmb-wt-famuni")
f <- "dev/famuni-findings.md"
x <- paste(readLines(f), collapse = "\n")
old <- "   holds `xi > -1`. Measured agreement with brms's scalar-`xi` fit is in
   section 6 (z = 0.101 for xi)."
new <- "   holds `xi > -1`. Measured agreement with brms's scalar-`xi` fit is in
   section 6 (z = 0.101 for xi). That argument is for the maximum only.
   For the SAMPLER the support's edge is a wall in the sampled space,
   and frm_sample() diverged on it (review fix 4: 1187, 263 and 3856 of
   4000 at n = 40). Fix round 1 samples a scalar `xi` through brms's
   `tmp_xi` instead (section 12)."
stopifnot(grepl(old, x, fixed = TRUE))
x <- sub(old, new, x, fixed = TRUE)
writeLines(strsplit(x, "\n")[[1]], f)
