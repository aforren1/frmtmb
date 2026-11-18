# Review of lane ordstarts: summary of dev/reviews/ordstarts-review-nore-*
# (no random effects; starts = 0 against the shipped default of 10).
#   Rscript dev/reviews/ordstarts-review-nore-sum.R
rd <- function(f) utils::read.delim(f, stringsAsFactors = FALSE,
                                    quote = "", comment.char = "")
N <- do.call(rbind, lapply(Sys.glob("dev/reviews/ordstarts-review-nore-*.tsv"),
                           rd))
by <- paste(N$family, N$n)
d <- N$D_ll - N$A0_ll
gain <- d > 1e-3 & !is.na(d)
med <- gain & N$A0_cls == "proper" & !N$A0_deg & N$D_deg
cat("fits", nrow(N), "\n")
print(data.frame(
  seeds = tapply(rep(1, nrow(N)), by, sum),
  up = tapply(gain, by, sum),
  down = tapply(d < -1e-3, by, sum, na.rm = TRUE),
  A0_deg = tapply(N$A0_deg, by, sum, na.rm = TRUE),
  D_deg = tapply(N$D_deg, by, sum, na.rm = TRUE),
  proper_silent_to_deg = tapply(med, by, sum, na.rm = TRUE),
  el_ratio = round(tapply(N$D_el, by, sum) / tapply(N$A0_el, by, sum), 2)))
cat("\nproper default without the warning moved up to a warned point:\n")
print(N[med, c("family", "n", "seed", "A0_ll", "D_ll", "A0_tse", "D_tse",
               "D_cls", "D_nwarn")], row.names = FALSE)
cat("\nevery D with the warning where A0 had none:\n")
print(N[N$D_deg & !N$A0_deg, c("family", "n", "seed", "A0_ll", "D_ll",
                                "A0_cls", "A0_tse", "D_tse")],
      row.names = FALSE)
