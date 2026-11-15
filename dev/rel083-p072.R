# The brief's survey item 1: pair p072 (asym_laplace x gr(cov = A)) is
# recorded as a decided refusal. Its adjudication row of 0.80.0 read the
# L1 identity of a model frmtmb then fitted; under the user's decision of
# 2026-11-13 (lane wcheck2) frm() refuses asym_laplace() with a random
# effect in mu under Laplace, so the measured verdict is now "defect,
# frmtmb refuses", and the decision makes it a divergence by decision.
f <- "dev/equiv-data/adjudication.tsv"
x <- readLines(f)
i <- grep("^p072\t", x)
stopifnot(length(i) == 1)
x[i] <- paste("p072", "divergence by decision", "FALSE", "", "",
  paste("asym_laplace x gr(cov = A): frm() refuses asym_laplace() with a",
        "random effect in mu under the Laplace approximation, by the user's",
        "decision of 2026-11-13 (NEWS 0.83.0, \"asym_laplace() with random",
        "effects\"): the refusal proves that the Laplace objective has no",
        "maximum there, and frm_sample() fits the model. brms builds it.",
        "Before 0.83.0 the L1 identity held at three points (constants",
        "-1.7e-08, -1.7e-05 and 0)."), sep = "\t")
con <- file(f, "wb"); writeLines(x, con, sep = "\r\n"); close(con)
cat("p072 adjudicated\n")
