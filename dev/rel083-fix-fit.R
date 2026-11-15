# refits x gpfield, the end of frm() in R/fit.R (a conflict beyond the
# brief's list): refits made the end `fit <- if (...) run(start, TRUE)
# else nl_fit_probed(...)` so that it can store refit_args on the fit;
# gpfield wraps nl_fit_probed() in gp_end_warn(), which says the gp()
# warnings of a refit the probe returns. The merge keeps refits' shape
# and puts gpfield's wrapper on the nonlinear branch, where it was.
f <- "R/fit.R"
x <- readLines(f)
a <- grep("^<<<<<<< rel$", x); m <- grep("^=======$", x)
b <- grep("^>>>>>>> gpfield$", x)
stopifnot(length(a) == 1, length(m) == 1, length(b) == 1)
stopifnot(identical(x[(a + 1):(m - 1)], "  fit"))
x <- c(x[seq_len(a - 1)], "  fit", x[-seq_len(b)])
i <- which(x == "    nl_fit_probed(run, frame, start, prior)")
stopifnot(length(i) == 1)
x[i] <- "    gp_end_warn(nl_fit_probed(run, frame, start, prior))"
x <- append(x, c(
  "    # the fit returned can be a refit (a sign probe's, a saddle",
  "    # restart's), which runs unannounced: its gp() warnings are said here"),
  after = i - 1)
con <- file(f, "wb"); writeLines(x, con, sep = "\r\n"); close(con)
cat("resolved\n")
