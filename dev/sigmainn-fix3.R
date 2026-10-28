# helper-brms.R (core and frmtmb.spatial): brms's sigma is frmtmb's
# sigma under ar(cov = TRUE) now, so the translator's scale rule goes.
for (f in c("tests/testthat/helper-brms.R",
            "extensions/frmtmb.spatial/tests/testthat/helper-brms.R")) {
  x <- readLines(f)
  i <- grep("brms_sigma_scale(fit, code, nm)", x, fixed = TRUE)
  stopifnot(length(i) == 1L,
            grepl("link$linkinv(fe[[\"(Intercept)\"]])) *", x[i - 1L],
                  fixed = TRUE))
  x[i - 1L] <- sub(" *", "", x[i - 1L], fixed = TRUE)
  x <- x[-i]
  s <- grep("^# brms's sigma is not always frmtmb's sigma[.]$", x)
  e <- grep("^brms_sigma_scale <- function", x)
  stopifnot(length(s) == 1L, length(e) == 1L)
  close <- e + which(x[(e + 1L):length(x)] == "}")[1L]
  # drop the comment block, the function and the blank line after it
  x <- x[-(s:(close + 1L))]
  writeLines(x, f)
}
