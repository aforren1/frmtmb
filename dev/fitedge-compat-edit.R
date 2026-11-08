# One-off: the two frm_compat() notes that state the mi() refusal beside a
# latent residual, changed to the refusal's new scope (item 4). Notes only;
# each row keeps its status.
f <- "C:/Users/adf44/source/r/frmtmb-wt-fitedge/R/compat.R"
x <- readLines(f, warn = FALSE)
a1 <- "and mi() on the response is refused because no measurement reads a latent response beside a latent residual."
b1 <- "and mi() on the response is refused where another density reads its missing values, because no measurement reads a latent response beside a latent residual; a missing value no density reads leaves the objective, and the fit is then exact on lognormal."
a2 <- "mi() on the response is refused, since no measurement reads a latent response beside a latent residual."
b2 <- "mi() on the response is refused where another density reads its missing values, since no measurement reads a latent response beside a latent residual; a missing value no density reads leaves the objective."
stopifnot(sum(grepl(a1, x, fixed = TRUE)) == 1L, sum(grepl(a2, x, fixed = TRUE)) == 1L)
x <- sub(a1, b1, x, fixed = TRUE)
x <- sub(a2, b2, x, fixed = TRUE)
writeLines(x, f)
cat("edited\n")
