# Which spellings reach mixture_mvn() on a build: the bare string, the
# brmsfamily() string, and the call with and without frmtmb.latent.
#
#   Rscript dev/mvnlatent-probe-strings.R <lib or "base">
a <- commandArgs(trailingOnly = TRUE)
base <- c("C:/Users/adf44/source/r/rellib-r11",
          "C:/Users/adf44/AppData/Local/R/win-library/4.6")
.libPaths(if (identical(a[1], "base")) base else c(a[1], base))
suppressMessages(library(frmtmb))
cat("frmtmb", as.character(packageVersion("frmtmb")),
    dirname(find.package("frmtmb")), "\n")
set.seed(1)
Y <- rbind(matrix(rnorm(60, 0), ncol = 2), matrix(rnorm(60, 4), ncol = 2))
dd <- data.frame(row = seq_len(nrow(Y)))
dd$Y <- Y
try_msg <- function(label, expr) {
  r <- tryCatch({
    force(expr)
    "OK"
  }, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(label, ":", r, "\n")
}
arm <- function(tag) {
  cat("--", tag, "--\n")
  try_msg("string", frm(bf(Y ~ 1), family = "mixture_mvn", data = dd))
  try_msg("brmsfamily", brmsfamily("mixture_mvn"))
  try_msg("brmsfamily K D", brmsfamily("mixture_mvn", K = 2, D = 2))
  try_msg("call", frm(bf(Y ~ 1) + mixture_mvn(K = 2, D = 2), data = dd))
  try_msg("exists", stopifnot(exists("mixture_mvn")))
  cat("frmtmb exports mixture_mvn:",
      "mixture_mvn" %in% getNamespaceExports("frmtmb"), "\n")
}
arm("frmtmb only")
suppressMessages(library(frmtmb.latent))
cat("frmtmb.latent", as.character(packageVersion("frmtmb.latent")),
    dirname(find.package("frmtmb.latent")), "\n")
arm("with frmtmb.latent")
cat("search():", paste(search(), collapse = " "), "\n")
