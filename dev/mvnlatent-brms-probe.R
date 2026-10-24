# What brms 2.23.0 offers for a multivariate normal mixture, so the
# findings can say which brms spelling, if any, reaches mixture_mvn().
# No Stan compile: stancode() only.
.libPaths(c("C:/Users/adf44/AppData/Local/R/win-library/4.6"))
suppressMessages(library(brms))
cat("brms", as.character(packageVersion("brms")), "\n")
ex <- getNamespaceExports("brms")
cat("exports matching mvn|mvnorm|multi_normal|mixture:",
    paste(sort(grep("mvn|mvnorm|multi_normal|mixture", ex, value = TRUE)),
          collapse = ", "), "\n")
fams <- tryCatch(brms:::.family_names(), error = function(e) NULL)
if (is.null(fams)) {
  fams <- tryCatch(ls(asNamespace("brms"), pattern = "^\\.family_"),
                   error = function(e) character(0))
}
cat("family names matching mvn|multi:",
    paste(grep("mvn|multi|mvnorm", fams, value = TRUE), collapse = ", "),
    "\n")
try_msg <- function(label, expr) {
  r <- tryCatch({
    force(expr)
    "OK"
  }, error = function(e) paste("ERROR:", conditionMessage(e)))
  cat(label, ":", r, "\n")
}
set.seed(1)
d <- data.frame(y1 = rnorm(40), y2 = rnorm(40))
try_msg("brmsfamily('mixture_mvn')", brmsfamily("mixture_mvn"))
try_msg("brmsfamily('mvnormal')", brmsfamily("mvnormal"))
# a per-response mixture in a multivariate model, with and without
# residual correlation: the nearest brms has to a joint MVN mixture
try_msg("mvbind + mixture, rescor = FALSE",
        stancode(bf(mvbind(y1, y2) ~ 1) + set_rescor(FALSE),
                 family = mixture(gaussian, gaussian), data = d))
try_msg("mvbind + mixture, rescor = TRUE",
        stancode(bf(mvbind(y1, y2) ~ 1) + set_rescor(TRUE),
                 family = mixture(gaussian, gaussian), data = d))
