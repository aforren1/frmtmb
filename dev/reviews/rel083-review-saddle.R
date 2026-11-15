# Release review of 0.83.0: under frmtmb_control(nl_starts = FALSE) the
# saddle restart still runs and still says so (frmtmb_nl_restart), so the
# fit is not "the one the optimizer reaches from the start". Case of
# test-se-limit.R, "a family limit on a nonlinear fit still triggers the
# restart".
.libPaths(c("C:/Users/adf44/source/r/rellib-r22",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
Sys.setenv(R_MAKEVARS_USER = "C:/Users/adf44/Documents/.R/Makevars.win")
suppressPackageStartupMessages(library(frmtmb))
cat("lib:", dirname(find.package("frmtmb")), as.character(packageVersion("frmtmb")), "\n")
n <- 120L
set.seed(81100 + 1000 * 4 + 100 * 2 + n + 10 * 6)
x <- sort(round(stats::runif(n, 0, 10), 3))
z <- round(stats::rnorm(n), 3)
w <- round(stats::runif(n, 0.2, 5), 3)
f <- 2 * sin(3 * x) + x
fs <- f / max(abs(f))
y <- (2 + fs) * w^0.5 + 0.1 * stats::rt(n, 3)
d <- data.frame(x, z, w, y)
for (ns in c(TRUE, FALSE)) {
  cl <- character()
  r <- withCallingHandlers(
    frm(bf(y ~ a * w^b, a ~ gp(x), b ~ 1, nl = TRUE), data = d,
        family = student(), control = frmtmb_control(nl_starts = ns)),
    warning = function(w) { cl <<- c(cl, paste0("W:", class(w)[1])); invokeRestart("muffleWarning") },
    message = function(m) { cl <<- c(cl, paste0("M:", class(m)[1])); invokeRestart("muffleMessage") })
  cat("nl_starts", ns, "logLik", format(as.numeric(logLik(r)), digits = 10),
      "conds", paste(cl, collapse = " "), "\n")
}
