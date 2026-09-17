# Reviewer, lane wt-priorform: the duplicate-prior tracer of
# dev/priorform-rev-priortrace-run1.R over every Rd example of core and
# the extensions (frmtmb.sample with its donttest sections commented).
#   Rscript dev/priorform-rev-priortrace-rd.R
.libPaths(c("C:/Users/adf44/source/r/priorform-lib",
            "C:/Users/adf44/source/r/rellib-r3",
            "C:/Users/adf44/source/r/pinlib",
            "C:/Users/adf44/AppData/Local/R/win-library/4.6"))
root <- "C:/Users/adf44/source/r/frmtmb-wt-priorform"
Sys.setenv(FRMTMB_STAN_CACHE = file.path(root, "dev", "stan-cache"))
src <- readLines(file.path(root, "dev", "priorform-rev-priortrace-run1.R"))
# reuse the tracer definition verbatim: everything from the logf line to
# the trace() call
a <- grep("^logf <- ", src); b <- grep("^suppressMessages[(]trace", src)
suppressMessages(library(testthat)); suppressMessages(library(frmtmb))
ns <- asNamespace("frmtmb")
pkgdir <- "Rd"; filt <- ""
eval(parse(text = src[a:(b + 1L)]))
pkgs <- c(frmtmb = ".", frmtmb.coupling = "extensions/frmtmb.coupling",
          frmtmb.eam = "extensions/frmtmb.eam", frmtmb.latent = "extensions/frmtmb.latent",
          frmtmb.learn = "extensions/frmtmb.learn", frmtmb.ode = "extensions/frmtmb.ode",
          frmtmb.spline = "extensions/frmtmb.spline", frmtmb.sample = "extensions/frmtmb.sample")
n_ex <- 0L
for (p in names(pkgs)) {
  suppressMessages(library(p, character.only = TRUE))
  for (rd in list.files(file.path(root, pkgs[[p]], "man"), "[.]Rd$", full.names = TRUE)) {
    ex <- tempfile(fileext = ".R")
    tools::Rd2ex(rd, ex, commentDontrun = TRUE,
                 commentDonttest = identical(p, "frmtmb.sample"))
    if (!file.exists(ex)) next
    assign(".rev_label", paste0("Rd ", p, "::", basename(rd)), envir = globalenv())
    set.seed(1); n_ex <- n_ex + 1L
    tryCatch(suppressWarnings(suppressMessages(utils::capture.output(
      sys.source(ex, envir = new.env(parent = globalenv()))))),
      error = function(e) cat("ERR", basename(rd), substr(conditionMessage(e), 1, 80), "\n"))
  }
}
cat("examples run:", n_ex, "\n")
