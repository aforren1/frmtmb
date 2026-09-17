# Lane wt-conditions: write tests/testthat/test-conditions-census.R into
# frmtmb and each extension from ONE template, so the eight copies
# cannot drift. Each copy differs only in its package name, its anchor
# file and its exemption list.
#
#   Rscript dev/conditions-census-gen.R
#
# Run from the worktree root. Edit this file and
# dev/conditions-census-template.R, not the generated tests.

exempt <- list(
  frmtmb = c(
    # the helpers themselves raise through the base functions
    'R/conditions.R | frm_stop | stop(cnd)',
    'R/conditions.R | frm_warning | warning(cnd)',
    'R/conditions.R | frm_message | message(cnd)',
    # rethrows of a caught condition, which keeps its own class: the
    # optimizer's error after every restart has failed, the same error
    # when there is no better point to restart from, and a
    # frmtmb_fit_error from the inner autoscale fit, which is already
    # classed and must not be wrapped twice
    'R/fit.R | fit_assembled | stop(e)',
    'R/fit.R | optimizer_from_best | stop(e)',
    'R/fit.R | fit_error_context | stop(e)'),
  frmtmb.coupling = character(0),
  frmtmb.eam = character(0),
  frmtmb.latent = character(0),
  frmtmb.learn = character(0),
  frmtmb.ode = character(0),
  frmtmb.sample = character(0),
  frmtmb.spline = character(0))

anchor <- c(frmtmb = "objective.R",
            frmtmb.coupling = "frmtmb.coupling-package.R",
            frmtmb.eam = "frmtmb.eam-package.R",
            frmtmb.latent = "frmtmb.latent-package.R",
            frmtmb.learn = "frmtmb.learn-package.R",
            frmtmb.ode = "frmtmb.ode-package.R",
            frmtmb.sample = "frmtmb.sample-package.R",
            frmtmb.spline = "frmtmb.spline-package.R")

# The template is a file of its own so that it is plain R, with no
# second level of string escaping to get wrong.
template <- paste(readLines("dev/conditions-census-template.R"),
                  collapse = "\n")

for (pkg in names(exempt)) {
  dir <- if (pkg == "frmtmb") "." else file.path("extensions", pkg)
  ex <- exempt[[pkg]]
  ex_txt <- if (length(ex)) {
    paste0("  \"", gsub("\"", "\\\\\"", ex), "\"", collapse = ",\n")
  } else "  character(0)"
  txt <- gsub("@PKG@", pkg, template, fixed = TRUE)
  txt <- gsub("@ANCHOR@", anchor[[pkg]], txt, fixed = TRUE)
  txt <- sub("@EXEMPT@", ex_txt, txt, fixed = TRUE)
  out <- file.path(dir, "tests", "testthat", "test-conditions-census.R")
  con <- file(out, open = "wb")
  writeLines(txt, con, sep = "\n")
  close(con)
  cat("wrote", out, "with", length(ex), "exemptions\n")
}
