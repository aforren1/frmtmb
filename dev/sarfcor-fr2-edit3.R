# Fix round 2: predict() draws on the block over every row of the
# term's matrix; the brms row's expected NA warnings. Part A:
# Rscript <this> <root>
root <- commandArgs(TRUE)[1]
edit <- function(p, pairs) {
  p <- file.path(root, p)
  txt <- paste(readLines(p), collapse = "\n")
  for (pr in pairs) {
    n <- lengths(regmatches(txt, gregexpr(pr[1], txt, fixed = TRUE)))
    if (n != 1L) stop(p, " matches: ", n, " for ", substr(pr[1], 1, 60))
    txt <- sub(pr[1], pr[2], txt, fixed = TRUE)
  }
  writeLines(strsplit(txt, "\n", fixed = TRUE)[[1]], p)
}
edit("R/predict-brms.R", list(c(
  "      autocor_seam_for_newdata(object, ctx0[[\"autocor\"]], rspec, newdata,
                               \"predict()\")
      ctx0[[\"autocor\"]] <- NULL",
  "      # the replicates draw through the block over newdata's rows,
      # which the fit's own frame then carries into every draw
      object$frame[[\"autocor\"]][[rspec$resp_name]] <-
        autocor_seam_for_newdata(object, ctx0[[\"autocor\"]], rspec, newdata,
                                 \"predict()\")
      ctx0[[\"autocor\"]] <- NULL")))
edit("extensions/frmtmb.spatial/tests/testthat/test-sar.R", list(
  c("  prior <- brms_flat_prior(bform, data = dd, family = fam,
                           data2 = list(Ws = Ws))
  code <- brms::make_stancode(bform, data = dd, family = fam,
                              prior = prior, data2 = list(Ws = Ws))
  sdat <- brms_standata(bform, data = dd, family = fam, prior = prior,
                        data2 = list(Ws = Ws))
  bd <- suppressMessages(brms::brm(bform, data = dd, family = fam,
                                   prior = prior, data2 = list(Ws = Ws),
                                   empty = TRUE))",
    "  # brms says so each time it drops the two rows
  na_note <- \"Rows containing NAs were excluded from the model\"
  prior <- allow_warnings(brms_flat_prior(bform, data = dd, family = fam,
                                          data2 = list(Ws = Ws)), na_note)
  code <- allow_warnings(brms::make_stancode(bform, data = dd,
                                             family = fam, prior = prior,
                                             data2 = list(Ws = Ws)),
                         na_note)
  sdat <- allow_warnings(brms_standata(bform, data = dd, family = fam,
                                       prior = prior,
                                       data2 = list(Ws = Ws)), na_note)
  bd <- allow_warnings(suppressMessages(
    brms::brm(bform, data = dd, family = fam, prior = prior,
              data2 = list(Ws = Ws), empty = TRUE)), na_note)")))
cat("edited", root, "\n")
