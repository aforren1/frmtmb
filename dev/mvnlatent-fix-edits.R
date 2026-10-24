# Text edits of fix round 1 that are easier as exact replacements; each
# old string must occur exactly once.
rep_in <- function(f, old, new) {
  x <- gsub("\r\n", "\n", readChar(f, file.info(f)$size, useBytes = TRUE),
            fixed = TRUE)
  n <- lengths(regmatches(x, gregexpr(old, x, fixed = TRUE)))
  if (n != 1L) stop(f, ": ", n, " matches for: ", old)
  writeChar(sub(old, new, x, fixed = TRUE), f, eos = NULL, useBytes = TRUE)
}
rep_in("tests/testthat/test-core-boundary.R",
"  # The family moved to frmtmb.latent, whose own test-mixture-mvn.R is
  # the arm with that package loaded. No string spelling reached the
  # family before the move (brms has no such family), so a string is
  # the generic refusal, and the call is an unknown function.
  expect_false(\"mixture_mvn\" %in% getNamespaceExports(\"frmtmb\"))
  skip_if(isNamespaceLoaded(\"frmtmb.latent\"))
  d <- data.frame(row = 1:4)",
"  # The family moved to frmtmb.latent, whose own test-mixture-mvn.R is
  # the arm with that package loaded. No string spelling reached the
  # family before the move (brms has no such family), so a string is
  # the generic refusal whether or not that package is loaded, and the
  # call is an unknown function unless it is attached.
  expect_false(\"mixture_mvn\" %in% getNamespaceExports(\"frmtmb\"))
  d <- data.frame(row = 1:4)")
rep_in("tests/testthat/test-core-boundary.R",
"  expect_error(brmsfamily(\"mixture_mvn\"), \"is not a supported family\")
  expect_false(exists(\"mixture_mvn\"))",
"  expect_error(brmsfamily(\"mixture_mvn\"), \"is not a supported family\")
  skip_if(\"package:frmtmb.latent\" %in% search())
  expect_false(exists(\"mixture_mvn\"))")
cat("done\n")
