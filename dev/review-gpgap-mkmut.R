# Review of lane gpgap, item 7: one mutant build, switched at run time by
# REV_MUT, in C:/Users/adf44/source/r/rev-gpgap-src. Writes nothing in the
# worktree but this script. Every replacement is asserted to have matched.
src <- "C:/Users/adf44/source/r/rev-gpgap-src"
rep1 <- function(txt, old, new) {
  stopifnot(length(grep(old, txt, fixed = TRUE)) >= 1L)
  sub(old, new, txt, fixed = TRUE)
}
f <- file.path(src, "R/covstruct.R")
x <- paste(readLines(f), collapse = "\n")
x <- rep1(x, "gp_nugget_message <- function(fit) {\n  fr <- fit$frame",
          paste0("gp_nugget_message <- function(fit) {\n",
                 "  mut <- Sys.getenv(\"REV_MUT\")\n",
                 "  if (mut == \"nomsg\") return(NULL)\n",
                 "  fr <- fit$frame"))
x <- rep1(x, "if (!identical(lp[[\"dpar\"]], \"mu\")) next",
          "if (mut != \"no_dpar\" && !identical(lp[[\"dpar\"]], \"mu\")) next")
x <- rep1(x, "exact <- Filter(function(g) identical(g[[\"type\"]], \"exact\"), lp[[\"gps\"]])",
          paste0("exact <- Filter(function(g) mut == \"all_types\" || ",
                 "identical(g[[\"type\"]], \"exact\"), lp[[\"gps\"]])"))
x <- rep1(x, "if (!is.null(av[[\"se\"]]) && !isTRUE(av[[\"se_sigma\"]])) next",
          paste0("if (mut != \"no_se\" && !is.null(av[[\"se\"]]) && ",
                 "!isTRUE(av[[\"se_sigma\"]])) next"))
x <- rep1(x, "if (is.null(sg) || !is.null(sg[[\"Z\"]])) next",
          "if (is.null(sg) || (mut != \"no_z\" && !is.null(sg[[\"Z\"]]))) next")
x <- rep1(x, "min(sg[[\"link\"]]$linkinv(eta))",
          "(if (mut == \"sig_max\") max else min)(sg[[\"link\"]]$linkinv(eta))")
x <- rep1(x, "noise <- sqrt(gp_nugget) * sdgp",
          paste0("noise <- (if (mut == \"noise_wrong\") gp_nugget else ",
                 "sqrt(gp_nugget)) * sdgp"))
x <- rep1(x, "if (isTRUE(noise^2 >= 0.01 * sigma^2)) {",
          paste0("thr <- switch(mut, thr_1e3 = 1e-3, thr_1e4 = 1e-4, ",
                 "thr_05 = 0.05, thr_2 = 0.2, 0.01)\n",
                 "      if (isTRUE(noise^2 >= thr * sigma^2)) {"))
# the kernel, on the tape
x <- rep1(x, "exp(-Q) + diag(gp_nugget, nrow(blk[[\"aux_D2\"]][[1]]))",
          paste0("exp(-switch(Sys.getenv(\"REV_MUT\"), kern_102 = 1.02, 1) * Q) + ",
                 "diag(switch(Sys.getenv(\"REV_MUT\"), nug_2e6 = 2e-6, ",
                 "nug_1e7 = 1e-7, gp_nugget), nrow(blk[[\"aux_D2\"]][[1]]))"))
writeLines(x, f)
g <- file.path(src, "R/fit.R")
y <- paste(readLines(g), collapse = "\n")
y <- rep1(y, "frm_warning(gp_msg, call. = FALSE, class = \"frmtmb_gp_nugget\")",
          paste0("frm_warning(gp_msg, call. = FALSE, class = if ",
                 "(Sys.getenv(\"REV_MUT\") != \"class_drop\") ",
                 "\"frmtmb_gp_nugget\")"))
writeLines(y, g)
h <- file.path(src, "tests/testthat/helper-brms.R")
z <- paste(readLines(h), collapse = "\n")
z <- rep1(z, "order(do.call(order, as.data.frame(as.matrix(xb))))",
          paste0("if (Sys.getenv(\"REV_MUT\") == \"ord_none\") return(seq_len(NROW(xb)))\n",
                 "  if (Sys.getenv(\"REV_MUT\") == \"ord_direct\") return(do.call(order, ",
                 "as.data.frame(as.matrix(xb))))\n",
                 "  if (Sys.getenv(\"REV_MUT\") == \"ord_col1\") return(order(order(as.matrix(xb)[, 1L])))\n",
                 "  order(do.call(order, as.data.frame(as.matrix(xb))))"))
z <- rep1(z, "if (is.null(sdat[[sub(\"^zgp_\", \"Jgp_\", zname)]])) {",
          paste0("if (Sys.getenv(\"REV_MUT\") != \"ord_always\" && ",
                 "is.null(sdat[[sub(\"^zgp_\", \"Jgp_\", zname)]])) {"))
writeLines(z, h)
cat("patched\n")
