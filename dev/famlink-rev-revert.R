## Reviewer check for lane wt-famlink, priority 5: does each new test
## FAIL when only the fix it guards is reverted, in the lane build?
##
## Each revert patches ONE piece of the lane's namespace back to the
## base behavior (or to "no check"), then runs
## tests/testthat/test-brms-families.R and lists the test blocks that
## fail or error. A guard whose block still passes under its own revert
## fails open.
##
## Usage: Rscript dev/famlink-rev-revert.R <revert>
##   none | no_rederive | partial_dollar | no_mu_set | no_mix_check |
##   no_bern_msg | ord_warn | no_unquoted | no_brms_dpar_links |
##   exact_names_only
rev <- commandArgs(trailingOnly = TRUE)[1]
ARM <- "lane"
source("dev/famlink-rev-common.R")
suppressMessages(library(testthat))
ns <- asNamespace("frmtmb")
put <- function(name, value) {
  unlockBinding(name, ns); assign(name, value, envir = ns); lockBinding(name, ns)
}
switch(rev,
  none = NULL,
  no_rederive = put("family_link_sources", character()),
  partial_dollar = { pf <- function(x, name) unclass(x)[[name, exact = FALSE]]; put("$.frmtmb_family", pf); registerS3method("$", "frmtmb_family",
    pf, envir = ns) },
  no_mu_set = put("mu_link", function(value, family, choices = NULL) {
    if (is.list(value)) return(frmtmb:::get_link(value, dpar = "mu", family = family))
    frmtmb:::get_link(value, dpar = "mu")
  }),
  no_mix_check = put("mixture_check_components", function(comps) invisible(NULL)),
  no_bern_msg = put("suggest_bernoulli", function(spec, frame) invisible(NULL)),
  ord_warn = {
    src <- deparse(get("extract_y", ns))
    i <- grep("requires either \", \"positive integers", src, fixed = TRUE)
    if (length(i) != 1L) i <- grep("requires either", src)[1]
    j <- max(grep("stop(\"Family '\"", src[seq_len(i)], fixed = TRUE))
    stopifnot(length(j) == 1L)
    src[j] <- sub("stop(", "warning(", src[j], fixed = TRUE)
    f <- eval(parse(text = src), envir = ns)
    environment(f) <- ns
    put("extract_y", f)
  },
  no_unquoted = put("link_arg_value", function(slink, link, choices, default) {
    if (is.null(link) || identical(link, NA)) default else link
  }),
  no_brms_dpar_links = {
    src <- deparse(get("as_frmtmb_family", ns))
    k <- grep("extra <- intersect", src, fixed = TRUE)
    stopifnot(length(k) == 1L)
    src[k] <- "extra <- character(0); if (FALSE) intersect("
    f <- eval(parse(text = src), envir = ns); environment(f) <- ns
    put("as_frmtmb_family", f)
  },
  exact_names_only = put("family_ctor", function(name) {
    ctor <- frmtmb:::family_registry[[name]]
    if (is.null(ctor)) stop(name, " is not a supported family", call. = FALSE)
    ctor
  }),
  stop("unknown revert ", rev)
)
res <- testthat::test_file("tests/testthat/test-brms-families.R",
                           reporter = "silent", package = "frmtmb",
                           stop_on_failure = FALSE)
d <- as.data.frame(res)
bad <- d[d$failed > 0 | d$error, ]
cat(sprintf("REVERT %-20s blocks %d failing %d | pass %d fail %d error %d\n",
            rev, nrow(d), nrow(bad), sum(d$passed), sum(d$failed), sum(d$error)))
for (t in bad$test) cat("   FAILS:", t, "\n")
