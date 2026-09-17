## Reviewer check for lane wt-famlink, priority 3: the mixture support
## check reads frmtmb's `type` as brms's support. Print both per family,
## and independently recompute brms's accepted mean links by CALLING
## brms::brmsfamily(name, link) on every link name frmtmb and brms know,
## rather than reading .family_<name>()$links as the generator does.
ARM <- "lane"
source("dev/famlink-rev-common.R")
suppressMessages(library(brms))
reg <- frmtmb:::family_registry
src <- frmtmb:::brms_mu_link_source
cat("---- type: frmtmb vs brms ----\n")
for (nm in names(src)) {
  f <- tryCatch(reg[[nm]](), error = function(e) NULL)
  if (is.null(f)) f <- tryCatch(get(nm)(), error = function(e) NULL)
  ft <- if (is.null(f)) "<not built>" else (f[["type"]] %||% "<NULL>")
  bt <- tryCatch(get(paste0(".family_", src[[nm]]), asNamespace("brms"))()$type,
                 error = function(e) "?")
  cat(sprintf("%-28s frmtmb %-12s brms %s\n", nm, ft, paste(bt, collapse = ",")))
}
cat("---- independent mean-link sets ----\n")
roster <- unique(c(names(frmtmb:::frmtmb_links), "1/mu^2", "inverse", "sqrt",
                   "logm1", "tan_half", "softit", "probit_approx", "squareplus",
                   "softplus", "cauchit", "cloglog", "log1p", "logc", "inv",
                   "1/mu", "identity", "log", "logit", "probit"))
n_diff <- 0L
for (nm in names(src)) {
  if (identical(nm, "multinomial")) next
  ok <- roster[vapply(roster, function(lk) {
    !inherits(tryCatch(brms::brmsfamily(src[[nm]], link = lk),
                       error = function(e) e), "error")
  }, TRUE)]
  gen <- frmtmb:::brms_mu_links[[nm]]
  in_roster_gen <- intersect(gen, roster)
  same <- setequal(ok, in_roster_gen)
  first_ok <- identical(gen[1], get(paste0(".family_", src[[nm]]),
                                    asNamespace("brms"))()$links[1])
  if (!same || !first_ok) {
    n_diff <- n_diff + 1L
    cat(nm, ": brms by call:", paste(sort(ok), collapse = " "),
        "| generated:", paste(sort(in_roster_gen), collapse = " "), "\n")
  }
  notroster <- setdiff(gen, roster)
  if (length(notroster)) cat(nm, ": generated links outside the probe roster:",
                             notroster, "\n")
}
cat("families whose set differs:", n_diff, "of", length(src) - 1L, "\n")
cat("links in generated sets that frmtmb's roster lacks:",
    paste(setdiff(unlist(frmtmb:::brms_mu_links), names(frmtmb:::frmtmb_links)),
          collapse = " "), "\n")
