source("C:/Users/adf44/source/r/frmtmb-wt-adefects/dev/adefects-rev-prelude.R")
suppressPackageStartupMessages({ library(brms) })

cat("== brms itself carries allow_new_levels through ... and ANSWERS\n")
bf1 <- get("brmsfit_example1", envir = asNamespace("brms"))
nd <- bf1$data[1:3, ]
nd$visit <- NULL
for (fn in c("posterior_epred", "posterior_predict", "log_lik")) {
  g <- get(fn, envir = asNamespace("brms"))
  for (anl in c(TRUE, FALSE)) {
    r <- tryCatch(g(bf1, newdata = nd, allow_new_levels = anl),
                  error = function(e) e)
    cat(sprintf("  %-18s allow_new_levels = %-5s -> %s\n", fn, anl,
                if (inherits(r, "condition"))
                  paste("ERR:", substr(gsub("[\r\n]+", " ",
                                            conditionMessage(r)), 1, 66))
                else paste("OK", paste(dim(r), collapse = "x"))))
  }
}
cat("  brms formals have allow_new_levels:",
    "allow_new_levels" %in% names(formals(
      get("posterior_epred.brmsfit", envir = asNamespace("brms")))), "\n")

cat("\n== the positional guard is sensitive to the WRONG placement\n")
pos <- function(f) {
  a <- names(formals(f))
  a[seq_len(match("...", a, nomatch = length(a) + 1L) - 1L)]
}
fdiv <- function(b, o) {
  k <- min(length(b), length(o)); if (!k) return(NA_integer_)
  d <- which(b[seq_len(k)] != o[seq_len(k)]); if (length(d)) d[[1L]] else NA
}
b <- pos(get("posterior_epred.brmsfit", envir = asNamespace("brms")))
cat("  brms positional slots:", paste(b, collapse = ", "), "\n")
suppressPackageStartupMessages(library(frmtmb.sample))
o_now <- pos(get("posterior_epred.frmtmb_draws",
                 envir = asNamespace("frmtmb.sample")))
cat("  lane positional slots:", paste(o_now, collapse = ", "), "\n")
cat("  divergence as shipped (want NA)   :", fdiv(b, o_now), "\n")
wrong <- function(object, newdata = NULL, re_formula = NULL,
                  re.form = NULL, resp = NULL, dpar = NULL, nlpar = NULL,
                  ndraws = NULL, draw_ids = NULL,
                  allow_new_levels = FALSE, ...) NULL
cat("  divergence if placed BEFORE dots  :", fdiv(b, pos(wrong)),
    "(want 10)\n")

cat("\n== the lane methods: named works, positional does not steal a slot\n")
ds <- structure(list(), class = "frmtmb_draws")
probe <- function(expr) {
  r <- tryCatch(eval(expr), error = function(e) conditionMessage(e))
  substr(gsub("[\r\n]+", " ", if (is.character(r)) r else "NOERROR"), 1, 95)
}
cat("  epred, allow_new_levels named   :",
    probe(quote(posterior_epred(ds, newdata = data.frame(x = 1),
                                allow_new_levels = TRUE))), "\n")
cat("  epred, brms's 10th positional   :",
    probe(quote(posterior_epred(ds, NULL, NULL, NULL, NULL, NULL, NULL,
                                NULL, NULL, TRUE))), "\n")
cat("  epred, a brms dots name brms takes:",
    probe(quote(posterior_epred(ds, sample_new_levels = "gaussian"))), "\n")
cat("  epred, a misspelling            :",
    probe(quote(posterior_epred(ds, allow_new_level = TRUE))), "\n")
cat("  ppred, allow_new_levels named   :",
    probe(quote(posterior_predict(ds, newdata = data.frame(x = 1),
                                  allow_new_levels = TRUE))), "\n")
cat("  log_lik, allow_new_levels named :",
    probe(quote(log_lik(ds, allow_new_levels = TRUE))), "\n")
cat("  epred, allow_new_levels = 'yes' :",
    probe(quote(posterior_epred(ds, allow_new_levels = "yes"))), "\n")
